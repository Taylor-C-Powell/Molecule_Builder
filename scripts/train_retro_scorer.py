#!/usr/bin/env python3
"""Train ML retrosynthetic disconnection scorer from heuristic-labeled data.

Produces ``molbuilder/data/retro_scorer.pkl`` -- a joblib dict
containing a GradientBoostingRegressor plus metadata.

Usage:
    python scripts/train_retro_scorer.py

Requirements:
    pip install molbuilder[ml]   # scikit-learn >= 1.3, joblib >= 1.3
"""

from __future__ import annotations

import hashlib
import math
import os
import random
import sys
import time

import numpy as np

# Ensure the project root is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import (
    retrosynthesis,
    score_disconnection,
    is_purchasable,
    Precursor,
)
from molbuilder.reactions.functional_group_detect import detect_functional_groups
from molbuilder.reactions.retro_features import (
    extract_retro_features,
    ALL_RETRO_FEATURE_NAMES,
)

# ---------------------------------------------------------------------------
#  Target SMILES for training data generation
# ---------------------------------------------------------------------------
# Drug-like molecules, building blocks, and heterocycles that produce
# diverse disconnection sets.

_TARGET_SMILES: list[str] = [
    # Common drugs (small MW, well-covered by templates)
    "CC(=O)Oc1ccccc1C(=O)O",           # aspirin
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",      # ibuprofen
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O",      # caffeine
    "CC(=O)Nc1ccc(O)cc1",              # acetaminophen
    "OC(=O)c1ccccc1O",                 # salicylic acid
    "c1ccc2c(c1)cc1ccccc12",           # fluorene
    "c1ccc(-c2ccccc2)cc1",             # biphenyl
    "c1ccc(Oc2ccccc2)cc1",             # diphenyl ether
    "CC(=O)c1ccccc1",                  # acetophenone
    "c1ccc(C=O)cc1",                   # benzaldehyde
    "CC(O)c1ccccc1",                   # 1-phenylethan-1-ol
    "OCC(O)CO",                        # glycerol
    "CCOC(=O)CC(=O)OCC",              # diethyl malonate
    "c1ccc(N)cc1",                     # aniline
    "c1ccc(O)cc1",                     # phenol
    "CC(=O)OCC",                       # ethyl acetate
    "CC#N",                            # acetonitrile
    "c1ccc(C(=O)O)cc1",               # benzoic acid
    "CC(=O)C",                         # acetone
    "CCCC(=O)O",                       # butyric acid
    "c1ccc(CC(=O)O)cc1",              # phenylacetic acid
    "OC(=O)CC(O)(CC(=O)O)C(=O)O",     # citric acid
    "c1ccc(NC(=O)C)cc1",              # acetanilide
    "CC(=O)c1ccc(O)cc1",              # 4-hydroxyacetophenone
    "CCOCC",                           # diethyl ether
    "CCO",                             # ethanol
    "CC(O)C",                          # isopropanol
    "CC=O",                            # acetaldehyde
    "OC(c1ccccc1)c1ccccc1",           # benzhydrol
    "c1ccc(C(=O)c2ccccc2)cc1",        # benzophenone
    # Slightly larger drug-like molecules
    "CC(C)(C)c1ccc(C(O)=O)cc1",       # 4-tert-butylbenzoic acid
    "Oc1ccc(O)cc1",                    # hydroquinone
    "CC(=O)Oc1ccc(OC(=O)C)cc1",       # hydroquinone diacetate
    "c1cc(O)c(O)cc1",                  # catechol
    "OC(=O)c1cc(O)c(O)c(O)c1",        # gallic acid
    "C(=O)c1ccc(Cl)cc1",              # 4-chlorobenzaldehyde
    "CC(=O)c1ccc(Cl)cc1",             # 4-chloroacetophenone
    "c1ccc(Br)cc1",                    # bromobenzene
    "c1ccc(I)cc1",                     # iodobenzene
    "CC(Br)c1ccccc1",                  # 1-bromoethylbenzene
    # More complex heterocycles
    "c1ccncc1",                        # pyridine
    "c1cc[nH]c1",                      # pyrrole
    "c1ccoc1",                         # furan
    "c1ccsc1",                         # thiophene
    "c1cnc2ccccc2n1",                  # quinazoline
    "c1ccc2[nH]ccc2c1",               # indole
    "c1ccc2ncccc2c1",                  # quinoline
    "c1ccc2cnccc2c1",                  # isoquinoline
    "c1cnc2[nH]cnc2c1",               # purine skeleton
    "O=c1ccoc2ccccc12",               # coumarin
    # Additional building blocks
    "c1ccc(B(O)O)cc1",                # phenylboronic acid
    "C=Cc1ccccc1",                     # styrene
    "C#Cc1ccccc1",                     # phenylacetylene
    "CC(=O)CC(=O)C",                   # acetylacetone
    "NCC(=O)O",                        # glycine
    "CC(N)C(=O)O",                     # alanine
    "CC(=O)NCC(=O)O",                 # N-acetylglycine
    "OC(=O)C=Cc1ccccc1",              # cinnamic acid
    "c1ccc(C=Cc2ccccc2)cc1",          # stilbene
    "C(=O)c1ccco1",                    # furfural
]


def _generate_training_data() -> tuple[list[list[float]], list[float]]:
    """Generate training data from heuristic-scored retro disconnections.

    For each target SMILES, runs a single-depth retrosynthesis to collect
    template-precursor-score tuples, then extracts features.
    """
    X: list[list[float]] = []
    y: list[float] = []

    random.seed(42)
    np.random.seed(42)

    total_targets = len(_TARGET_SMILES)
    skipped = 0

    for idx, smi in enumerate(_TARGET_SMILES):
        try:
            mol = parse(smi)
        except Exception:
            skipped += 1
            continue

        fgs = detect_functional_groups(mol)

        # Run retrosynthesis at depth 1 with wide beam
        try:
            tree = retrosynthesis(mol, max_depth=1, beam_width=10)
        except Exception:
            skipped += 1
            continue

        disconnections = tree.target.disconnections
        if not disconnections:
            skipped += 1
            continue

        for disc in disconnections:
            # Compute heuristic score (ground truth label)
            heuristic_score = score_disconnection(
                disc.template, disc.precursors, mol
            )

            # Label smoothing: bonus when all purchasable and high score
            label = heuristic_score
            all_purch = all(is_purchasable(p.smiles) for p in disc.precursors)
            if all_purch and label > 60:
                label += 5.0
            label = min(100.0, max(0.0, label))

            try:
                features = extract_retro_features(
                    mol, smi, disc.template, disc.precursors,
                    depth=0, target_fgs=fgs,
                )
                fvec = [features[k] for k in ALL_RETRO_FEATURE_NAMES]
            except Exception:
                continue

            X.append(fvec)
            y.append(label)

            # Data augmentation: 5 copies with Gaussian noise + depth variation
            for aug_i in range(5):
                noisy = list(fvec)
                for j in range(len(noisy)):
                    if noisy[j] != 0.0 and noisy[j] != 1.0:
                        # 5% Gaussian noise on continuous features
                        noisy[j] *= (1.0 + np.random.normal(0, 0.05))
                # Vary depth
                depth_feat_idx = ALL_RETRO_FEATURE_NAMES.index("depth_in_tree")
                noisy[depth_feat_idx] = float(random.randint(0, 4))

                # Slightly noisy label
                noisy_label = label + np.random.normal(0, 1.5)
                noisy_label = min(100.0, max(0.0, noisy_label))

                X.append(noisy)
                y.append(noisy_label)

        if (idx + 1) % 10 == 0:
            print(f"  Processed {idx + 1}/{total_targets} targets "
                  f"({len(X)} samples so far)")

    print(f"\nData generation complete: {len(X)} samples from "
          f"{total_targets - skipped}/{total_targets} targets")
    return X, y


def main():
    """Train the retro scorer model and save it."""
    print("=" * 60)
    print("Training Retrosynthetic Disconnection Scorer")
    print("=" * 60)

    # Check dependencies
    try:
        import joblib
        from sklearn.ensemble import GradientBoostingRegressor
        from sklearn.model_selection import cross_val_score
        from sklearn.metrics import mean_absolute_error
    except ImportError:
        print("ERROR: scikit-learn/joblib required. "
              "Install with: pip install molbuilder[ml]")
        sys.exit(1)

    # Generate training data
    print("\n--- Generating training data ---")
    t0 = time.time()
    X_list, y_list = _generate_training_data()
    t1 = time.time()
    print(f"Data generation took {t1 - t0:.1f}s")

    if len(X_list) < 50:
        print(f"WARNING: Only {len(X_list)} samples generated. "
              "Model may underperform.")

    X = np.array(X_list, dtype=np.float64)
    y = np.array(y_list, dtype=np.float64)

    print(f"\nTraining set: {X.shape[0]} samples, {X.shape[1]} features")
    print(f"Label range: [{y.min():.1f}, {y.max():.1f}], "
          f"mean={y.mean():.1f}, std={y.std():.1f}")

    # Train model
    print("\n--- Training GradientBoostingRegressor ---")
    model = GradientBoostingRegressor(
        n_estimators=150,
        max_depth=5,
        learning_rate=0.1,
        subsample=0.8,
        min_samples_leaf=5,
        random_state=42,
    )

    # 3-fold cross-validation
    print("Running 3-fold cross-validation...")
    cv_scores = cross_val_score(model, X, y, cv=3, scoring="r2")
    print(f"CV R-squared: {cv_scores.mean():.4f} (+/- {cv_scores.std():.4f})")
    print(f"  Per fold: {[f'{s:.4f}' for s in cv_scores]}")

    cv_mae = cross_val_score(
        model, X, y, cv=3, scoring="neg_mean_absolute_error"
    )
    mae_vals = -cv_mae
    print(f"CV MAE: {mae_vals.mean():.2f} (+/- {mae_vals.std():.2f})")

    if cv_scores.mean() < 0.3:
        print("\nWARNING: R-squared < 0.3 -- model may underperform vs heuristic")

    # Final fit on all data
    print("\nFitting final model on all data...")
    model.fit(X, y)

    # Feature importance
    importances = model.feature_importances_
    top_indices = np.argsort(importances)[::-1][:15]
    print("\nTop 15 features by importance:")
    for rank, idx in enumerate(top_indices, 1):
        print(f"  {rank:2d}. {ALL_RETRO_FEATURE_NAMES[idx]:35s} "
              f"{importances[idx]:.4f}")

    # Training set metrics
    y_pred = model.predict(X)
    train_r2 = 1 - np.sum((y - y_pred) ** 2) / np.sum((y - y.mean()) ** 2)
    train_mae = mean_absolute_error(y, y_pred)
    print(f"\nFull training set R-squared: {train_r2:.4f}")
    print(f"Full training set MAE: {train_mae:.2f}")

    # Save model
    model_dict = {
        "score_model": model,
        "feature_names": list(ALL_RETRO_FEATURE_NAMES),
        "version": "1.0",
        "r2_score": float(cv_scores.mean()),
        "mae": float(mae_vals.mean()),
    }

    output_dir = os.path.join(
        os.path.dirname(__file__), "..", "molbuilder", "data"
    )
    output_path = os.path.join(output_dir, "retro_scorer.pkl")

    print(f"\nSaving model to {output_path}")
    joblib.dump(model_dict, output_path, compress=("zlib", 6))

    # Generate SHA-256 sidecar
    sha = hashlib.sha256()
    with open(output_path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            sha.update(chunk)
    sha_path = output_path + ".sha256"
    with open(sha_path, "w") as f:
        f.write(sha.hexdigest() + "\n")
    print(f"SHA-256 sidecar: {sha_path}")

    # Verify round-trip
    loaded = joblib.load(output_path)
    assert list(loaded["feature_names"]) == list(ALL_RETRO_FEATURE_NAMES)
    print(f"\nModel verified: {loaded['version']}, "
          f"R2={loaded['r2_score']:.4f}, MAE={loaded['mae']:.2f}")

    file_size = os.path.getsize(output_path)
    print(f"File size: {file_size / 1024:.1f} KB")
    print("\nDone!")


if __name__ == "__main__":
    main()
