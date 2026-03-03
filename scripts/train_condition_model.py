#!/usr/bin/env python3
"""Train ML condition prediction model from ORD aggregate data.

Produces ``molbuilder/data/condition_model.pkl`` -- a joblib dict
containing four gradient-boosting estimators (temperature, solvent,
catalyst, yield) plus class labels and metadata.

Usage:
    python scripts/train_condition_model.py

Requirements:
    pip install molbuilder[ml]   # scikit-learn >= 1.3, joblib >= 1.3
"""

from __future__ import annotations

import json
import math
import os
import random
import sys

import numpy as np

# Ensure the project root is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from molbuilder.data import load_ord_conditions
from molbuilder.process.ml_features import (
    extract_features,
    ALL_FEATURE_NAMES,
    _CATEGORY_NAMES,
)

# ---------------------------------------------------------------------------
# Representative SMILES per reaction category
# ---------------------------------------------------------------------------
# Each entry: (smiles, category_hint) -- used to build training features.
# Multiple SMILES per ORD reaction type generate diverse feature vectors.
_REPRESENTATIVE_SMILES: dict[str, list[str]] = {
    # -- Coupling --
    "suzuki_coupling": ["c1ccc(B(O)O)cc1", "c1ccc(Br)cc1", "CC(=O)c1ccc(B(O)O)cc1"],
    "heck_coupling": ["C=Cc1ccccc1", "c1ccc(Br)cc1", "C=CC(=O)O"],
    "sonogashira_coupling": ["C#Cc1ccccc1", "c1ccc(Br)cc1", "C#CC(O)C"],
    "stille_coupling": ["c1ccc(Br)cc1", "C=Cc1ccccc1"],
    "negishi_coupling": ["c1ccc(Br)cc1", "C=Cc1ccccc1"],
    "buchwald_hartwig_amination": ["Nc1ccccc1", "c1ccc(Br)cc1", "c1cc(N)ccc1C"],
    "chan_lam_coupling": ["Nc1ccccc1", "c1ccc(B(O)O)cc1"],
    "ullmann_coupling": ["c1ccc(I)cc1", "Oc1ccccc1"],
    "kumada_coupling": ["c1ccc(Br)cc1", "C=Cc1ccccc1"],
    "hiyama_coupling": ["c1ccc(Br)cc1", "C=Cc1ccccc1"],
    "miyaura_borylation": ["c1ccc(Br)cc1", "c1ccc(I)cc1"],
    "ir_ch_borylation": ["c1ccccc1", "c1ccncc1"],
    "tsuji_trost_reaction": ["C=CCC(OC(=O)C)C", "C=CCOC(=O)C"],
    "metallaphotoredox_coupling": ["c1ccc(Br)cc1", "CC(=O)C"],
    "liebeskind_srogl_coupling": ["c1ccc(Br)cc1", "CC=O"],
    "cadiot_chodkiewicz_coupling": ["C#CC", "C#Cc1ccccc1"],
    "pd_ch_activation": ["c1ccccc1", "c1ccncc1"],
    "ni_reductive_coupling": ["c1ccc(Br)cc1", "CC=O"],
    "cuaac_click": ["C#CC", "N=[N+]=[N-]"],
    "photocatalytic_ch_coupling": ["c1ccccc1", "CC(=O)C"],
    # -- Carbonyl --
    "amide_coupling": ["CC(=O)O", "NCC", "c1ccc(C(=O)O)cc1"],
    "evans_aldol": ["CC(=O)CC", "CC=O"],
    "grignard_reaction": ["CC=O", "CC(=O)C", "c1ccc(C=O)cc1"],
    "wittig_olefination": ["CC=O", "c1ccc(C=O)cc1"],
    "horner_wadsworth_emmons": ["CC=O", "c1ccc(C=O)cc1"],
    "aldol_condensation": ["CC(=O)C", "CC=O"],
    "claisen_condensation": ["CC(=O)OCC", "CC(=O)C"],
    "schotten_baumann_acylation": ["CC(=O)Cl", "NCC"],
    "fischer_esterification": ["CC(=O)O", "CCO"],
    "cannizzaro_reaction": ["c1ccc(C=O)cc1"],
    "tishchenko_reaction": ["CC=O"],
    "michael_addition": ["C=CC(=O)C", "CC(=O)CC(=O)C"],
    "robinson_annulation": ["CC(=O)C", "C=CC(=O)C"],
    "julia_kocienski_olefination": ["CC=O", "c1ccc(C=O)cc1"],
    "peterson_olefination": ["CC=O"],
    "tebbe_olefination": ["CC(=O)C"],
    "aza_wittig_reaction": ["CC=O"],
    "knoevenagel_condensation": ["CC=O", "CC(=O)CC(=O)OCC"],
    "weinreb_amide": ["CC(=O)Cl"],
    "reformatsky_reaction": ["CC=O", "BrCC(=O)OCC"],
    "lactamisation": ["NCCCC(=O)O"],
    # -- Oxidation --
    "swern_oxidation": ["OCC", "OC(C)C", "OCc1ccccc1"],
    "dess_martin_oxidation": ["OCC", "OC(C)C"],
    "jones_oxidation": ["OCC", "OC(C)C"],
    "pcc_oxidation": ["OCC"],
    "permanganate_oxidation": ["C=Cc1ccccc1", "OCC"],
    "baeyer_villiger_oxidation": ["CC(=O)C", "c1ccc(C(=O)C)cc1"],
    "sharpless_epoxidation": ["C=CC(O)C"],
    "wacker_oxidation": ["C=CC", "C=Cc1ccccc1"],
    "tpap_oxidation": ["OCC"],
    "ibx_oxidation": ["OCC", "OC(C)C"],
    "pinnick_oxidation": ["CC=O", "c1ccc(C=O)cc1"],
    "riley_oxidation": ["CC(=O)C"],
    "dakin_oxidation": ["c1ccc(C=O)cc1O"],
    "tempo_oxidation": ["OCC"],
    "oppenauer_oxidation": ["OC(C)C"],
    "lemieux_johnson_oxidation": ["C=Cc1ccccc1"],
    "jacobsen_epoxidation": ["C=CC(C)C"],
    "shi_epoxidation": ["C=Cc1ccccc1"],
    "anodic_oxidation": ["c1ccccc1"],
    # -- Reduction --
    "nabh4_reduction": ["CC=O", "CC(=O)C", "c1ccc(C=O)cc1"],
    "lialh4_reduction": ["CC(=O)O", "CC(=O)OCC"],
    "dibal_reduction": ["CC(=O)OCC", "CC#N"],
    "catalytic_hydrogenation": ["C=Cc1ccccc1", "c1ccc([N+](=O)[O-])cc1"],
    "birch_reduction": ["c1ccccc1", "c1ccncc1"],
    "wolff_kishner_reduction": ["CC(=O)C"],
    "clemmensen_reduction": ["CC(=O)C"],
    "luche_reduction": ["C=CC(=O)C"],
    "asymmetric_hydrogenation": ["C=CC(=O)O"],
    "rosenmund_reduction": ["CC(=O)Cl"],
    "lindlar_hydrogenation": ["C#CC", "C#Cc1ccccc1"],
    "cbs_reduction": ["CC(=O)c1ccccc1"],
    "noyori_asymmetric_hydrogenation": ["CC(=O)c1ccccc1"],
    "noyori_transfer_hydrogenation": ["CC(=O)c1ccccc1"],
    "ketoreductase_reduction": ["CC(=O)c1ccccc1"],
    # -- Substitution --
    "sn2_substitution": ["CCBr", "CCCl", "CCI"],
    "sn1_substitution": ["CC(C)(C)Br"],
    "gabriel_synthesis": ["CCBr"],
    "finkelstein_reaction": ["CCCl"],
    "mitsunobu_reaction": ["OCC", "c1ccc(O)cc1"],
    "williamson_ether_synthesis": ["CCBr", "OCC"],
    "n_alkylation": ["NCC", "CCBr"],
    "n_methylation": ["NCC"],
    "friedel_crafts_alkylation": ["c1ccccc1", "CCCl"],
    "friedel_crafts_acylation": ["c1ccccc1", "CC(=O)Cl"],
    "vilsmeier_haack_reaction": ["c1ccccc1"],
    "sandmeyer_reaction": ["c1ccc(N)cc1"],
    "kolbe_schmitt_reaction": ["c1ccc(O)cc1"],
    "duff_reaction": ["c1ccc(O)cc1"],
    "gattermann_koch_reaction": ["c1ccccc1"],
    "birch_alkylation": ["c1ccccc1"],
    "balz_schiemann_reaction": ["c1ccc(N)cc1"],
    "reimer_tiemann_reaction": ["c1ccc(O)cc1"],
    # -- Elimination --
    "e2_elimination": ["CC(Br)CC", "CCCBr"],
    "e1_elimination": ["CC(C)(C)Br"],
    "hofmann_elimination": ["CC(C)N(C)(C)C"],
    "cope_elimination": ["CC([N+](C)(C)[O-])CC"],
    "shapiro_reaction": ["CC(=NNS(=O)(=O)c1ccc(C)cc1)C"],
    "bamford_stevens_reaction": ["CC(=NNS(=O)(=O)c1ccc(C)cc1)C"],
    # -- Protection / Deprotection --
    "boc_protection": ["NCC", "Nc1ccccc1"],
    "boc_deprotection": ["CC(C)(C)OC(=O)NCC"],
    "fmoc_protection": ["NCC"],
    "fmoc_deprotection": ["NCC"],
    "tbs_protection": ["OCC"],
    "tbs_deprotection": ["OCC"],
    "benzyl_protection": ["OCC"],
    "hydrogenolysis": ["OCc1ccccc1"],
    "pmb_protection": ["OCC"],
    "pmb_deprotection": ["OCC"],
    "acetonide_protection": ["OCC(O)C"],
    "thp_protection": ["OCC"],
    "acetylation": ["OCC", "NCC"],
    # -- Rearrangement --
    "claisen_rearrangement": ["C=CCOC=C"],
    "beckmann_rearrangement": ["CC(=NO)C"],
    "cope_rearrangement": ["C=CCC=CC"],
    "pinacol_rearrangement": ["OC(C)(C)C(O)(C)C"],
    "curtius_rearrangement": ["CC(=O)Cl"],
    "hofmann_rearrangement": ["CC(=O)N"],
    "wolff_rearrangement": ["CC(=O)C=[N+]=[N-]"],
    "favorskii_rearrangement": ["CC(=O)C(Cl)C"],
    "fries_rearrangement": ["c1ccc(OC(=O)C)cc1"],
    "schmidt_reaction": ["CC(=O)O"],
    "overman_rearrangement": ["C=CC(OC(=N)CC)C"],
    # -- Pericyclic --
    "diels_alder": ["C=CC=CC", "C=CC(=O)C"],
    # -- Misc / Heterocycle synthesis --
    "traube_synthesis": ["NC(=O)c1nc[nH]c1N"],
    "pyrimidine_synthesis": ["CC(=O)CC(=O)C"],
    "minisci_reaction": ["c1ccncc1"],
    "barton_decarboxylation": ["CC(=O)O"],
    "norrish_type_i": ["CC(=O)C"],
    "norrish_type_ii": ["CCCC(=O)C"],
    "wohl_ziegler_bromination": ["CC=CC"],
    "paterno_buchi_reaction": ["CC=O", "C=CC"],
    "appel_reaction": ["OCC"],
    "staudinger_reduction": ["N=[N+]=[N-]"],
    "corey_chaykovsky_reaction": ["CC=O"],
    "ring_closing_metathesis": ["C=CCC=C"],
    "cross_metathesis": ["C=CC", "C=Cc1ccccc1"],
    "romp": ["C1CC=CC1"],
    "hantzsch_pyridine_synthesis": ["CC=O", "CC(=O)CC(=O)OCC"],
    "skraup_synthesis": ["Nc1ccccc1"],
    "doebner_miller_synthesis": ["Nc1ccccc1", "CC=O"],
    "friedlander_synthesis": ["Nc1ccccc1C=O"],
    "paal_knorr_synthesis": ["CC(=O)CCC(=O)C"],
    "fischer_indole_synthesis": ["c1ccc(NN)cc1"],
    "madelung_indole_synthesis": ["c1ccc(NC(=O)C)cc1"],
    "leimgruber_batcho_synthesis": ["c1ccc([N+](=O)[O-])cc1C"],
    "chichibabin_synthesis": ["c1ccncc1"],
    "knorr_pyrazole_synthesis": ["CC(=O)CC(=O)C"],
    "gewald_reaction": ["CC(=O)C", "CC#N"],
    "combes_synthesis": ["Nc1ccccc1", "CC(=O)CC(=O)C"],
    "strecker_synthesis": ["CC=O"],
    "mannich_reaction": ["CC(=O)C", "NCC"],
    "passerini_reaction": ["CC=O", "CC(=O)O"],
    "ugi_reaction": ["CC=O", "NCC", "CC(=O)O"],
    "baylis_hillman_reaction": ["CC=O", "C=CC(=O)OC"],
    "henry_reaction": ["CC=O"],
    "enzymatic_resolution": ["CC(O)C(=O)OCC"],
    "sharpless_dihydroxylation": ["C=Cc1ccccc1"],
    "prilezhaev_epoxidation": ["C=CC"],
    "electrophilic_addition": ["C=CC"],
    "radical_addition": ["C=CC"],
    "hydroboration_oxidation": ["C=CC"],
    "upjohn_dihydroxylation": ["C=Cc1ccccc1"],
    "halogenation": ["C=CC"],
    "acid_catalysed_hydration": ["C=CC"],
    "ozonolysis": ["C=Cc1ccccc1"],
    "simmons_smith_reaction": ["C=CC"],
    # -- General categories (these get broad representative SMILES) --
    "substitution_general": ["CCBr", "CCCl"],
    "elimination_general": ["CC(Br)CC"],
    "addition_general": ["C=CC"],
    "oxidation_general": ["OCC", "OC(C)C"],
    "reduction_general": ["CC=O", "CC(=O)C"],
    "cross_coupling_general": ["c1ccc(Br)cc1", "c1ccc(B(O)O)cc1"],
    "carbonyl_general": ["CC=O", "CC(=O)C"],
    "protection_general": ["OCC", "NCC"],
    "deprotection_general": ["CC(C)(C)OC(=O)NCC"],
    "rearrangement_general": ["C=CCOC=C"],
    "radical_general": ["CC=CC"],
    "pericyclic_general": ["C=CC=CC"],
    "polymerization_general": ["C=CC", "C=Cc1ccccc1"],
    "misc_general": ["c1ccccc1", "CC=O"],
}

# Map ORD reaction keys to MolBuilder category names for feature extraction
_ORD_KEY_TO_CATEGORY: dict[str, str] = {}
for key in _REPRESENTATIVE_SMILES:
    # Match against category names
    for cat_name in _CATEGORY_NAMES:
        if cat_name.lower() in key:
            _ORD_KEY_TO_CATEGORY[key] = cat_name
            break
    if key not in _ORD_KEY_TO_CATEGORY:
        # Derive from suffix patterns
        if "coupling" in key or "borylation" in key or "amination" in key:
            _ORD_KEY_TO_CATEGORY[key] = "COUPLING"
        elif "oxidation" in key or "epoxidation" in key:
            _ORD_KEY_TO_CATEGORY[key] = "OXIDATION"
        elif "reduction" in key or "hydrogenation" in key:
            _ORD_KEY_TO_CATEGORY[key] = "REDUCTION"
        elif "substitution" in key or "alkylation" in key or "acylation" in key:
            _ORD_KEY_TO_CATEGORY[key] = "SUBSTITUTION"
        elif "elimination" in key:
            _ORD_KEY_TO_CATEGORY[key] = "ELIMINATION"
        elif "protection" in key or "deprotection" in key or "acetylation" in key:
            _ORD_KEY_TO_CATEGORY[key] = "PROTECTION"
        elif "rearrangement" in key:
            _ORD_KEY_TO_CATEGORY[key] = "REARRANGEMENT"
        elif "addition" in key or "dihydroxylation" in key or "hydroboration" in key:
            _ORD_KEY_TO_CATEGORY[key] = "ADDITION"
        elif "aldol" in key or "condensation" in key or "olefination" in key:
            _ORD_KEY_TO_CATEGORY[key] = "CARBONYL"
        elif "diels_alder" in key:
            _ORD_KEY_TO_CATEGORY[key] = "PERICYCLIC"
        elif "metathesis" in key or "romp" in key or "polymerization" in key:
            _ORD_KEY_TO_CATEGORY[key] = "POLYMERIZATION"
        elif "radical" in key or "barton" in key or "minisci" in key:
            _ORD_KEY_TO_CATEGORY[key] = "RADICAL"
        else:
            _ORD_KEY_TO_CATEGORY[key] = "MISC"


def _build_training_data(
    ord_data: dict,
    min_n: int = 10,
    augment_per_smiles: int = 8,
    max_solvent_classes: int = 20,
    max_catalyst_classes: int = 25,
) -> tuple[list[list[float]], list[float], list[int], list[int], list[float]]:
    """Build feature matrix and target arrays from ORD data.

    Returns (X, y_temp, y_solvent, y_catalyst, y_yield,
             solvent_classes, catalyst_classes).
    """
    reactions = ord_data.get("reactions", {})

    # Collect solvents/catalysts ranked by total occurrence count
    solvent_counts: dict[str, int] = {}
    catalyst_counts: dict[str, int] = {}
    for rxn_key, stats in reactions.items():
        if stats["n"] < min_n:
            continue
        for s in stats.get("solvents", []):
            solvent_counts[s["name"]] = solvent_counts.get(s["name"], 0) + s["count"]
        for c in stats.get("catalysts", []):
            catalyst_counts[c["name"]] = catalyst_counts.get(c["name"], 0) + c["count"]

    # Keep top-N by frequency, map the rest to "other"
    sorted_solvents = sorted(solvent_counts, key=solvent_counts.get, reverse=True)
    all_solvents = sorted_solvents[:max_solvent_classes]
    if "other" not in all_solvents:
        all_solvents.append("other")

    sorted_catalysts = sorted(catalyst_counts, key=catalyst_counts.get, reverse=True)
    all_catalysts = sorted_catalysts[:max_catalyst_classes]
    if "none" not in all_catalysts:
        all_catalysts.append("none")
    if "other" not in all_catalysts:
        all_catalysts.append("other")

    solvent_to_idx = {s: i for i, s in enumerate(all_solvents)}
    catalyst_to_idx = {c: i for i, c in enumerate(all_catalysts)}
    other_solvent_idx = solvent_to_idx["other"]
    other_catalyst_idx = catalyst_to_idx["other"]

    X: list[list[float]] = []
    y_temp: list[float] = []
    y_solvent: list[int] = []
    y_catalyst: list[int] = []
    y_yield: list[float] = []

    random.seed(42)
    np.random.seed(42)

    for rxn_key, stats in reactions.items():
        if stats["n"] < min_n:
            continue

        smiles_list = _REPRESENTATIVE_SMILES.get(rxn_key, [])
        if not smiles_list:
            continue

        category = _ORD_KEY_TO_CATEGORY.get(rxn_key)

        # Extract target values from ORD
        temp_median = float(stats["temperature_C"]["median"])
        temp_iqr = float(
            stats["temperature_C"]["p75"] - stats["temperature_C"]["p25"]
        )

        # Top solvent and catalyst (mapped to known class or "other")
        solvents = stats.get("solvents", [])
        top_solvent = solvents[0]["name"] if solvents else "THF"
        top_solvent_idx = solvent_to_idx.get(top_solvent, other_solvent_idx)

        catalysts = stats.get("catalysts", [])
        top_catalyst = catalysts[0]["name"] if catalysts else "none"
        top_catalyst_idx = catalyst_to_idx.get(top_catalyst, other_catalyst_idx)

        yield_median = float(stats.get("yield_pct", {}).get("median", 70))
        yield_iqr = float(
            stats.get("yield_pct", {}).get("p75", 85)
            - stats.get("yield_pct", {}).get("p25", 55)
        )

        for smi in smiles_list:
            try:
                features = extract_features(smi, category, scale_kg=1.0)
                base_vec = [features[k] for k in ALL_FEATURE_NAMES]
            except Exception:
                continue

            # Generate augmented samples with noise
            for _ in range(augment_per_smiles):
                vec = list(base_vec)
                # Add small Gaussian noise to continuous descriptors (first 8)
                for j in range(8):
                    vec[j] += np.random.normal(0, abs(vec[j]) * 0.05 + 0.01)

                # Vary scale_kg for the last feature
                log_scale = math.log1p(np.random.lognormal(0, 1.5))
                vec[-1] = log_scale

                X.append(vec)

                # Temperature target with noise within IQR
                t_noise = np.random.normal(0, max(temp_iqr * 0.3, 1.0))
                y_temp.append(temp_median + t_noise)

                # Solvent: mostly top, sometimes second
                if solvents and len(solvents) > 1 and random.random() < 0.2:
                    alt = solvents[1]["name"]
                    y_solvent.append(solvent_to_idx.get(alt, other_solvent_idx))
                else:
                    y_solvent.append(top_solvent_idx)

                # Catalyst: mostly top, sometimes second
                if catalysts and len(catalysts) > 1 and random.random() < 0.15:
                    alt_cat = catalysts[1]["name"]
                    y_catalyst.append(
                        catalyst_to_idx.get(alt_cat, other_catalyst_idx)
                    )
                else:
                    y_catalyst.append(top_catalyst_idx)

                # Yield target with noise
                y_noise = np.random.normal(0, max(yield_iqr * 0.25, 2.0))
                y_yield.append(max(0.0, min(100.0, yield_median + y_noise)))

    return X, y_temp, y_solvent, y_catalyst, y_yield, all_solvents, all_catalysts


def main() -> None:
    """Train and save the condition prediction model."""
    from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier
    from sklearn.model_selection import cross_val_score
    import joblib

    print("Loading ORD condition data...")
    ord_data = load_ord_conditions()
    reactions = ord_data.get("reactions", {})
    n_reactions = sum(1 for v in reactions.values() if v["n"] >= 10)
    print(f"  {n_reactions} reaction types with n >= 10")

    print("Building training data...")
    result = _build_training_data(ord_data, min_n=10, augment_per_smiles=5)
    X, y_temp, y_solvent, y_catalyst, y_yield, solvent_classes, catalyst_classes = result
    X_arr = np.array(X, dtype=np.float64)
    y_temp_arr = np.array(y_temp, dtype=np.float64)
    y_solvent_arr = np.array(y_solvent, dtype=np.int32)
    y_catalyst_arr = np.array(y_catalyst, dtype=np.int32)
    y_yield_arr = np.array(y_yield, dtype=np.float64)

    print(f"  Training samples: {len(X)}")
    print(f"  Features: {len(ALL_FEATURE_NAMES)}")
    print(f"  Solvent classes: {len(solvent_classes)}")
    print(f"  Catalyst classes: {len(catalyst_classes)}")

    # -- Train temperature regressor --
    print("\nTraining temperature model...")
    temp_model = GradientBoostingRegressor(
        n_estimators=50, max_depth=3, learning_rate=0.15,
        subsample=0.8, random_state=42,
    )
    temp_scores = cross_val_score(temp_model, X_arr, y_temp_arr, cv=3, scoring="r2")
    print(f"  Temperature R2 (3-fold CV): {temp_scores.mean():.3f} +/- {temp_scores.std():.3f}")
    temp_model.fit(X_arr, y_temp_arr)

    # -- Train solvent classifier --
    print("Training solvent model...")
    solvent_model = GradientBoostingClassifier(
        n_estimators=50, max_depth=3, learning_rate=0.15,
        subsample=0.8, random_state=42,
    )
    solvent_scores = cross_val_score(solvent_model, X_arr, y_solvent_arr, cv=3, scoring="accuracy")
    print(f"  Solvent accuracy (3-fold CV): {solvent_scores.mean():.3f} +/- {solvent_scores.std():.3f}")
    solvent_model.fit(X_arr, y_solvent_arr)

    # -- Train catalyst classifier --
    print("Training catalyst model...")
    catalyst_model = GradientBoostingClassifier(
        n_estimators=50, max_depth=3, learning_rate=0.15,
        subsample=0.8, random_state=42,
    )
    catalyst_scores = cross_val_score(catalyst_model, X_arr, y_catalyst_arr, cv=3, scoring="accuracy")
    print(f"  Catalyst accuracy (3-fold CV): {catalyst_scores.mean():.3f} +/- {catalyst_scores.std():.3f}")
    catalyst_model.fit(X_arr, y_catalyst_arr)

    # -- Train yield regressor --
    print("Training yield model...")
    yield_model = GradientBoostingRegressor(
        n_estimators=50, max_depth=3, learning_rate=0.15,
        subsample=0.8, random_state=42,
    )
    yield_scores = cross_val_score(yield_model, X_arr, y_yield_arr, cv=3, scoring="r2")
    print(f"  Yield R2 (3-fold CV): {yield_scores.mean():.3f} +/- {yield_scores.std():.3f}")
    yield_model.fit(X_arr, y_yield_arr)

    # -- Save model --
    model_dict = {
        "temperature_model": temp_model,
        "solvent_model": solvent_model,
        "catalyst_model": catalyst_model,
        "yield_model": yield_model,
        "solvent_classes": solvent_classes,
        "catalyst_classes": catalyst_classes,
        "feature_names": ALL_FEATURE_NAMES,
        "version": "1.0",
    }

    out_path = os.path.join(
        os.path.dirname(__file__), "..", "molbuilder", "data", "condition_model.pkl"
    )
    out_path = os.path.normpath(out_path)
    joblib.dump(model_dict, out_path, compress=("zlib", 6))

    file_size = os.path.getsize(out_path)
    print(f"\nModel saved to {out_path}")
    print(f"File size: {file_size / 1024:.1f} KB")
    print("\nDone!")


if __name__ == "__main__":
    main()
