"""Generate a markdown validation report for the drug retrosynthesis suite.

Usage:
    python tests/validation_report.py

Outputs:
    docs/validation/drug_validation_results.md
"""

import sys
import time
from pathlib import Path

# Ensure the package is importable
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis


# Same drug list used in test_drug_validation.py
DRUG_DATA = [
    ("aspirin", "CC(=O)Oc1ccccc1C(O)=O", 1, "ester"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(O)=O", 1, None),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1", 1, "amide"),
    ("caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O", 1, None),
    ("naproxen", "COc1ccc2cc(CC(C)C(O)=O)ccc2c1", 1, None),
    ("atorvastatin", "CC(C)c1c(C(=O)Nc2ccccc2)c(c(c3ccc(F)cc3)n1CC(O)CC(O)CC(O)=O)c4ccccc4", 1, None),
    ("metformin", "CN(C)C(=N)NC(N)=N", 1, None),
    ("omeprazole", "COc1ccc2nc(CS(=O)c3ncc(C)c(OC)c3C)[nH]c2c1", 1, None),
    ("losartan", "CCCCc1nc(Cl)c(n1Cc2ccc(cc2)c3ccccc3C=O)CO", 1, None),
    ("lisinopril", "NCCCC(NC(CCc1ccccc1)C(O)=O)C(=O)N1CCCC1C(O)=O", 1, None),
    ("amoxicillin", "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(O)=O", 1, None),
    ("ciprofloxacin", "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O", 1, None),
    ("sildenafil", "CCCc1nn(C)c2c(=O)[nH]c(nc12)c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1", 1, None),
    ("sertraline", "CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc21", 1, None),
    ("warfarin", "CC(=O)CC(c1ccccc1)c2c(O)c3ccccc3oc2=O", 1, None),
    ("fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", 1, None),
    ("diazepam", "CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc31", 1, None),
    ("gabapentin", "NCC1(CC(O)=O)CCCCC1", 1, None),
    ("metoprolol", "COCCc1ccc(OCC(O)CNC(C)C)cc1", 1, None),
    ("amlodipine", "CCOC(=O)C1=C(COCCN)NC(C)=C(C1c1ccccc1Cl)C(=O)OC", 1, None),
    ("furosemide", "NS(=O)(=O)c1cc(C(O)=O)c(NCc2ccco2)cc1Cl", 1, None),
    ("cetirizine", "OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1", 1, None),
    ("loratadine", "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3ncccc32)CC1", 1, None),
    ("captopril", "CC(CS)C(=O)N1CCCC1C(O)=O", 1, None),
    ("enalapril", "CCOC(=O)C(CCc1ccccc1)NC(C)C(=O)N1CCCC1C(O)=O", 1, None),
    ("valsartan", "CCCCC(=O)N(Cc1ccc(c(c1)c1ccccc1C=O)c1nn[nH]n1)C(C(C)C)C(O)=O", 1, None),
    ("methotrexate", "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)NC(CCC(O)=O)C(O)=O)cc1", 1, None),
    ("tamoxifen", "CCC(=C(c1ccccc1)c1ccccc1)c1ccc(OCCN(C)C)cc1", 1, None),
    ("prednisone", "CC12CC(=O)C3C(CCC4=CC(=O)C=CC34C)C1CCC2(O)C(=O)CO", 1, None),
    ("dexamethasone", "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO", 1, None),
    ("citalopram", "N#Cc1ccc2c(c1)C(CCCN1CCCC1)(OC2)c1ccc(F)cc1", 1, None),
    ("venlafaxine", "COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1", 1, None),
    ("bupropion", "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1", 1, None),
    ("lamotrigine", "Nc1nnc(c(N)n1)c1cccc(Cl)c1Cl", 1, None),
    ("montelukast", "CC(C)(O)c1ccccc1CCC(SCC1(CC(O)=O)CC1)c1cccc(c1)C=Cc1ccc2ccc(Cl)cc2n1", 1, None),
    ("pantoprazole", "COc1ccnc(CS(=O)c2nc3cc(OC(F)F)ccc3[nH]2)c1OC", 1, None),
    ("duloxetine", "CNCC(Oc1cccc2ccccc12)c1cccs1", 1, None),
    ("aripiprazole", "O=c1[nH]c2ccccc2n1CCCCN1CCN(CC1)c1ccc(Cl)cc1OC1CCCC1", 1, None),
    ("rosuvastatin", "CC(c1nc(N(C)S(C)(=O)=O)nc1c1ccc(F)cc1)c1ccc(cc1)C(O)CC(O)CC(O)=O", 1, None),
    ("esomeprazole", "COc1ccc2nc(CS(=O)c3ncc(C)c(OC)c3C)[nH]c2c1", 1, None),
    ("clopidogrel", "COC(=O)C(c1ccccc1Cl)N1CCc2sccc2C1", 1, None),
    ("levothyroxine", "NC(Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O", 1, None),
    ("allopurinol", "O=c1[nH]cnc2[nH]ncc12", 1, None),
    ("ranitidine", "CNC(NCCSCc1ccc(CN(C)C)o1)=C[N+](=O)[O-]", 1, None),
    ("simvastatin", "CCC(C)(C)C(=O)OC1CC(O)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C21", 1, None),
    ("propranolol", "CC(C)NCC(O)COc1cccc2ccccc12", 1, None),
    ("methylphenidate", "OC(=O)C(c1ccccc1)C1CCCCN1C", 1, None),
    ("clonazepam", "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1", 1, None),
    ("azithromycin", "CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(O)CC(C)CN(C)C(C)C(O)C1(C)O", None, None),
    ("doxycycline", "CC1C2C(O)C3C(=C(O)c4c(O)cccc4C3(C)O)C(=O)C2(O)C(N(C)C)C1O", None, None),
]


def run_validation():
    """Run retrosynthesis on all drugs and collect results."""
    results = []
    for drug_name, smiles, min_steps, key_reaction in DRUG_DATA:
        entry = {
            "name": drug_name,
            "smiles": smiles,
            "min_steps": min_steps,
            "key_reaction": key_reaction,
            "parsed": False,
            "tree_built": False,
            "has_disconnections": False,
            "num_disconnections": 0,
            "num_functional_groups": 0,
            "is_purchasable": False,
            "key_reaction_found": False,
            "elapsed_s": 0.0,
            "error": None,
        }

        t0 = time.perf_counter()
        try:
            mol = parse(smiles)
            entry["parsed"] = True
        except Exception as exc:
            entry["error"] = f"Parse error: {exc}"
            entry["elapsed_s"] = time.perf_counter() - t0
            results.append(entry)
            continue

        try:
            tree = retrosynthesis(mol, max_depth=5, beam_width=5)
            entry["tree_built"] = True
            target = tree.target
            entry["is_purchasable"] = target.is_purchasable
            entry["num_functional_groups"] = len(target.functional_groups)
            entry["num_disconnections"] = len(target.disconnections)
            entry["has_disconnections"] = len(target.disconnections) > 0

            if key_reaction and target.disconnections:
                reaction_names = " ".join(
                    d.template.name.lower() for d in target.disconnections
                )
                entry["key_reaction_found"] = key_reaction.lower() in reaction_names
        except Exception as exc:
            entry["error"] = f"Retrosynthesis error: {exc}"

        entry["elapsed_s"] = time.perf_counter() - t0
        results.append(entry)

    return results


def generate_report(results):
    """Generate a markdown report from validation results."""
    total = len(results)
    parsed = sum(1 for r in results if r["parsed"])
    tree_built = sum(1 for r in results if r["tree_built"])
    has_disconnections = sum(1 for r in results if r["has_disconnections"])
    has_fg = sum(1 for r in results if r["num_functional_groups"] > 0)
    errors = sum(1 for r in results if r["error"])
    total_time = sum(r["elapsed_s"] for r in results)
    avg_time = total_time / total if total else 0

    # Success rate: tree built / total
    success_rate = tree_built / total * 100 if total else 0
    # Route rate: has disconnections / total non-macrocycle
    non_macro = [r for r in results if r["min_steps"] is not None]
    route_rate = (
        sum(1 for r in non_macro if r["has_disconnections"]) / len(non_macro) * 100
        if non_macro
        else 0
    )

    lines = [
        "# MolBuilder Drug Validation Report",
        "",
        f"**Generated**: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"**MolBuilder Version**: 1.1.1",
        f"**Drugs Tested**: {total}",
        "",
        "## Summary",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Total drugs | {total} |",
        f"| SMILES parsed | {parsed}/{total} ({parsed/total*100:.0f}%) |",
        f"| Retrosynthesis tree built | {tree_built}/{total} ({success_rate:.0f}%) |",
        f"| Disconnections found | {has_disconnections}/{total} ({has_disconnections/total*100:.0f}%) |",
        f"| Functional groups detected | {has_fg}/{total} ({has_fg/total*100:.0f}%) |",
        f"| Route coverage (non-macrocycle) | {route_rate:.0f}% |",
        f"| Errors | {errors} |",
        f"| Total runtime | {total_time:.1f}s |",
        f"| Average per drug | {avg_time:.2f}s |",
        "",
        "## Per-Drug Results",
        "",
        "| Drug | Parsed | Tree | FGs | Disconnections | Time (s) | Notes |",
        "|------|--------|------|-----|----------------|----------|-------|",
    ]

    for r in results:
        parsed_str = "Y" if r["parsed"] else "N"
        tree_str = "Y" if r["tree_built"] else "N"
        fg_count = str(r["num_functional_groups"])
        disc_count = str(r["num_disconnections"])
        time_str = f"{r['elapsed_s']:.2f}"

        notes = []
        if r["error"]:
            notes.append(f"Error: {r['error']}")
        if r["min_steps"] is None:
            notes.append("macrocycle")
        if r["key_reaction"] and r["key_reaction_found"]:
            notes.append(f"key: {r['key_reaction']}")
        if r["is_purchasable"]:
            notes.append("purchasable")

        note_str = "; ".join(notes) if notes else "-"
        lines.append(
            f"| {r['name']} | {parsed_str} | {tree_str} | {fg_count} | "
            f"{disc_count} | {time_str} | {note_str} |"
        )

    lines.extend([
        "",
        "## Methodology",
        "",
        "Each drug's canonical SMILES was parsed by MolBuilder's SMILES parser, "
        "then fed into the retrosynthetic analysis engine with `max_depth=5` and "
        "`beam_width=5`. The engine identifies functional groups, proposes "
        "disconnections using the reaction knowledge base (181 templates), and "
        "builds a retrosynthetic tree toward purchasable starting materials "
        "(250+ entries).",
        "",
        "Complex macrocycles (azithromycin, doxycycline) are expected to be "
        "challenging and are evaluated on whether the analysis completes without "
        "errors rather than on route quality.",
        "",
        "## Compliance Notes",
        "",
        "This validation report supports IQ/OQ documentation for 21 CFR Part 11 "
        "compliance. The test suite is fully automated and reproducible via:",
        "",
        "```bash",
        "python -m pytest tests/test_drug_validation.py -v",
        "python tests/validation_report.py",
        "```",
        "",
    ])

    return "\n".join(lines)


def main():
    print("Running drug validation suite...")
    results = run_validation()
    report = generate_report(results)

    out_path = ROOT / "docs" / "validation" / "drug_validation_results.md"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report, encoding="utf-8")
    print(f"Report written to {out_path}")

    # Print summary to stdout
    total = len(results)
    tree_built = sum(1 for r in results if r["tree_built"])
    print(f"\n  {tree_built}/{total} drugs: retrosynthesis tree built successfully")
    has_disc = sum(1 for r in results if r["has_disconnections"])
    print(f"  {has_disc}/{total} drugs: disconnections found")
    errors = sum(1 for r in results if r["error"])
    if errors:
        print(f"  {errors} errors encountered")
    else:
        print("  No errors encountered")


if __name__ == "__main__":
    main()
