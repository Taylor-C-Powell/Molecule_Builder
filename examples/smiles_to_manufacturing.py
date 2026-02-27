#!/usr/bin/env python3
"""SMILES to Manufacturing: end-to-end process chemistry demo.

Usage:
    python examples/smiles_to_manufacturing.py CCO
    python examples/smiles_to_manufacturing.py "c1ccccc1"
    python examples/smiles_to_manufacturing.py "CC(=O)O"

This script demonstrates the full MolBuilder pipeline:
1. Parse a SMILES string into a Molecule
2. Detect functional groups
3. Run retrosynthetic analysis
4. Extract the best synthesis route
5. Optimize reaction conditions
6. Estimate manufacturing cost
7. Assess safety hazards
8. Print formatted reports
"""

import sys

from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles
from molbuilder.reactions.functional_group_detect import detect_functional_groups
from molbuilder.reactions.retrosynthesis import retrosynthesis
from molbuilder.reactions.synthesis_route import extract_best_route
from molbuilder.process.costing import estimate_cost
from molbuilder.process.safety import assess_safety
from molbuilder.reports import (
    generate_molecule_report,
    generate_synthesis_report,
    generate_cost_report,
    generate_safety_report,
)


def main():
    smiles = sys.argv[1] if len(sys.argv) > 1 else "CC(=O)O"

    # 1. Parse SMILES
    print(f"Parsing SMILES: {smiles}")
    mol = parse(smiles)
    mol.name = smiles
    canonical = to_smiles(mol)
    print(f"Canonical SMILES: {canonical}")
    print(f"Atoms: {len(mol.atoms)}, Bonds: {len(mol.bonds)}")
    print()

    # 2. Detect functional groups
    fgs = detect_functional_groups(mol)
    if fgs:
        print(f"Functional groups: {', '.join(fg.name for fg in fgs)}")
    else:
        print("No functional groups detected (simple hydrocarbon)")
    print()

    # 3. Retrosynthesis
    print("Running retrosynthetic analysis...")
    tree = retrosynthesis(mol)
    route = extract_best_route(tree)
    if route is None:
        print("No synthesis route found (molecule may be a starting material).")
        print()
        print(generate_molecule_report(mol))
        return

    print(f"Found route: {route.total_steps} steps, "
          f"overall yield {route.overall_yield:.1%}")
    print()

    # 4. Synthesis report
    print(generate_synthesis_report(route))
    print()

    # 5. Cost estimate
    try:
        cost = estimate_cost(route.steps, scale_kg=1.0)
        print(generate_cost_report(cost))
        print()
    except Exception as e:
        print(f"Cost estimation skipped: {e}")
        print()

    # 6. Safety assessment
    try:
        safety = assess_safety(route.steps)
        print(generate_safety_report(safety))
    except Exception as e:
        print(f"Safety assessment skipped: {e}")


if __name__ == "__main__":
    main()
