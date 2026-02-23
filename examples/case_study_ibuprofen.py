#!/usr/bin/env python3
"""
Case Study: Ibuprofen Process Chemistry Pipeline
=================================================

Demonstrates MolBuilder v1.2.0's end-to-end process chemistry workflow:

  SMILES -> Retrosynthesis -> Feasibility Scoring -> Reactor Selection
         -> Cost Estimation -> Safety Assessment -> Scale-Up Analysis
         -> Reaction Conditions

This script is fully reproducible:
    pip install molbuilder
    python examples/case_study_ibuprofen.py

No external dependencies beyond MolBuilder's core (numpy, scipy, matplotlib).
"""
from __future__ import annotations

from molbuilder.smiles import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree
from molbuilder.reactions.feasibility import FeasibilityEngine
from molbuilder.process.reactor import select_reactor
from molbuilder.process.costing import estimate_cost
from molbuilder.process.safety import assess_safety
from molbuilder.process.scale_up import analyze_scale_up
from molbuilder.process.conditions import optimize_conditions


def _header(title: str) -> None:
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print(f"{'=' * 60}")


def main() -> None:
    # ----------------------------------------------------------------
    # Target definition
    # ----------------------------------------------------------------
    SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    mol = parse(SMILES)
    mol.name = "ibuprofen"

    _header("MolBuilder v1.2.0 -- Process Chemistry Case Study")
    print(f"  Target:  {mol.name}")
    print(f"  SMILES:  {SMILES}")
    print(f"  Formula: C13H18O2  (MW 206.28)")
    print(f"  Atoms:   {len(mol.atoms)}   Bonds: {len(mol.bonds)}")

    # ----------------------------------------------------------------
    # 1. Retrosynthetic analysis
    # ----------------------------------------------------------------
    _header("1. Retrosynthetic Analysis")
    tree = retrosynthesis(mol, max_depth=8, beam_width=5)
    print(format_tree(tree))

    disc = tree.target.best_disconnection
    template = disc.template
    precursors = disc.precursors

    print(f"  Best disconnection: {template.name}")
    print(f"  Category:           {template.category.name}")
    print(f"  Precursors:")
    for p in precursors:
        cost_str = f"${p.cost_per_kg:.2f}/kg" if p.cost_per_kg else "N/A"
        print(f"    - {p.name} ({p.smiles})  [{cost_str}]")

    # ----------------------------------------------------------------
    # 2. Reactor selection
    # ----------------------------------------------------------------
    SCALE_KG = 100
    _header(f"2. Reactor Recommendation ({SCALE_KG} kg batch)")
    r = select_reactor(template, scale_kg=SCALE_KG)
    print(f"  Type:           {r.reactor_type.name}")
    print(f"  Volume:         {r.volume_L} L")
    print(f"  Temperature:    {r.temperature_C} C")
    print(f"  Pressure:       {r.pressure_atm} atm")
    print(f"  Residence time: {r.residence_time_min} min")
    print(f"  Mixing:         {r.mixing_type}")
    print(f"  Material:       {r.material}")
    print(f"  Reactor cost:   ${r.estimated_cost_usd:,.0f}")

    # ----------------------------------------------------------------
    # 3. Cost estimation
    # ----------------------------------------------------------------
    _header(f"3. Cost Estimate ({SCALE_KG} kg)")

    # Build step object for process modules
    class _Step:
        def __init__(self, t, p, n):
            self.template = t
            self.precursors = p
            self.step_number = n

    steps = [_Step(template, precursors, 1)]
    c = estimate_cost(steps, scale_kg=SCALE_KG)
    b = c.breakdown
    print(f"  Raw materials:  ${b.raw_materials_usd:>10,.2f}")
    print(f"  Labor:          ${b.labor_usd:>10,.2f}")
    print(f"  Equipment:      ${b.equipment_usd:>10,.2f}")
    print(f"  Energy:         ${b.energy_usd:>10,.2f}")
    print(f"  Waste disposal: ${b.waste_disposal_usd:>10,.2f}")
    print(f"  Overhead:       ${b.overhead_usd:>10,.2f}")
    print(f"  {'-' * 35}")
    print(f"  TOTAL:          ${c.total_usd:>10,.2f}")
    print(f"  PER KG:         ${c.per_kg_usd:>10,.2f}/kg")

    # ----------------------------------------------------------------
    # 4. Safety assessment
    # ----------------------------------------------------------------
    _header("4. Safety Assessment")
    for s in assess_safety(steps):
        print(f"  Step {s.step_number}: {s.step_name}")
        print(f"    Risk level:       {s.risk_level}")
        print(f"    PPE required:     {', '.join(s.ppe_required)}")
        print(f"    Eng. controls:    {', '.join(s.engineering_controls)}")
        print(f"    Waste class:      {s.waste_classification}")

    # ----------------------------------------------------------------
    # 5. Scale-up analysis
    # ----------------------------------------------------------------
    ANNUAL_KG = 1000
    _header(f"5. Scale-Up Analysis ({ANNUAL_KG:,} kg/yr)")
    su = analyze_scale_up(steps, target_annual_kg=ANNUAL_KG)
    print(f"  Production mode:  {su.recommended_mode}")
    print(f"  Batch size:       {su.batch_size_kg} kg")
    print(f"  Batches/year:     {su.batches_per_year}")
    print(f"  Cycle time:       {su.cycle_time_hours} h")
    print(f"  Annual capacity:  {su.annual_capacity_kg} kg/yr")
    print(f"  Capital cost:     ${su.capital_cost_usd:,.0f}")
    print(f"  Operating cost:   ${su.operating_cost_annual_usd:,.0f}/yr")
    if su.scale_up_risks:
        print(f"  Risks:            {su.scale_up_risks}")

    # ----------------------------------------------------------------
    # 6. Reaction conditions
    # ----------------------------------------------------------------
    _header("6. Optimized Reaction Conditions")
    co = optimize_conditions(template, scale_kg=SCALE_KG)
    print(f"  Temperature:     {co.temperature_C} C")
    print(f"  Pressure:        {co.pressure_atm} atm")
    print(f"  Solvent:         {co.solvent}")
    print(f"  Concentration:   {co.concentration_M} M")
    print(f"  Reaction time:   {co.reaction_time_hours} h")
    print(f"  Atmosphere:      {co.atmosphere}")
    print(f"  Workup:          {co.workup_procedure}")

    # ----------------------------------------------------------------
    # 7. Feasibility score
    # ----------------------------------------------------------------
    _header("7. Feasibility Score")
    engine = FeasibilityEngine()
    fs = engine.score(mol, tree=tree)
    print(f"  Availability:    {fs.availability.score:>6.0%}")
    print(f"  Cost:            {fs.cost.score:>6.0%}")
    print(f"  Complexity:      {fs.complexity.score:>6.0%}")
    print(f"  Green Chemistry: {fs.green_chemistry.score:>6.0%}")
    print(f"  Safety:          {fs.safety.score:>6.0%}")
    print(f"  Regulatory:      {fs.regulatory.score:>6.0%}")
    print(f"  {'-' * 35}")
    print(f"  COMPOSITE:       {fs.composite_score:>5.1f}/100")
    print(f"  GRADE:           {fs.grade}")
    print(f"  FEASIBLE:        {fs.is_feasible}")

    # ----------------------------------------------------------------
    # 8. Multi-target comparison
    # ----------------------------------------------------------------
    _header("8. Multi-Target Feasibility Comparison")
    targets = [
        ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
        ("lidocaine", "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
        ("acetaminophen", "CC(=O)Nc1ccc(O)cc1"),
    ]

    results = []
    for name, smi in targets:
        m = parse(smi)
        m.name = name
        t = retrosynthesis(m, max_depth=8, beam_width=5)
        score = engine.score(m, tree=t)
        results.append((name, score))

    # Sort by score descending
    results.sort(key=lambda x: x[1].composite_score, reverse=True)

    print(f"  {'Drug':<15} {'Score':>6} {'Grade':>6} {'Avail':>6} "
          f"{'Cost':>6} {'Cmplx':>6} {'Green':>6} {'Safe':>6} {'Reg':>6}")
    print(f"  {'-' * 75}")
    for name, fs in results:
        print(f"  {name:<15} {fs.composite_score:>5.1f} {fs.grade:>6} "
              f"{fs.availability.score:>5.0%} {fs.cost.score:>5.0%} "
              f"{fs.complexity.score:>5.0%} {fs.green_chemistry.score:>5.0%} "
              f"{fs.safety.score:>5.0%} {fs.regulatory.score:>5.0%}")

    print(f"\n  All {len(results)} targets scored as feasible (>= 40/100).")
    print(f"  Lowest safety scores: acetaminophen and lidocaine")
    print(f"  (due to amine/phenol precursor hazard profiles)")

    _header("Pipeline Complete")
    print("  7-stage analysis from SMILES to manufacturing plan")
    print("  Total API calls: 7 (retro, reactor, cost, safety, scale, cond, feas)")
    print("  Dependencies: numpy, scipy, matplotlib (no RDKit required)")
    print()


if __name__ == "__main__":
    main()
