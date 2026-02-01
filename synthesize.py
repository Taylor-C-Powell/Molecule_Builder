"""Quick synthesis planner - enter a SMILES string, get a full pathway."""
import sys

def main():
    print("=" * 60)
    print("  MOLECULE SYNTHESIS PLANNER")
    print("=" * 60)
    print()
    print("  Example SMILES strings:")
    print("    CCO           - ethanol")
    print("    CC(=O)O       - acetic acid")
    print("    c1ccccc1      - benzene")
    print("    CC(=O)Oc1ccccc1C(=O)O  - aspirin")
    print("    C1CCCCC1O     - cyclohexanol")
    print("    CC(O)=O       - acetic acid")
    print("    CCN           - ethylamine")
    print("    CC=O          - acetaldehyde")
    print("    CCCC          - butane")
    print("    CC(=O)CC      - 2-butanone (MEK)")
    print()

    smiles = input("  Enter SMILES string: ").strip()
    if not smiles:
        print("  No input. Exiting.")
        return

    name = input("  Molecule name (optional): ").strip()
    if not name:
        name = smiles

    scale_input = input("  Production scale in kg [default: 1.0]: ").strip()
    try:
        scale_kg = float(scale_input) if scale_input else 1.0
    except ValueError:
        scale_kg = 1.0

    print()
    print("Working...")
    print()

    # 1. Parse molecule
    from molbuilder.smiles import parse, to_smiles
    try:
        mol = parse(smiles)
        mol.name = name
    except Exception as e:
        print(f"  ERROR: Could not parse SMILES '{smiles}': {e}")
        return

    # 2. Molecule report
    from molbuilder.reports import generate_molecule_report
    print(generate_molecule_report(mol))

    # 3. Functional group detection
    from molbuilder.reactions.functional_group_detect import detect_functional_groups
    fgs = detect_functional_groups(mol)
    if fgs:
        print(f"  Detected functional groups: {', '.join(fg.name for fg in fgs)}")
    else:
        print("  No functional groups detected.")
    print()

    # 4. Retrosynthesis
    from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree, is_purchasable
    if is_purchasable(to_smiles(mol)):
        print(f"  '{name}' is directly purchasable as a starting material!")
        print("  No synthesis route needed.")
        return

    tree = retrosynthesis(mol, max_depth=5, beam_width=5)
    print(format_tree(tree))

    # 5. Forward synthesis route
    from molbuilder.reactions.synthesis_route import extract_best_route, format_route
    if tree.target.best_disconnection:
        route = extract_best_route(tree)
        print(format_route(route))

        # 6. Process engineering for each step
        from molbuilder.process.reactor import select_reactor
        from molbuilder.process.conditions import optimize_conditions
        from molbuilder.process.purification import recommend_purification
        from molbuilder.process.safety import assess_safety
        from molbuilder.process.costing import estimate_cost
        from molbuilder.process.scale_up import analyze_scale_up

        print("=" * 60)
        print("  PROCESS ENGINEERING DETAILS")
        print("=" * 60)
        for step in route.steps:
            t = step.template
            reactor = select_reactor(t, scale_kg)
            conditions = optimize_conditions(t, scale_kg)
            purif = recommend_purification(t, scale_kg)
            print(f"\n  Step {step.step_number}: {t.name}")
            print(f"    Reactor: {reactor.reactor_type.name} "
                  f"({reactor.volume_L:.1f} L, {reactor.material})")
            print(f"    Temperature: {conditions.temperature_C:.0f} C")
            print(f"    Solvent: {conditions.solvent}")
            print(f"    Atmosphere: {conditions.atmosphere}")
            print(f"    Reaction time: {conditions.reaction_time_hours:.1f} h")
            print(f"    Addition: {conditions.addition_rate}")
            print(f"    Workup: {conditions.workup_procedure}")
            if purif:
                print(f"    Purification:")
                for p in purif:
                    print(f"      - {p.method.name}: {p.description} "
                          f"(recovery ~{p.estimated_recovery:.0f}%, "
                          f"purity ~{p.estimated_purity:.0f}%)")

        # 7. Safety assessment
        assessments = assess_safety(route.steps)
        print()
        print("=" * 60)
        print("  SAFETY ASSESSMENT")
        print("=" * 60)
        for a in assessments:
            print(f"\n  Step {a.step_number}: {a.step_name}")
            print(f"    Risk level: {a.risk_level.upper()}")
            print(f"    PPE: {', '.join(a.ppe_required)}")
            if a.hazards:
                for h in a.hazards[:5]:
                    descs = h.hazard_descriptions[:2]
                    print(f"    {h.reagent_name}: {'; '.join(descs)}")
            if a.incompatible_materials:
                print(f"    Incompatible: {', '.join(a.incompatible_materials)}")

        # 8. Cost estimate
        cost = estimate_cost(route.steps, scale_kg)
        print()
        print("=" * 60)
        print(f"  COST ESTIMATE ({scale_kg:.1f} kg scale)")
        print("=" * 60)
        print(f"    Total cost:     ${cost.total_usd:>10,.2f}")
        print(f"    Cost per kg:    ${cost.per_kg_usd:>10,.2f}")
        print(f"    Raw materials:  ${cost.breakdown.raw_materials_usd:>10,.2f}")
        print(f"    Labor:          ${cost.breakdown.labor_usd:>10,.2f}")
        print(f"    Equipment:      ${cost.breakdown.equipment_usd:>10,.2f}")
        print(f"    Energy:         ${cost.breakdown.energy_usd:>10,.2f}")
        print(f"    Waste disposal: ${cost.breakdown.waste_disposal_usd:>10,.2f}")
        print(f"    Overhead:       ${cost.breakdown.overhead_usd:>10,.2f}")

        # 9. Scale-up analysis
        annual_kg = scale_kg * 100
        scaleup = analyze_scale_up(route.steps, annual_kg)
        print()
        print("=" * 60)
        print(f"  SCALE-UP ANALYSIS ({annual_kg:.0f} kg/year)")
        print("=" * 60)
        print(f"    Mode: {scaleup.recommended_mode}")
        if scaleup.batch_size_kg:
            print(f"    Batch size: {scaleup.batch_size_kg:.1f} kg")
        if scaleup.batches_per_year:
            print(f"    Batches/year: {scaleup.batches_per_year}")
        print(f"    Cycle time: {scaleup.cycle_time_hours:.1f} h")
        print(f"    Annual capacity: {scaleup.annual_capacity_kg:.0f} kg")
        print(f"    Capital cost: ${scaleup.capital_cost_usd:,.0f}")
        print(f"    Operating cost: ${scaleup.operating_cost_annual_usd:,.0f}/year")
        if scaleup.scale_up_risks:
            print(f"    Risks:")
            for r in scaleup.scale_up_risks:
                print(f"      - {r}")
    else:
        print("  No synthesis route found for this molecule.")
        print("  The molecule may be too simple or no matching reaction templates exist.")

    print()
    print("=" * 60)
    print("  DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
