# Process Engineering: From Lab to Plant

**By Wedge_Dev | Published on Teapot Commons**

## What You'll Learn

Take a synthesis route and evaluate it for manufacturing: select the right reactor, assess chemical safety hazards, estimate production costs, and analyze scale-up to 100 kg -- all with MolBuilder's process engineering module.

## Prerequisites

- **Python** 3.11 or later
- **MolBuilder** 1.1.1 (`pip install molbuilder`)
- **MolBuilder Cloud SDK** 0.1.1 (`pip install molbuilder-client`) -- optional, for cloud API examples
- Familiarity with retrosynthetic analysis (see [Retrosynthetic Analysis](../Retrosynthetic%20Analysis%20Planning%20a%20Synthesis%20Route/))

---

## 1. From Route to Process

In the retrosynthesis tutorial, we generated a synthesis plan from a target molecule. Process engineering answers the next question: *how do we actually make this at scale?*

MolBuilder's process module evaluates four dimensions:

| Module | Question | Key output |
|--------|----------|------------|
| `reactor` | What equipment do we need? | Reactor type, volume, materials, cost |
| `safety` | What are the hazards? | GHS hazards, PPE, emergency procedures |
| `costing` | How much will it cost? | Per-kg cost breakdown |
| `scale_up` | Can we scale to production? | Batch size, capital cost, risks |

## 2. Generating a Synthesis Route

First, we need a retrosynthetic route to evaluate:

```python
from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis

target = parse("CCOC(=O)C")  # ethyl acetate
tree = retrosynthesis(target)

# Extract the best route's steps for process evaluation
# Each node in the tree has a best_disconnection with template and precursors
node = tree.target
best = node.best_disconnection
print(f"Best route: {best.template.name}")
print(f"Category:   {best.template.category}")
```

## 3. Reactor Selection

The `select_reactor()` function recommends equipment based on the reaction type and production scale:

```python
from molbuilder.process.reactor import select_reactor

reactor = select_reactor(best.template, scale_kg=100.0)

print(f"Reactor type:      {reactor.reactor_type}")
print(f"Volume:            {reactor.volume_L:.0f} L")
print(f"Temperature:       {reactor.temperature_C:.0f} C")
print(f"Pressure:          {reactor.pressure_atm:.1f} atm")
print(f"Residence time:    {reactor.residence_time_min:.0f} min")
print(f"Mixing:            {reactor.mixing_type}")
print(f"Heat transfer:     {reactor.heat_transfer}")
print(f"Material:          {reactor.material}")
print(f"Estimated cost:    ${reactor.estimated_cost_usd:,.0f}")
```

The function considers:
- **Reaction temperature** -- cryogenic reactions need specialized vessels
- **Scale** -- microreactors for < 1 kg, batch reactors for 1-1000 kg, continuous flow for > 1000 kg
- **Exothermicity** -- highly exothermic reactions need better heat removal (jacketed or coil)
- **Catalyst type** -- heterogeneous catalysts may need a fixed-bed reactor

### Reactor Types

| Type | Typical scale | Best for |
|------|--------------|----------|
| BATCH | 1-1000 kg | Most synthetic chemistry |
| SEMI_BATCH | 10-500 kg | Controlled reagent addition |
| CSTR | 100+ kg/day | Continuous homogeneous reactions |
| PFR | 100+ kg/day | Fast reactions, good heat transfer |
| MICROREACTOR | 0.1-10 kg | Hazardous or highly exothermic |
| FIXED_BED | 100+ kg/day | Heterogeneous catalysis |

## 4. Safety Assessment

Chemical safety is non-negotiable. The `assess_safety()` function evaluates each synthesis step for GHS hazards, required PPE, and emergency procedures:

```python
from molbuilder.process.safety import assess_safety

# Create a step-like object with template and precursors
class Step:
    def __init__(self, template, precursors):
        self.template = template
        self.precursors = precursors

steps = [Step(best.template, best.precursors)]
safety_reports = assess_safety(steps)

for report in safety_reports:
    print(f"\nStep {report.step_number}: {report.step_name}")
    print(f"  Risk level: {report.risk_level}")
    print(f"  PPE required: {', '.join(report.ppe_required)}")
    print(f"  Engineering controls: {', '.join(report.engineering_controls)}")
    print(f"  Waste classification: {report.waste_classification}")

    for hazard in report.hazards:
        print(f"  Hazard: {hazard.reagent_name}")
        for desc in hazard.hazard_descriptions:
            print(f"    - {desc}")
```

Always review safety assessments with a qualified chemist before scaling up. MolBuilder's assessments are computational estimates, not substitutes for MSDS review and institutional safety protocols.

## 5. Cost Estimation

The `estimate_cost()` function breaks down manufacturing costs:

```python
from molbuilder.process.costing import estimate_cost

cost = estimate_cost(steps, scale_kg=100.0)

print(f"Total cost:        ${cost.total_usd:,.0f}")
print(f"Cost per kg:       ${cost.per_kg_usd:,.0f}/kg")
print(f"\nBreakdown:")
print(f"  Raw materials:   ${cost.breakdown.raw_materials_usd:,.0f}")
print(f"  Labor:           ${cost.breakdown.labor_usd:,.0f}")
print(f"  Equipment:       ${cost.breakdown.equipment_usd:,.0f}")
print(f"  Energy:          ${cost.breakdown.energy_usd:,.0f}")
print(f"  Waste disposal:  ${cost.breakdown.waste_disposal_usd:,.0f}")
print(f"  Overhead:        ${cost.breakdown.overhead_usd:,.0f}")

for note in cost.notes:
    print(f"  Note: {note}")
```

Cost drivers to watch:
- **Raw materials** dominate for simple chemistry; **labor** dominates for complex multi-step routes
- **Waste disposal** costs increase sharply with halogenated solvents
- **Equipment** costs are amortized across batches (check the notes for amortization assumptions)

## 6. Scale-Up Analysis

The `analyze_scale_up()` function evaluates production feasibility at a target annual volume:

```python
from molbuilder.process.scale_up import analyze_scale_up

scale = analyze_scale_up(steps, target_annual_kg=1000.0)

print(f"Target:            {scale.target_annual_kg:.0f} kg/year")
print(f"Recommended mode:  {scale.recommended_mode}")
print(f"Batch size:        {scale.batch_size_kg:.0f} kg")
print(f"Batches/year:      {scale.batches_per_year}")
print(f"Cycle time:        {scale.cycle_time_hours:.0f} hours")
print(f"Annual capacity:   {scale.annual_capacity_kg:.0f} kg/year")
print(f"Capital cost:      ${scale.capital_cost_usd:,.0f}")
print(f"Operating cost:    ${scale.operating_cost_annual_usd:,.0f}/year")

print(f"\nScale-up risks:")
for risk in scale.scale_up_risks:
    print(f"  - {risk}")

print(f"\nRecommendations:")
for rec in scale.recommendations:
    print(f"  - {rec}")
```

Key scale-up considerations:
- **Mixing efficiency** -- what works in a 1 L flask may not mix properly in a 500 L vessel
- **Heat removal** -- surface-to-volume ratio drops as vessels get larger
- **Mass transfer** -- gas-liquid reactions are particularly sensitive to scale

---

## Complete Example

```python
"""Process Engineering: From Lab to Plant -- complete working example.

Requirements: pip install molbuilder==1.1.1
"""
from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis
from molbuilder.process.reactor import select_reactor
from molbuilder.process.safety import assess_safety
from molbuilder.process.costing import estimate_cost
from molbuilder.process.scale_up import analyze_scale_up

# --- 1. Get a synthesis route ---
target = parse("CCOC(=O)C")  # ethyl acetate
tree = retrosynthesis(target)
best = tree.target.best_disconnection

print(f"Route: {best.template.name} (score {best.score:.1f})")
print(f"Precursors: {[p.smiles for p in best.precursors]}")

# --- 2. Reactor selection at 100 kg ---
reactor = select_reactor(best.template, scale_kg=100.0)
print(f"\nReactor: {reactor.reactor_type}, {reactor.volume_L:.0f} L, "
      f"${reactor.estimated_cost_usd:,.0f}")

# --- 3. Safety assessment ---
class Step:
    def __init__(self, template, precursors):
        self.template = template
        self.precursors = precursors

steps = [Step(best.template, best.precursors)]
safety = assess_safety(steps)
for s in safety:
    print(f"\nSafety [{s.risk_level}]: {s.step_name}")
    print(f"  PPE: {', '.join(s.ppe_required)}")

# --- 4. Cost estimate ---
cost = estimate_cost(steps, scale_kg=100.0)
print(f"\nCost: ${cost.total_usd:,.0f} total (${cost.per_kg_usd:,.0f}/kg)")

# --- 5. Scale-up analysis ---
scale = analyze_scale_up(steps, target_annual_kg=1000.0)
print(f"\nScale-up: {scale.recommended_mode}, {scale.batch_size_kg:.0f} kg batches")
print(f"  Capital: ${scale.capital_cost_usd:,.0f}")
print(f"  Operating: ${scale.operating_cost_annual_usd:,.0f}/year")
```

---

## Using the Cloud API

The full process evaluation pipeline is available as a single endpoint:

```python
"""Cloud API equivalent -- requires pip install molbuilder-client==0.1.1"""
from molbuilder_client import MolBuilder

client = MolBuilder(
    api_key="mb_YOUR_KEY_HERE",
    base_url="https://molbuilder-api-production.up.railway.app",
)

with client:
    result = client.process_evaluate(
        "CCOC(=O)C",  # ethyl acetate
        scale_kg=100.0,
        max_depth=5,
        beam_width=5,
    )

    print(f"Route found:    {result.route_found}")
    print(f"Total steps:    {result.total_steps}")
    print(f"Overall yield:  {result.overall_yield:.0%}")
    print(f"Cost per kg:    ${result.cost.per_kg_usd:,.0f}")
    print(f"Scale-up mode:  {result.scale_up.recommended_mode}")

    for step in result.step_details:
        print(f"\n  Step {step.step_number}: {step.reaction_name}")
        print(f"    Reactor: {step.reactor.reactor_type}, "
              f"{step.reactor.volume_L:.0f} L")

    for s in result.safety:
        print(f"\n  Safety [{s.risk_level}]: {s.step_name}")
        print(f"    PPE: {', '.join(s.ppe_required)}")
```

The API combines retrosynthesis, reactor selection, safety, costing, and scale-up into a single request -- useful for batch screening of candidate molecules.

---

## Next Steps

- [Retrosynthetic Analysis](../Retrosynthetic%20Analysis%20Planning%20a%20Synthesis%20Route/) -- understand how the synthesis route was generated
- [Drug-Likeness Screening](../Drug-Likeness%20Screening%20with%20Lipinski%20Rule%20of%20Five/) -- screen target molecules before committing to process design
- [From SMILES to 3D](../From%20SMILES%20to%203D%20Building%20Molecules%20with%20MolBuilder/) -- learn the fundamentals of molecular representation

---

## References

1. Silla Santos, M. H. Process Intensification: Transforming Chemical Engineering. *Chem. Eng. Prog.* **2000**, *96* (1), 22--34.
2. Levenspiel, O. *Chemical Reaction Engineering*, 3rd ed.; Wiley: New York, 1999.
3. Etchells, J. C. Why Reactions Run Away. *Org. Process Res. Dev.* **1997**, *1* (6), 435--437. DOI: 10.1021/op970027w

---

*Tested with MolBuilder 1.1.1, molbuilder-client 0.1.1, Python 3.13. This article was written with AI assistance using Claude (Anthropic) and reviewed by the author.*
