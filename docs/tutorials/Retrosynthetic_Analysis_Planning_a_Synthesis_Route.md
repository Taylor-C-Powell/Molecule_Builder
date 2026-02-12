# Retrosynthetic Analysis: Planning a Synthesis Route

**By Wedge_Dev | Published on Teapot Commons**

## What You'll Learn

Use MolBuilder's retrosynthesis engine to work backward from a target molecule to purchasable starting materials: detect functional groups, score disconnection strategies, and generate a multi-step synthesis plan with a single function call.

## Prerequisites

- **Python** 3.11 or later
- **MolBuilder** 1.1.1 (`pip install molbuilder`)
- **MolBuilder Cloud SDK** 0.1.1 (`pip install molbuilder-client`) -- optional, for cloud API examples
- Familiarity with SMILES notation (see [From SMILES to 3D](../From%20SMILES%20to%203D%20Building%20Molecules%20with%20MolBuilder/))
- Basic organic chemistry (functional groups, reaction types)

---

## 1. What Is Retrosynthetic Analysis?

Retrosynthesis is the strategy of working backward from a target molecule: at each step, you identify a bond to "disconnect" and ask "what simpler precursors could form this product?" This process repeats until every precursor is commercially available.

E. J. Corey formalized this approach in the 1960s and received the Nobel Prize in Chemistry in 1990 for it. MolBuilder automates the key steps:

1. **Detect functional groups** in the target molecule
2. **Propose disconnections** using a library of named reaction templates
3. **Score** each disconnection by feasibility, selectivity, and precursor availability
4. **Recurse** on the precursor fragments until purchasable materials are reached

## 2. A Simple Example: Ethyl Acetate

Ethyl acetate is a common solvent and flavoring agent. Its SMILES is `CCOC(=O)C`.

```python
from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree

target = parse("CCOC(=O)C")
tree = retrosynthesis(target)

print(f"Target:       ethyl acetate")
print(f"Routes found: {tree.routes_found}")
print(f"Max depth:    {tree.max_depth}")
print(f"Beam width:   {tree.beam_width}")
```

### Examining the Disconnection Tree

```python
node = tree.target
print(f"\nFunctional groups: {[fg.name for fg in node.functional_groups]}")
print(f"Is purchasable:    {node.is_purchasable}")
print(f"\nTop disconnections:")
for disc in node.disconnections:
    print(f"  {disc.template.name:<35s} ({disc.template.category})  score={disc.score:.1f}")
```

Output:

```
Functional groups: ['ester']
Is purchasable:    False

Top disconnections:
  Grignard addition to ester          (COUPLING)  score=63.1
  Liebeskind-Srogl coupling           (COUPLING)  score=62.5
  Claisen condensation                (COUPLING)  score=60.6
  Tischenko reaction                  (CARBONYL)  score=52.8
  Baeyer-Villiger oxidation           (CARBONYL)  score=52.8
```

The engine identified the ester functional group, proposed five disconnection strategies, and ranked them by score. The Grignard route scores highest because its precursors are inexpensive and widely available.

### ASCII Tree Visualization

```python
print(format_tree(tree))
```

The `format_tree()` function renders the retrosynthetic tree as an indented ASCII diagram, showing each intermediate and whether it is purchasable.

## 3. Tuning the Search

The `retrosynthesis()` function takes two optional parameters:

| Parameter | Default | Effect |
|-----------|---------|--------|
| `max_depth` | 8 | Maximum number of retrosynthetic steps to explore |
| `beam_width` | 5 | Number of disconnections retained at each level |

Increasing `beam_width` explores more routes but takes longer. For complex targets, start with defaults and increase if no satisfactory route is found:

```python
# Broader search for a complex target
deep_tree = retrosynthesis(target, max_depth=10, beam_width=8)
print(f"Routes found with wider search: {deep_tree.routes_found}")
```

## 4. Working with Disconnection Details

Each disconnection carries a `ReactionTemplate` with conditions:

```python
from molbuilder.reactions.retrosynthesis import get_purchasable

best = tree.target.best_disconnection
print(f"Best route:     {best.template.name}")
print(f"Category:       {best.template.category}")
print(f"Score:          {best.score:.1f}")
print(f"Precursors:     {[p.smiles for p in best.precursors]}")
print(f"Purchasable?    {[p.name for p in best.precursors if get_purchasable(p.smiles)]}")
```

The score (0-100) incorporates:
- **Reaction reliability** -- well-known reactions score higher
- **Precursor availability** -- purchasable materials boost the score
- **Selectivity** -- reactions with fewer side products score higher

## 5. A More Complex Target: Ibuprofen-Like Molecule

```python
from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis

ibuprofen = parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
tree = retrosynthesis(ibuprofen, max_depth=8, beam_width=5)

print(f"Ibuprofen retrosynthesis:")
print(f"  Functional groups: {[fg.name for fg in tree.target.functional_groups]}")
print(f"  Routes found:      {tree.routes_found}")
print(f"  Best disconnection: {tree.target.best_disconnection.template.name}")
print(f"  Score:             {tree.target.best_disconnection.score:.1f}")
```

For a more complex molecule, the tree will have multiple levels, with intermediates that themselves require further disconnection.

---

## Complete Example

```python
"""Retrosynthetic Analysis -- complete working example.

Requirements: pip install molbuilder==1.1.1
"""
from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree

targets = {
    "ethyl acetate": "CCOC(=O)C",
    "phenol":        "c1ccc(cc1)O",
    "ibuprofen":     "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
}

for name, smiles in targets.items():
    mol = parse(smiles)
    tree = retrosynthesis(mol)
    node = tree.target

    print(f"\n{'='*60}")
    print(f"Target: {name} ({smiles})")
    print(f"Functional groups: {[fg.name for fg in node.functional_groups]}")
    print(f"Purchasable:       {node.is_purchasable}")
    print(f"Routes found:      {tree.routes_found}")

    if node.disconnections:
        print(f"Top 3 disconnections:")
        for disc in node.disconnections[:3]:
            print(f"  {disc.template.name:<35s} score={disc.score:.1f}")

    if node.best_disconnection:
        best = node.best_disconnection
        print(f"Recommended:       {best.template.name} (score {best.score:.1f})")
```

---

## Using the Cloud API

The retrosynthesis engine is available as a REST endpoint:

```python
"""Cloud API equivalent -- requires pip install molbuilder-client==0.1.1"""
from molbuilder_client import MolBuilder

client = MolBuilder(
    api_key="mb_YOUR_KEY_HERE",
    base_url="https://molbuilder-api-production.up.railway.app",
)

with client:
    plan = client.retrosynthesis(
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # ibuprofen
        max_depth=5,
        beam_width=5,
    )

    print(f"Routes found: {plan.routes_found}")
    if plan.best_route:
        print(f"Total steps:  {plan.best_route.total_steps}")
        print(f"Overall yield: {plan.best_route.overall_yield:.0%}")
        for step in plan.best_route.steps:
            print(f"  Step {step.step_number}: {step.reaction_name} "
                  f"(yield {step.expected_yield:.0%})")
```

The API response includes a `best_route` with step-by-step conditions, expected yields, and starting material costs.

---

## Next Steps

- [From SMILES to 3D](../From%20SMILES%20to%203D%20Building%20Molecules%20with%20MolBuilder/) -- parse the SMILES strings used in this tutorial
- [Drug-Likeness Screening](../Drug-Likeness%20Screening%20with%20Lipinski%20Rule%20of%20Five/) -- check if your target passes Lipinski screening before planning synthesis
- [Process Engineering: From Lab to Plant](../Process%20Engineering%20From%20Lab%20to%20Plant/) -- take your synthesis route and evaluate reactor selection, safety, and cost at scale

---

## References

1. Corey, E. J. The Logic of Chemical Synthesis: Multistep Synthesis of Complex Carbogenic Molecules (Nobel Lecture). *Angew. Chem. Int. Ed. Engl.* **1991**, *30* (5), 455--465. DOI: 10.1002/anie.199104553
2. Corey, E. J.; Cheng, X.-M. *The Logic of Chemical Synthesis*; Wiley: New York, 1989.
3. Szymkuc, S.; Gajewska, E. P.; Klucznik, T.; Molga, K.; Dittwald, P.; Startek, M.; Bajczyk, M.; Grzybowski, B. A. Computer-Assisted Synthetic Planning: The End of the Beginning. *Angew. Chem. Int. Ed.* **2016**, *55* (20), 5904--5937. DOI: 10.1002/anie.201506101

---

*Tested with MolBuilder 1.1.1, molbuilder-client 0.1.1, Python 3.13. This article was written with AI assistance using Claude (Anthropic) and reviewed by the author.*
