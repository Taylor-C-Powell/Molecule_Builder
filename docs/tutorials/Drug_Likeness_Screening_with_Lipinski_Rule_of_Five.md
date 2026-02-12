# Drug-Likeness Screening with Lipinski Rule-of-5

**By Wedge_Dev | Published on Teapot Commons**

## What You'll Learn

Screen a set of molecules for oral bioavailability using Lipinski's Rule-of-Five, interpret the results, and identify compounds that are likely to fail in early-stage drug development -- all from SMILES strings.

## Prerequisites

- **Python** 3.11 or later
- **MolBuilder** 1.1.1 (`pip install molbuilder`)
- **MolBuilder Cloud SDK** 0.1.1 (`pip install molbuilder-client`) -- optional, for cloud API examples
- Familiarity with SMILES notation (see [From SMILES to 3D](../From%20SMILES%20to%203D%20Building%20Molecules%20with%20MolBuilder/))

---

## 1. What Is the Rule-of-Five?

In 1997, Christopher Lipinski and colleagues at Pfizer analyzed compounds that had reached Phase II clinical trials and found that poor oral absorption is associated with molecules that violate two or more of these thresholds:

| Property | Threshold | Why it matters |
|----------|-----------|----------------|
| Molecular weight | <= 500 Da | Larger molecules struggle to cross cell membranes |
| Calculated logP  | <= 5      | Too lipophilic = poor solubility |
| H-bond donors    | <= 5      | Too many donors reduce permeability |
| H-bond acceptors | <= 10     | Too many acceptors reduce permeability |

The name "Rule-of-Five" comes from the thresholds all being multiples of five. A molecule **passes** Lipinski screening if it has at most one violation.

## 2. Computing Properties for a Single Molecule

MolBuilder's `lipinski_properties()` function returns all Ro5 descriptors plus additional drug-likeness indicators:

```python
from molbuilder.smiles.parser import parse
from molbuilder.molecule.properties import lipinski_properties

aspirin = parse("CC(=O)Oc1ccccc1C(=O)O")
props = lipinski_properties(aspirin)

print(f"Molecular weight:    {props.molecular_weight:.1f} Da")
print(f"Calculated logP:     {props.logp:.2f}")
print(f"H-bond donors:       {props.hbd}")
print(f"H-bond acceptors:    {props.hba}")
print(f"Rotatable bonds:     {props.rotatable_bonds}")
print(f"TPSA:                {props.tpsa:.1f} A^2")
print(f"Heavy atom count:    {props.heavy_atom_count}")
print(f"Lipinski violations: {props.lipinski_violations}")
print(f"Passes Ro5:          {props.lipinski_pass}")
```

Output:

```
Molecular weight:    180.2 Da
Calculated logP:     -0.05
H-bond donors:       1
H-bond acceptors:    4
Rotatable bonds:     3
TPSA:                63.6 A^2
Heavy atom count:    13
Lipinski violations: 0
Passes Ro5:          True
```

Aspirin passes easily -- it is a small, well-behaved oral drug.

## 3. Understanding Each Property

**Molecular weight** -- Sum of atomic masses for all atoms (including hydrogens). Computed from the full 3D structure, so implicit hydrogens are included.

**Calculated logP** -- Uses the Wildman-Crippen atom-contribution method. Each heavy atom contributes a fragment-based logP increment; the total estimates the octanol-water partition coefficient. Positive values indicate lipophilicity.

**H-bond donors (HBD)** -- Count of N-H and O-H groups. These form hydrogen bonds with biological targets but also with water, which can limit membrane permeation.

**H-bond acceptors (HBA)** -- Count of nitrogen and oxygen atoms. The same trade-off applies: useful for target binding, costly for permeability.

**Rotatable bonds** -- Single bonds not in a ring. High flexibility increases entropy loss upon binding and correlates with poor oral bioavailability. Veber's guideline: <= 10 rotatable bonds.

**TPSA (Topological Polar Surface Area)** -- Sum of surface contributions from polar atoms (N, O, and their hydrogens), calculated by the Ertl fragment method. TPSA < 140 A^2 is a common oral absorption threshold.

## 4. Screen a Compound Library

```python
from molbuilder.smiles.parser import parse
from molbuilder.molecule.properties import lipinski_properties

library = {
    "Aspirin":     "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen":   "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Caffeine":    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "Metformin":   "CN(C)C(=N)NC(=N)N",
    "Cyclosporine A (violator)": "CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC"
        "(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1C)C(C)CC=CC)C)CC(C)C)C)"
        "CC(C)C)C)C)CC(C)C)C(C)C)CC(C)C)C)C(C)C)CC(C)C)C)C",
}

print(f"{'Compound':<30s} {'MW':>6s} {'logP':>6s} {'HBD':>4s} {'HBA':>4s} "
      f"{'Rot':>4s} {'TPSA':>6s} {'Viol':>4s} {'Pass':>5s}")
print("-" * 78)

for name, smiles in library.items():
    mol = parse(smiles)
    p = lipinski_properties(mol)
    print(f"{name:<30s} {p.molecular_weight:6.1f} {p.logp:6.2f} {p.hbd:4d} "
          f"{p.hba:4d} {p.rotatable_bonds:4d} {p.tpsa:6.1f} "
          f"{p.lipinski_violations:4d} {'Yes' if p.lipinski_pass else 'NO':>5s}")
```

Expected output:

```
Compound                           MW   logP  HBD  HBA  Rot   TPSA Viol  Pass
------------------------------------------------------------------------------
Aspirin                          180.2  -0.05    1    4    3   63.6    0   Yes
Ibuprofen                        206.3   2.97    1    2    4   37.3    0   Yes
Caffeine                         194.2  -1.03    0    6    0   58.4    0   Yes
Metformin                        129.2  -2.64    3    5    2  115.0    0   Yes
Cyclosporine A (violator)       1202.6  13.92    5   12    8  278.8    3    NO
```

Cyclosporine A violates three rules (MW, logP, and HBA) -- yet it *is* an approved oral drug. This illustrates an important caveat: the Rule-of-Five is a guideline, not an absolute filter. Cyclosporine achieves oral absorption through active transport mechanisms that bypass passive diffusion.

## 5. Interpreting Results

When screening hits:

- **0 violations**: Strong candidate for oral delivery. Proceed to further optimization.
- **1 violation**: Still likely oral. Most approved drugs have 0-1 violations.
- **2+ violations**: Flag for review. Consider alternative delivery routes (IV, topical, inhaled) or structural modifications to improve drug-likeness.

Common fixes for violations:

| Violation | Strategy |
|-----------|----------|
| MW too high | Remove non-essential substituents; use bioisosteric replacements |
| logP too high | Add polar groups (OH, NH2); replace phenyl with pyridine |
| Too many HBD | Methylate NH groups; replace OH with F |
| Too many HBA | Replace ether O with CH2; reduce ring heteroatoms |

---

## Complete Example

```python
"""Drug-likeness screening -- complete working example.

Requirements: pip install molbuilder==1.1.1
"""
from molbuilder.smiles.parser import parse
from molbuilder.molecule.properties import lipinski_properties

DRUG_LIBRARY = {
    "Aspirin":    "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen":  "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Caffeine":   "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "Metformin":  "CN(C)C(=N)NC(=N)N",
}

passes = 0
for name, smiles in DRUG_LIBRARY.items():
    mol = parse(smiles)
    props = lipinski_properties(mol)
    status = "PASS" if props.lipinski_pass else "FAIL"
    print(f"[{status}] {name}: MW={props.molecular_weight:.1f}, "
          f"logP={props.logp:.2f}, HBD={props.hbd}, HBA={props.hba}, "
          f"violations={props.lipinski_violations}")
    if props.lipinski_pass:
        passes += 1

print(f"\n{passes}/{len(DRUG_LIBRARY)} compounds pass Lipinski screening.")
```

---

## Using the Cloud API

The same screening is available without a local install:

```python
"""Cloud API equivalent -- requires pip install molbuilder-client==0.1.1"""
from molbuilder_client import MolBuilder

client = MolBuilder(
    api_key="mb_YOUR_KEY_HERE",
    base_url="https://molbuilder-api-production.up.railway.app",
)

with client:
    mol_info = client.from_smiles("CC(=O)Oc1ccccc1C(=O)O", name="aspirin")
    props = client.get_properties(mol_info.id)

    print(f"Molecular weight:    {props.molecular_weight:.1f} Da")
    print(f"Calculated logP:     {props.logp:.2f}")
    print(f"H-bond donors:       {props.hbd}")
    print(f"H-bond acceptors:    {props.hba}")
    print(f"Lipinski violations: {props.lipinski_violations}")
    print(f"Passes Ro5:          {props.lipinski_pass}")
```

The API response includes all the same Lipinski properties plus `formula`, `functional_groups`, and `heavy_atom_count`.

---

## Next Steps

- [From SMILES to 3D](../From%20SMILES%20to%203D%20Building%20Molecules%20with%20MolBuilder/) -- learn SMILES parsing and file export
- [Retrosynthetic Analysis](../Retrosynthetic%20Analysis%20Planning%20a%20Synthesis%20Route/) -- plan how to synthesize your drug candidates
- [Process Engineering](../Process%20Engineering%20From%20Lab%20to%20Plant/) -- evaluate manufacturing cost and safety at scale

---

## References

1. Lipinski, C. A.; Lombardo, F.; Dominy, B. W.; Feeney, P. J. Experimental and Computational Approaches to Estimate Solubility and Permeability in Drug Discovery and Development Settings. *Adv. Drug Delivery Rev.* **1997**, *23* (1--3), 3--25. DOI: 10.1016/S0169-409X(96)00423-1
2. Wildman, S. A.; Crippen, G. M. Prediction of Physicochemical Parameters by Atomic Contributions. *J. Chem. Inf. Comput. Sci.* **1999**, *39* (5), 868--873. DOI: 10.1021/ci990307l
3. Ertl, P.; Rohde, B.; Selzer, P. Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-Based Contributions and Its Application to the Prediction of Drug Transport Properties. *J. Med. Chem.* **2000**, *43* (20), 3714--3717. DOI: 10.1021/jm000942e
4. Veber, D. F.; Johnson, S. R.; Cheng, H.-Y.; Smith, B. R.; Ward, K. W.; Kopple, K. D. Molecular Properties That Influence the Oral Bioavailability of Drug Candidates. *J. Med. Chem.* **2002**, *45* (12), 2615--2623. DOI: 10.1021/jm020017n

---

*Tested with MolBuilder 1.1.1, molbuilder-client 0.1.1, Python 3.13. This article was written with AI assistance using Claude (Anthropic) and reviewed by the author.*
