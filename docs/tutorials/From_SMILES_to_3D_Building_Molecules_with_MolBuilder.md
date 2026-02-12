# From SMILES to 3D: Building Molecules with MolBuilder

**By Wedge_Dev | Published on Teapot Commons**

## What You'll Learn

Parse any SMILES string into a fully three-dimensional molecular structure, inspect its atoms and bonds programmatically, and export to standard chemistry file formats -- all in a few lines of Python.

## Prerequisites

- **Python** 3.11 or later
- **MolBuilder** 1.1.1 (`pip install molbuilder`)
- **MolBuilder Cloud SDK** 0.1.1 (`pip install molbuilder-client`) -- optional, for cloud API examples
- Basic familiarity with chemical structures (what atoms and bonds are)

---

## 1. A Quick SMILES Primer

SMILES (Simplified Molecular-Input Line-Entry System) encodes molecular structure as a compact string. Here are the essentials:

| SMILES | Molecule | Key feature |
|--------|----------|-------------|
| `C`    | Methane  | Single atom; hydrogens are implicit |
| `CCO`  | Ethanol  | Chain of atoms; bonds implied between adjacent letters |
| `C=O`  | Formaldehyde | `=` for double bond, `#` for triple |
| `c1ccccc1` | Benzene | Lowercase = aromatic; digits open/close rings |
| `Cn1cnc2c1c(=O)n(c(=O)n2C)C` | Caffeine | Branches in parentheses |

SMILES is the *lingua franca* of cheminformatics. If a molecule has a PubChem or ChEMBL entry, it has a SMILES string.

## 2. Parse and Inspect a Molecule

MolBuilder's SMILES parser reads a SMILES string, adds implicit hydrogens, assigns hybridization states, and places every atom in 3D space.

```python
from molbuilder.smiles.parser import parse

mol = parse("CCO")  # ethanol

print(f"Atoms: {len(mol.atoms)}")   # 9 (2 C + 1 O + 6 H)
print(f"Bonds: {len(mol.bonds)}")   # 8

for atom in mol.atoms:
    x, y, z = atom.position
    print(f"  {atom.index:2d}  {atom.symbol:2s}  ({x:7.3f}, {y:7.3f}, {z:7.3f})  {atom.hybridization}")
```

Output (abbreviated):

```
Atoms: 9
Bonds: 8
   0  C   (  0.000,   0.000,   0.000)  Hybridization.SP3
   1  C   (  0.000,   0.000,   1.540)  Hybridization.SP3
   2  O   ( -1.168,  -0.674,   2.017)  Hybridization.SP3
   ...
```

Every `Atom` object carries:

- `symbol` -- element symbol (`"C"`, `"O"`, `"H"`, ...)
- `position` -- NumPy array of 3D coordinates in Angstroms
- `index` -- integer index in the molecule
- `hybridization` -- `SP3`, `SP2`, or `SP`
- `formal_charge` -- integer (0 for neutral atoms)

Every `Bond` object carries:

- `atom_i`, `atom_j` -- indices of the bonded atoms
- `order` -- 1 (single), 2 (double), or 3 (triple)
- `rotatable` -- `True` for single bonds not in a ring

## 3. Geometry Queries

The `Molecule` object provides methods for common geometry measurements:

```python
# Bond length between atoms 0 and 1 (C-C)
print(f"C-C distance: {mol.distance(0, 1):.3f} A")  # ~1.540 A

# Bond angle: H-C-C
print(f"H-C-C angle: {mol.bond_angle(3, 0, 1):.1f} deg")  # ~109.5 deg

# Neighbors of atom 1
print(f"Neighbors of C1: {mol.neighbors(1)}")  # [0, 2, 6, 7]
```

## 4. Try Larger Molecules

The parser handles aromatic rings, branching, and complex heterocycles:

```python
caffeine = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
print(f"Caffeine: {len(caffeine.atoms)} atoms, {len(caffeine.bonds)} bonds")

aspirin = parse("CC(=O)Oc1ccccc1C(=O)O")
print(f"Aspirin:  {len(aspirin.atoms)} atoms, {len(aspirin.bonds)} bonds")
```

## 5. Export to Standard File Formats

MolBuilder writes three widely-used formats:

```python
from molbuilder.io.xyz import write_xyz
from molbuilder.io.mol_sdf import write_mol
from molbuilder.io.json_io import write_json

mol = parse("CCO")

write_xyz(mol, "ethanol.xyz")    # XYZ -- simple coordinate table
write_mol(mol, "ethanol.mol")    # MDL MOL -- bond connectivity included
write_json(mol, "ethanol.json")  # JSON -- full MolBuilder data model
```

**XYZ** is read by virtually every molecular viewer (Avogadro, VMD, PyMOL). **MOL/SDF** is the pharmaceutical industry standard. **JSON** round-trips all MolBuilder metadata including hybridization and chirality.

You can also read these files back:

```python
from molbuilder.io.xyz import read_xyz
from molbuilder.io.mol_sdf import read_mol
from molbuilder.io.json_io import read_json

mol2 = read_xyz("ethanol.xyz")
mol3 = read_mol("ethanol.mol")
mol4 = read_json("ethanol.json")
```

---

## Complete Example

```python
"""From SMILES to 3D -- complete working example.

Requirements: pip install molbuilder==1.1.1
"""
from molbuilder.smiles.parser import parse
from molbuilder.io.xyz import write_xyz
from molbuilder.io.mol_sdf import write_mol
from molbuilder.io.json_io import write_json

# --- Parse four molecules ---
molecules = {
    "methane":  "C",
    "ethanol":  "CCO",
    "caffeine": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "aspirin":  "CC(=O)Oc1ccccc1C(=O)O",
}

for name, smiles in molecules.items():
    mol = parse(smiles)
    mol.name = name
    print(f"\n{name} ({smiles})")
    print(f"  Atoms: {len(mol.atoms)}, Bonds: {len(mol.bonds)}")
    for atom in mol.atoms[:4]:
        x, y, z = atom.position
        print(f"    {atom.symbol:2s}  ({x:7.3f}, {y:7.3f}, {z:7.3f})")
    if len(mol.atoms) > 4:
        print(f"    ... and {len(mol.atoms) - 4} more atoms")

# --- Export ethanol to all three formats ---
ethanol = parse("CCO")
ethanol.name = "ethanol"

write_xyz(ethanol, "ethanol.xyz")
write_mol(ethanol, "ethanol.mol")
write_json(ethanol, "ethanol.json")

print("\nExported ethanol.xyz, ethanol.mol, ethanol.json")
```

---

## Using the Cloud API

If you prefer not to install MolBuilder locally, the same functionality is available through the cloud API:

```python
"""Cloud API equivalent -- requires pip install molbuilder-client==0.1.1"""
from molbuilder_client import MolBuilder

client = MolBuilder(
    api_key="mb_YOUR_KEY_HERE",
    base_url="https://molbuilder-api-production.up.railway.app",
)

with client:
    # Parse SMILES on the server
    mol_info = client.from_smiles("CCO", name="ethanol")
    print(f"ID: {mol_info.id}")
    print(f"Atoms: {mol_info.num_atoms}, Bonds: {mol_info.num_bonds}")

    # Get full 3D coordinates
    structure = client.get_3d(mol_info.id)
    for atom in structure.atoms:
        print(f"  {atom.symbol} {atom.position}")
```

Register for a free API key at the `/api/v1/auth/register` endpoint.

---

## Next Steps

- [Drug-Likeness Screening with Lipinski Rule-of-5](../Drug-Likeness%20Screening%20with%20Lipinski%20Rule%20of%20Five/) -- screen molecules for oral bioavailability
- [Retrosynthetic Analysis: Planning a Synthesis Route](../Retrosynthetic%20Analysis%20Planning%20a%20Synthesis%20Route/) -- plan how to make a molecule from purchasable starting materials
- [3D Visualization and Conformational Analysis](../3D%20Visualization%20and%20Conformational%20Analysis/) -- render molecules in 3D and explore conformations

---

## References

1. Weininger, D. SMILES, a Chemical Language and Information System. 1. Introduction to Methodology and Encoding Rules. *J. Chem. Inf. Comput. Sci.* **1988**, *28* (1), 31--36. DOI: 10.1021/ci00057a005
2. O'Boyle, N. M.; Banck, M.; James, C. A.; Morley, C.; Vandermeersch, T.; Hutchison, G. R. Open Babel: An Open Chemical Toolbox. *J. Cheminf.* **2011**, *3*, 33. DOI: 10.1186/1758-2946-3-33

---

*Tested with MolBuilder 1.1.1, molbuilder-client 0.1.1, Python 3.13. This article was written with AI assistance using Claude (Anthropic) and reviewed by the author.*
