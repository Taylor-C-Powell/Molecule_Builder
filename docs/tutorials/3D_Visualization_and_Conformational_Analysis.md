# 3D Visualization and Conformational Analysis

**By Wedge_Dev | Published on Teapot Commons**

## What You'll Learn

Render molecules as interactive 3D ball-and-stick models, measure dihedral angles, scan torsional energy profiles, and compare conformations of ethane, butane, and cyclohexane -- all with MolBuilder and matplotlib.

## Prerequisites

- **Python** 3.11 or later
- **MolBuilder** 1.1.1 (`pip install molbuilder`)
- **matplotlib** (`pip install matplotlib`) -- included as a MolBuilder dependency
- **NumPy** (`pip install numpy`) -- included as a MolBuilder dependency
- Familiarity with SMILES notation (see [From SMILES to 3D](../From%20SMILES%20to%203D%20Building%20Molecules%20with%20MolBuilder/))

> **Note**: This tutorial covers local Python visualization only. The cloud API returns raw 3D coordinates but does not render images server-side.

---

## 1. Ball-and-Stick Rendering

MolBuilder includes a matplotlib-based 3D renderer that displays atoms as colored spheres and bonds as cylinders:

```python
from molbuilder.visualization.molecule_viz import visualize_molecule
from molbuilder.bonding.vsepr import VSEPRMolecule

# Build a VSEPR-geometry molecule for visualization
water = VSEPRMolecule("H2O")
visualize_molecule(water, show_lone_pairs=True, show_labels=True, show_angles=True)
```

The `visualize_molecule()` function renders:
- Atoms colored by element (CPK convention: C = dark gray, O = red, N = blue, H = white)
- Bond sticks connecting bonded atoms
- Optional lone-pair lobes (for VSEPR geometries)
- Optional atom labels and bond angle arcs

For comparing multiple geometries side-by-side:

```python
from molbuilder.visualization.molecule_viz import visualize_gallery
from molbuilder.bonding.vsepr import VSEPRMolecule

gallery = [
    VSEPRMolecule("CH4"),
    VSEPRMolecule("H2O"),
    VSEPRMolecule("NH3"),
    VSEPRMolecule("BF3"),
]
visualize_gallery(gallery, cols=2)
```

## 2. Dihedral Angles

A dihedral (torsion) angle describes rotation around a central bond. It is defined by four consecutive atoms i-j-k-l, where j-k is the bond being rotated:

```python
from molbuilder.smiles.parser import parse

butane = parse("CCCC")  # n-butane

# Measure the C0-C1-C2-C3 dihedral
angle = butane.dihedral_angle(0, 1, 2, 3)
print(f"Initial dihedral: {angle:.1f} deg")
```

The dihedral convention:
- **0 deg** = eclipsed (substituents aligned, highest energy)
- **60 deg** = gauche (partially staggered)
- **180 deg** = anti (fully staggered, lowest energy for butane)

You can set a specific dihedral angle:

```python
# Set to anti conformation
butane.set_dihedral(0, 1, 2, 3, 180.0)
print(f"Anti dihedral: {butane.dihedral_angle(0, 1, 2, 3):.1f} deg")

# Set to gauche conformation
butane.set_dihedral(0, 1, 2, 3, 60.0)
print(f"Gauche dihedral: {butane.dihedral_angle(0, 1, 2, 3):.1f} deg")
```

## 3. Torsional Energy Scan

MolBuilder can compute the OPLS-AA torsional strain energy for any rotatable bond. Combined with `rotate_dihedral()`, you can scan the full 360-degree profile:

```python
from molbuilder.smiles.parser import parse
import matplotlib.pyplot as plt

butane = parse("CCCC")
angles = []
energies = []

for deg in range(0, 360, 5):
    butane.set_dihedral(0, 1, 2, 3, float(deg))
    strain = butane.torsional_energy(1, 2)
    angles.append(deg)
    energies.append(strain.total_kj_per_mol)

plt.figure(figsize=(8, 4))
plt.plot(angles, energies, "b-", linewidth=2)
plt.xlabel("Dihedral angle (degrees)")
plt.ylabel("Torsional energy (kJ/mol)")
plt.title("Butane C1-C2 torsional energy profile")
plt.axvline(60, color="gray", linestyle="--", alpha=0.5, label="gauche")
plt.axvline(180, color="gray", linestyle=":", alpha=0.5, label="anti")
plt.axvline(0, color="red", linestyle="--", alpha=0.5, label="eclipsed")
plt.legend()
plt.tight_layout()
plt.savefig("butane_torsion_scan.png", dpi=150)
plt.show()
```

The resulting plot shows:
- **Minima** at 60 deg (gauche) and 180 deg (anti) -- staggered conformations
- **Maxima** at 0 deg and 120 deg -- eclipsed conformations
- The anti minimum is the global minimum (lowest energy)
- The energy difference between anti and gauche is the "gauche penalty" (~3.8 kJ/mol)

## 4. Newman Projections

Newman projections are the classic way to visualize dihedral relationships. MolBuilder generates projection data for any bond:

```python
from molbuilder.smiles.parser import parse

ethane = parse("CC")
newman = ethane.newman_projection(0, 1)

print(f"Front atom: {newman.front_atom}, Back atom: {newman.back_atom}")
print(f"Dihedral: {newman.dihedral_deg:.1f} deg")
print("Front substituents:")
for idx, symbol, angle in newman.front_substituents:
    print(f"  atom {idx} ({symbol}) at {angle:.1f} deg")
print("Back substituents:")
for idx, symbol, angle in newman.back_substituents:
    print(f"  atom {idx} ({symbol}) at {angle:.1f} deg")
```

The `NewmanProjection` dataclass gives you:
- `front_atom` / `back_atom` -- the two atoms defining the bond being viewed
- `front_substituents` / `back_substituents` -- lists of `(atom_index, symbol, angle_degrees)` for each substituent, measured clockwise from 12 o'clock
- `dihedral_deg` -- the overall dihedral between the first front and first back substituent

## 5. Chair vs. Boat Cyclohexane

Cyclohexane is the classic example of conformational analysis. The chair form is ~29 kJ/mol more stable than the boat:

```python
from molbuilder.smiles.parser import parse

cyclohexane = parse("C1CCCCC1")

# Measure ring torsions
for i in range(6):
    j = (i + 1) % 6
    k = (i + 2) % 6
    l = (i + 3) % 6
    angle = cyclohexane.dihedral_angle(i, j, k, l)
    print(f"  C{i}-C{j}-C{k}-C{l}: {angle:7.1f} deg")
```

In a perfect chair, the ring dihedrals alternate between approximately +55 deg and -55 deg. Deviations indicate a twist-boat or other intermediate conformation.

---

## Complete Example

```python
"""3D Visualization and Conformational Analysis -- complete working example.

Requirements: pip install molbuilder==1.1.1 matplotlib
"""
from molbuilder.smiles.parser import parse
import matplotlib.pyplot as plt

# --- Ethane: staggered vs eclipsed ---
ethane = parse("CC")
staggered = ethane.torsional_energy(0, 1)
print(f"Ethane staggered energy: {staggered.total_kj_per_mol:.2f} kJ/mol")

ethane.set_dihedral(2, 0, 1, 5, 0.0)  # Force eclipsed
eclipsed = ethane.torsional_energy(0, 1)
print(f"Ethane eclipsed energy:  {eclipsed.total_kj_per_mol:.2f} kJ/mol")

# --- Butane: full torsion scan ---
butane = parse("CCCC")
angles, energies = [], []
for deg in range(0, 360, 5):
    butane.set_dihedral(0, 1, 2, 3, float(deg))
    strain = butane.torsional_energy(1, 2)
    angles.append(deg)
    energies.append(strain.total_kj_per_mol)

plt.figure(figsize=(8, 4))
plt.plot(angles, energies, "b-", linewidth=2)
plt.xlabel("Dihedral angle (degrees)")
plt.ylabel("Torsional energy (kJ/mol)")
plt.title("Butane C1-C2 torsional energy profile")
plt.tight_layout()
plt.savefig("butane_torsion_scan.png", dpi=150)
print("Saved butane_torsion_scan.png")

# --- Newman projection ---
newman = butane.newman_projection(1, 2)
print(f"\nNewman projection of butane C1-C2:")
print(f"  Dihedral: {newman.dihedral_deg:.1f} deg")
print(f"  Front substituents: {len(newman.front_substituents)}")
print(f"  Back substituents:  {len(newman.back_substituents)}")

# --- Cyclohexane ring torsions ---
cyclohexane = parse("C1CCCCC1")
print("\nCyclohexane ring dihedrals:")
for i in range(6):
    j, k, l = (i+1) % 6, (i+2) % 6, (i+3) % 6
    angle = cyclohexane.dihedral_angle(i, j, k, l)
    print(f"  C{i}-C{j}-C{k}-C{l}: {angle:7.1f} deg")
```

---

## Next Steps

- [From SMILES to 3D](../From%20SMILES%20to%203D%20Building%20Molecules%20with%20MolBuilder/) -- learn SMILES parsing and file I/O
- [Drug-Likeness Screening](../Drug-Likeness%20Screening%20with%20Lipinski%20Rule%20of%20Five/) -- screen your conformations for drug-likeness
- [Quantum Mechanics: Atomic Structure and Orbitals](../Quantum%20Mechanics%20Atomic%20Structure%20and%20Orbitals/) -- dive deeper into the electronic structure of atoms

---

## References

1. Allinger, N. L. Conformational Analysis. 130. MM2. A Hydrocarbon Force Field Utilizing V1 and V2 Torsional Terms. *J. Am. Chem. Soc.* **1977**, *99* (25), 8127--8134. DOI: 10.1021/ja00467a001
2. Jorgensen, W. L.; Maxwell, D. S.; Tirado-Rives, J. Development and Testing of the OPLS All-Atom Force Field on Conformational Energetics and Properties of Organic Liquids. *J. Am. Chem. Soc.* **1996**, *118* (45), 11225--11236. DOI: 10.1021/ja9621760
3. Barton, D. H. R. The Conformation of the Steroid Nucleus. *Experientia* **1950**, *6*, 316--320. DOI: 10.1007/BF02170915

---

*Tested with MolBuilder 1.1.1, Python 3.13, matplotlib 3.9. This article was written with AI assistance using Claude (Anthropic) and reviewed by the author.*
