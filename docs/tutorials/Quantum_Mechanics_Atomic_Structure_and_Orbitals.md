# Quantum Mechanics: Atomic Structure and Orbitals

**By Wedge_Dev | Published on Teapot Commons**

## What You'll Learn

Explore the quantum mechanical structure of atoms using MolBuilder: compute Bohr model energy levels, predict spectral emission lines, generate electron configurations with Aufbau filling, and visualize atomic orbitals in 3D.

## Prerequisites

- **Python** 3.11 or later
- **MolBuilder** 1.1.1 (`pip install molbuilder`)
- **matplotlib** (`pip install matplotlib`) -- included as a MolBuilder dependency
- Introductory chemistry (what electrons, orbitals, and quantum numbers are)

> **Note**: This tutorial covers local Python computation and visualization only. The cloud API does not currently expose quantum mechanics endpoints.

---

## 1. The Bohr Model

Niels Bohr's 1913 model treats electrons as orbiting the nucleus at discrete energy levels. While superseded by quantum mechanics for accurate calculations, the Bohr model still correctly predicts hydrogen's spectral lines and provides an excellent pedagogical starting point.

```python
from molbuilder.atomic.bohr import BohrAtom

hydrogen = BohrAtom(1)
print(hydrogen.summary())
```

This prints a complete summary: atomic number, mass, charge, shell configuration, and the ground-state energy level.

### Energy Levels

```python
for n in range(1, 6):
    E = hydrogen.energy_level(n)
    r = hydrogen.orbital_radius_pm(n)
    print(f"  n={n}: E = {E:8.3f} eV, r = {r:8.1f} pm")
```

Output:

```
  n=1: E =  -13.606 eV, r =     52.9 pm
  n=2: E =   -3.401 eV, r =    211.7 pm
  n=3: E =   -1.512 eV, r =    476.3 pm
  n=4: E =   -0.850 eV, r =    846.8 pm
  n=5: E =   -0.544 eV, r =   1323.1 pm
```

The energy scales as -13.6/n^2 eV, and the radius scales as n^2 * 52.9 pm.

### Ionization Energy

```python
IE = hydrogen.ionization_energy()
print(f"Ionization energy of H: {IE:.3f} eV")  # 13.606 eV
```

## 2. Spectral Transitions and the Balmer Series

When an electron drops from level n_i to n_f, it emits a photon. The Balmer series (n_f = 2) produces visible light:

```python
hydrogen = BohrAtom(1)

print("Balmer series (transitions to n=2):")
for n_i in range(3, 8):
    wavelength = hydrogen.transition_wavelength_nm(n_i, 2)
    energy = hydrogen.transition_energy(n_i, 2)
    print(f"  n={n_i} -> n=2: {wavelength:7.1f} nm  ({energy:.3f} eV)")
```

Output:

```
Balmer series (transitions to n=2):
  n=3 -> n=2:   656.5 nm  (1.889 eV)   # H-alpha, red
  n=4 -> n=2:   486.3 nm  (2.551 eV)   # H-beta, cyan
  n=5 -> n=2:   434.2 nm  (2.857 eV)   # H-gamma, violet
  n=6 -> n=2:   410.3 nm  (3.023 eV)   # H-delta, violet
  n=7 -> n=2:   397.1 nm  (3.124 eV)   # near-UV
```

The 656.5 nm H-alpha line is the reason hydrogen emission appears red -- this is the dominant visible transition.

## 3. Electron Configurations

For multi-electron atoms, MolBuilder uses the full quantum mechanical model with Aufbau filling, Hund's rule, and Pauli exclusion:

```python
from molbuilder.atomic.quantum_atom import QuantumAtom

carbon = QuantumAtom(6)
print(f"Carbon electron configuration: {carbon.electron_configuration_string()}")
print(f"Noble gas notation:            {carbon.noble_gas_notation()}")
print(f"Valence electrons:             {carbon.valence_electrons}")
print(f"Ionization energy:             {carbon.ionization_energy_eV():.2f} eV")
```

Output:

```
Carbon electron configuration: 1s2 2s2 2p2
Noble gas notation:            [He] 2s2 2p2
Valence electrons:             4
Ionization energy:             11.79 eV
```

### Exploring Subshell Details

Each subshell carries quantum numbers and occupation:

```python
for sub in carbon.subshells:
    zeff = carbon.effective_nuclear_charge(sub.n, sub.l)
    E = carbon.orbital_energy_eV(sub.n, sub.l)
    print(f"  {sub.n}{['s','p','d','f'][sub.l]}{sub.electron_count}: "
          f"Z_eff = {zeff:.2f}, E = {E:.2f} eV")
```

The effective nuclear charge (Slater's rules) explains why 2s electrons are lower in energy than 2p: they penetrate closer to the nucleus.

## 4. Aufbau Exceptions

Some elements deviate from the expected filling order. Iron is a classic example:

```python
iron = QuantumAtom(26)
print(f"Iron: {iron.electron_configuration_string()}")
print(f"Noble gas: {iron.noble_gas_notation()}")
print(f"Spin multiplicity: {iron.spin_multiplicity()}")

# Copper has the well-known [Ar] 3d10 4s1 exception
copper = QuantumAtom(29)
print(f"\nCopper: {copper.electron_configuration_string()}")
print(f"Noble gas: {copper.noble_gas_notation()}")
```

The spin multiplicity (2S + 1) tells you the number of unpaired electrons. High-spin Fe^2+ complexes have a multiplicity of 5 (four unpaired d-electrons), which is why iron compounds are often paramagnetic.

## 5. Orbital Visualization

MolBuilder renders atomic orbitals as 3D probability clouds using rejection sampling of |psi(r, theta, phi)|^2:

```python
from molbuilder.visualization.quantum_viz import (
    plot_orbital_3d,
    plot_radial_probability,
    plot_electron_configuration,
    visualize_atom,
)
from molbuilder.atomic.quantum_atom import QuantumAtom

# --- 3D orbital cloud for hydrogen 2p_z ---
plot_orbital_3d(n=2, l=1, m=0, Z=1, num_points=30000)

# --- Radial probability for 1s, 2s, 2p ---
plot_radial_probability([(1, 0), (2, 0), (2, 1)], Z=1)

# --- Electron configuration box diagram ---
carbon = QuantumAtom(6)
plot_electron_configuration(carbon)
```

### Multi-Panel Overview

The `visualize_atom()` function generates a four-panel figure showing radial probability, energy levels, configuration diagram, and a 3D orbital cloud -- all in one call:

```python
from molbuilder.visualization.quantum_viz import visualize_atom
from molbuilder.atomic.quantum_atom import QuantumAtom

carbon = QuantumAtom(6)
visualize_atom(carbon, orbital_nlm=(2, 1, 0))  # Show the 2p_z orbital
```

### Available Visualization Functions

| Function | What it shows |
|----------|---------------|
| `plot_orbital_3d(n, l, m)` | 3D probability cloud for a single orbital |
| `plot_radial_wavefunction(pairs)` | R_nl(r) radial wave function curves |
| `plot_radial_probability(pairs)` | r^2 |R_nl(r)|^2 radial probability |
| `plot_angular_distribution(l, m)` | Polar plot of angular probability |
| `plot_energy_levels(atom)` | Energy level diagram for all subshells |
| `plot_electron_configuration(atom)` | Box diagram with spin arrows |
| `visualize_atom(atom)` | Combined four-panel overview |

---

## Complete Example

```python
"""Quantum Mechanics: Atomic Structure and Orbitals -- complete working example.

Requirements: pip install molbuilder==1.1.1 matplotlib
"""
from molbuilder.atomic.bohr import BohrAtom
from molbuilder.atomic.quantum_atom import QuantumAtom
from molbuilder.visualization.quantum_viz import (
    plot_orbital_3d,
    plot_radial_probability,
    visualize_atom,
)

# --- 1. Hydrogen Bohr model ---
H = BohrAtom(1)
print("=== Hydrogen (Bohr model) ===")
for n in range(1, 5):
    print(f"  n={n}: E = {H.energy_level(n):.3f} eV, "
          f"r = {H.orbital_radius_pm(n):.1f} pm")

print("\nBalmer series:")
for ni in range(3, 7):
    wl = H.transition_wavelength_nm(ni, 2)
    print(f"  {ni} -> 2: {wl:.1f} nm")

# --- 2. Carbon quantum model ---
C = QuantumAtom(6)
print(f"\n=== Carbon ===")
print(f"Configuration: {C.electron_configuration_string()}")
print(f"Noble gas:     {C.noble_gas_notation()}")
print(f"Valence e-:    {C.valence_electrons}")
print(f"IE:            {C.ionization_energy_eV():.2f} eV")

# --- 3. Iron ---
Fe = QuantumAtom(26)
print(f"\n=== Iron ===")
print(f"Configuration: {Fe.electron_configuration_string()}")
print(f"Noble gas:     {Fe.noble_gas_notation()}")
print(f"Spin mult.:    {Fe.spin_multiplicity()}")

# --- 4. Visualizations (uncomment to display) ---
# plot_orbital_3d(n=2, l=1, m=0, Z=1)           # Hydrogen 2p_z
# plot_radial_probability([(1,0), (2,0), (2,1)]) # Compare s and p orbitals
# visualize_atom(C, orbital_nlm=(2, 1, 0))        # Carbon overview
print("\nUncomment visualization lines to display plots.")
```

---

## Next Steps

- [From SMILES to 3D](../From%20SMILES%20to%203D%20Building%20Molecules%20with%20MolBuilder/) -- move from atoms to molecules
- [3D Visualization and Conformational Analysis](../3D%20Visualization%20and%20Conformational%20Analysis/) -- visualize molecular geometry and torsional energy
- [Drug-Likeness Screening](../Drug-Likeness%20Screening%20with%20Lipinski%20Rule%20of%20Five/) -- apply molecular properties to drug design

---

## References

1. Bohr, N. On the Constitution of Atoms and Molecules. *Philos. Mag.* **1913**, *26* (151), 1--25. DOI: 10.1080/14786441308634955
2. Slater, J. C. Atomic Shielding Constants. *Phys. Rev.* **1930**, *36* (1), 57--64. DOI: 10.1103/PhysRev.36.57
3. Griffiths, D. J. *Introduction to Quantum Mechanics*, 3rd ed.; Cambridge University Press: Cambridge, 2018.
4. Atkins, P. W.; Friedman, R. S. *Molecular Quantum Mechanics*, 5th ed.; Oxford University Press: Oxford, 2011.

---

*Tested with MolBuilder 1.1.1, Python 3.13, matplotlib 3.9. This article was written with AI assistance using Claude (Anthropic) and reviewed by the author.*
