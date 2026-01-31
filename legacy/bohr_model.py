"""
Bohr Model of the Atom

Represents atoms using the Bohr model with quantized energy levels,
electron shell configurations, and orbital mechanics. Provides both
computational physics and animated visualization.

Physics reference:
    - Orbital radius:  r_n = n^2 * a_0 / Z
    - Energy level:    E_n = -13.6 eV * Z^2 / n^2
    - Photon energy:   delta_E = E_final - E_initial
    - Wavelength:      lambda = h*c / |delta_E|
    - Orbital velocity: v_n = (Z * e^2) / (n * hbar)
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
BOHR_RADIUS_M = 5.29177210903e-11       # metres
BOHR_RADIUS_PM = 52.9177210903          # picometres
PLANCK_CONSTANT = 6.62607015e-34         # J s
HBAR = PLANCK_CONSTANT / (2 * math.pi)  # reduced Planck constant
SPEED_OF_LIGHT = 2.99792458e8           # m/s
ELECTRON_CHARGE = 1.602176634e-19       # C
ELECTRON_MASS = 9.1093837015e-31        # kg
COULOMB_CONSTANT = 8.9875517873681764e9 # N m^2 / C^2
EV_TO_JOULES = 1.602176634e-19
RYDBERG_ENERGY_EV = 13.605693122994     # eV

# ---------------------------------------------------------------------------
# Element data  (atomic_number -> symbol, name, standard_atomic_weight)
# ---------------------------------------------------------------------------
ELEMENTS = {
    1:  ("H",  "Hydrogen",      1.008),
    2:  ("He", "Helium",        4.003),
    3:  ("Li", "Lithium",       6.941),
    4:  ("Be", "Beryllium",     9.012),
    5:  ("B",  "Boron",        10.811),
    6:  ("C",  "Carbon",       12.011),
    7:  ("N",  "Nitrogen",     14.007),
    8:  ("O",  "Oxygen",       15.999),
    9:  ("F",  "Fluorine",     18.998),
    10: ("Ne", "Neon",         20.180),
    11: ("Na", "Sodium",       22.990),
    12: ("Mg", "Magnesium",    24.305),
    13: ("Al", "Aluminium",    26.982),
    14: ("Si", "Silicon",      28.086),
    15: ("P",  "Phosphorus",   30.974),
    16: ("S",  "Sulfur",       32.065),
    17: ("Cl", "Chlorine",     35.453),
    18: ("Ar", "Argon",        39.948),
    19: ("K",  "Potassium",    39.098),
    20: ("Ca", "Calcium",      40.078),
    21: ("Sc", "Scandium",     44.956),
    22: ("Ti", "Titanium",     47.867),
    23: ("V",  "Vanadium",     50.942),
    24: ("Cr", "Chromium",     51.996),
    25: ("Mn", "Manganese",    54.938),
    26: ("Fe", "Iron",         55.845),
    27: ("Co", "Cobalt",       58.933),
    28: ("Ni", "Nickel",       58.693),
    29: ("Cu", "Copper",       63.546),
    30: ("Zn", "Zinc",         65.380),
    31: ("Ga", "Gallium",      69.723),
    32: ("Ge", "Germanium",    72.630),
    33: ("As", "Arsenic",      74.922),
    34: ("Se", "Selenium",     78.971),
    35: ("Br", "Bromine",      79.904),
    36: ("Kr", "Krypton",      83.798),
    37: ("Rb", "Rubidium",     85.468),
    38: ("Sr", "Strontium",    87.620),
    39: ("Y",  "Yttrium",      88.906),
    40: ("Zr", "Zirconium",    91.224),
    41: ("Nb", "Niobium",      92.906),
    42: ("Mo", "Molybdenum",   95.950),
    43: ("Tc", "Technetium",   98.000),
    44: ("Ru", "Ruthenium",   101.070),
    45: ("Rh", "Rhodium",     102.906),
    46: ("Pd", "Palladium",   106.420),
    47: ("Ag", "Silver",      107.868),
    48: ("Cd", "Cadmium",     112.414),
    49: ("In", "Indium",      114.818),
    50: ("Sn", "Tin",         118.710),
    51: ("Sb", "Antimony",    121.760),
    52: ("Te", "Tellurium",   127.600),
    53: ("I",  "Iodine",      126.904),
    54: ("Xe", "Xenon",       131.293),
    55: ("Cs", "Caesium",     132.905),
    56: ("Ba", "Barium",      137.327),
    57: ("La", "Lanthanum",   138.905),
    58: ("Ce", "Cerium",      140.116),
    59: ("Pr", "Praseodymium",140.908),
    60: ("Nd", "Neodymium",   144.242),
    61: ("Pm", "Promethium",  145.000),
    62: ("Sm", "Samarium",    150.360),
    63: ("Eu", "Europium",    151.964),
    64: ("Gd", "Gadolinium",  157.250),
    65: ("Tb", "Terbium",     158.925),
    66: ("Dy", "Dysprosium",  162.500),
    67: ("Ho", "Holmium",     164.930),
    68: ("Er", "Erbium",      167.259),
    69: ("Tm", "Thulium",     168.934),
    70: ("Yb", "Ytterbium",   173.045),
    71: ("Lu", "Lutetium",    174.967),
    72: ("Hf", "Hafnium",     178.490),
    73: ("Ta", "Tantalum",    180.948),
    74: ("W",  "Tungsten",    183.840),
    75: ("Re", "Rhenium",     186.207),
    76: ("Os", "Osmium",      190.230),
    77: ("Ir", "Iridium",     192.217),
    78: ("Pt", "Platinum",    195.084),
    79: ("Au", "Gold",        196.967),
    80: ("Hg", "Mercury",     200.592),
    81: ("Tl", "Thallium",    204.383),
    82: ("Pb", "Lead",        207.200),
    83: ("Bi", "Bismuth",     208.980),
    84: ("Po", "Polonium",    209.000),
    85: ("At", "Astatine",    210.000),
    86: ("Rn", "Radon",       222.000),
    87: ("Fr", "Francium",    223.000),
    88: ("Ra", "Radium",      226.000),
    89: ("Ac", "Actinium",    227.000),
    90: ("Th", "Thorium",     232.038),
    91: ("Pa", "Protactinium",231.036),
    92: ("U",  "Uranium",     238.029),
    93: ("Np", "Neptunium",   237.000),
    94: ("Pu", "Plutonium",   244.000),
    95: ("Am", "Americium",   243.000),
    96: ("Cm", "Curium",      247.000),
    97: ("Bk", "Berkelium",   247.000),
    98: ("Cf", "Californium", 251.000),
    99: ("Es", "Einsteinium", 252.000),
    100:("Fm", "Fermium",     257.000),
    101:("Md", "Mendelevium", 258.000),
    102:("No", "Nobelium",    259.000),
    103:("Lr", "Lawrencium",  266.000),
    104:("Rf", "Rutherfordium",267.000),
    105:("Db", "Dubnium",     268.000),
    106:("Sg", "Seaborgium",  269.000),
    107:("Bh", "Bohrium",     270.000),
    108:("Hs", "Hassium",     277.000),
    109:("Mt", "Meitnerium",  278.000),
    110:("Ds", "Darmstadtium",281.000),
    111:("Rg", "Roentgenium", 282.000),
    112:("Cn", "Copernicium", 285.000),
    113:("Nh", "Nihonium",    286.000),
    114:("Fl", "Flerovium",   289.000),
    115:("Mc", "Moscovium",   290.000),
    116:("Lv", "Livermorium", 293.000),
    117:("Ts", "Tennessine",  294.000),
    118:("Og", "Oganesson",   294.000),
}

SYMBOL_TO_Z = {v[0]: k for k, v in ELEMENTS.items()}

# Maximum electrons per shell: 2n^2
MAX_ELECTRONS_PER_SHELL = [2 * n**2 for n in range(1, 8)]  # shells 1-7


# ---------------------------------------------------------------------------
# Core Bohr model class
# ---------------------------------------------------------------------------
class BohrAtom:
    """Represents an atom using the Bohr model.

    Parameters
    ----------
    atomic_number : int
        Number of protons (Z). Determines the element.
    mass_number : int or None
        Protons + neutrons. Defaults to rounded standard atomic weight.
    charge : int
        Net ionic charge (0 = neutral, +1 = lost one electron, etc.).
    """

    def __init__(self, atomic_number: int, mass_number: int = None, charge: int = 0):
        if atomic_number < 1 or atomic_number > 118:
            raise ValueError(f"Atomic number must be 1-118, got {atomic_number}")

        self.atomic_number = atomic_number  # Z
        self.protons = atomic_number
        symbol, name, weight = ELEMENTS[atomic_number]
        self.symbol = symbol
        self.name = name
        self.standard_atomic_weight = weight

        if mass_number is None:
            self.mass_number = round(weight)
        else:
            self.mass_number = mass_number

        self.neutrons = self.mass_number - self.protons
        self.charge = charge
        self.num_electrons = self.protons - charge

        if self.num_electrons < 0:
            raise ValueError("Charge exceeds number of protons (no electrons left)")

        self.shell_config = self._compute_shell_config()
        self.num_shells = len(self.shell_config)

    # ----- electron configuration -----

    def _compute_shell_config(self) -> list[int]:
        """Fill electron shells following the 2n^2 rule."""
        remaining = self.num_electrons
        shells = []
        n = 1
        while remaining > 0 and n <= 7:
            capacity = 2 * n**2
            electrons_in_shell = min(remaining, capacity)
            shells.append(electrons_in_shell)
            remaining -= electrons_in_shell
            n += 1
        return shells

    # ----- Bohr radius & energy -----

    def orbital_radius(self, n: int) -> float:
        """Radius of the n-th Bohr orbit in metres.

        r_n = n^2 * a_0 / Z
        """
        if n < 1:
            raise ValueError("Principal quantum number n must be >= 1")
        return (n**2 * BOHR_RADIUS_M) / self.atomic_number

    def orbital_radius_pm(self, n: int) -> float:
        """Radius of the n-th Bohr orbit in picometres."""
        return self.orbital_radius(n) * 1e12

    def energy_level(self, n: int) -> float:
        """Energy of the n-th level in eV.

        E_n = -13.6 eV * Z^2 / n^2
        """
        if n < 1:
            raise ValueError("Principal quantum number n must be >= 1")
        return -RYDBERG_ENERGY_EV * (self.atomic_number**2) / (n**2)

    def energy_level_joules(self, n: int) -> float:
        """Energy of the n-th level in joules."""
        return self.energy_level(n) * EV_TO_JOULES

    def orbital_velocity(self, n: int) -> float:
        """Velocity of the electron in the n-th orbit in m/s.

        v_n = Z * e^2 / (n * 4*pi*eps_0 * hbar)
            = Z * alpha * c / n
        where alpha ~ 1/137 is the fine-structure constant.
        """
        alpha = ELECTRON_CHARGE**2 / (4 * math.pi * 8.8541878128e-12 * HBAR * SPEED_OF_LIGHT)
        return self.atomic_number * alpha * SPEED_OF_LIGHT / n

    def orbital_period(self, n: int) -> float:
        """Orbital period of the electron in shell n, in seconds."""
        r = self.orbital_radius(n)
        v = self.orbital_velocity(n)
        return 2 * math.pi * r / v

    # ----- transitions -----

    def transition_energy(self, n_initial: int, n_final: int) -> float:
        """Energy of a photon emitted (positive) or absorbed (negative)
        during a transition, in eV."""
        return self.energy_level(n_initial) - self.energy_level(n_final)

    def transition_wavelength(self, n_initial: int, n_final: int) -> float:
        """Wavelength in metres of a photon from a transition between levels."""
        delta_e = abs(self.transition_energy(n_initial, n_final)) * EV_TO_JOULES
        if delta_e == 0:
            return float('inf')
        return PLANCK_CONSTANT * SPEED_OF_LIGHT / delta_e

    def transition_wavelength_nm(self, n_initial: int, n_final: int) -> float:
        """Wavelength in nanometres."""
        return self.transition_wavelength(n_initial, n_final) * 1e9

    def ionization_energy(self) -> float:
        """Minimum energy (eV) to remove the outermost electron (Bohr approx)."""
        if self.num_shells == 0:
            return 0.0
        return -self.energy_level(self.num_shells)

    # ----- electrostatic force on electron -----

    def coulomb_force_on_electron(self, n: int) -> float:
        """Attractive Coulomb force on an electron in shell n, in newtons."""
        r = self.orbital_radius(n)
        return COULOMB_CONSTANT * (self.protons * ELECTRON_CHARGE**2) / r**2

    # ----- representation -----

    def __repr__(self):
        iso = f"-{self.mass_number}" if self.mass_number != round(self.standard_atomic_weight) else ""
        charge_str = ""
        if self.charge > 0:
            charge_str = f" (+{self.charge})"
        elif self.charge < 0:
            charge_str = f" ({self.charge})"
        return (f"BohrAtom({self.symbol}{iso}, Z={self.atomic_number}, "
                f"e={self.num_electrons}{charge_str}, "
                f"shells={self.shell_config})")

    def summary(self) -> str:
        """Return a human-readable summary of the atom."""
        lines = [
            f"{'='*50}",
            f"  {self.name} ({self.symbol})",
            f"  Atomic number (Z): {self.atomic_number}",
            f"  Mass number (A):   {self.mass_number}",
            f"  Protons:  {self.protons}   Neutrons: {self.neutrons}",
            f"  Electrons: {self.num_electrons}  Charge: {self.charge:+d}",
            f"  Shell configuration: {self.shell_config}",
            f"{'='*50}",
            f"  Shell   Electrons   Radius (pm)     Energy (eV)",
            f"  {'-'*46}",
        ]
        for i, count in enumerate(self.shell_config, start=1):
            r = self.orbital_radius_pm(i)
            e = self.energy_level(i)
            lines.append(f"  n={i:<4}  {count:<10}  {r:>10.2f}      {e:>10.4f}")
        lines.append(f"  {'-'*46}")
        lines.append(f"  Ionization energy (Bohr): {self.ionization_energy():.4f} eV")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------
def wavelength_to_rgb(wavelength_nm: float) -> tuple:
    """Convert a visible-light wavelength (380-780 nm) to an RGB tuple.
    Returns white for wavelengths outside visible range."""
    if wavelength_nm < 380 or wavelength_nm > 780:
        return (1.0, 1.0, 1.0)

    if wavelength_nm < 440:
        r = -(wavelength_nm - 440) / (440 - 380)
        g = 0.0
        b = 1.0
    elif wavelength_nm < 490:
        r = 0.0
        g = (wavelength_nm - 440) / (490 - 440)
        b = 1.0
    elif wavelength_nm < 510:
        r = 0.0
        g = 1.0
        b = -(wavelength_nm - 510) / (510 - 490)
    elif wavelength_nm < 580:
        r = (wavelength_nm - 510) / (580 - 510)
        g = 1.0
        b = 0.0
    elif wavelength_nm < 645:
        r = 1.0
        g = -(wavelength_nm - 645) / (645 - 580)
        b = 0.0
    else:
        r = 1.0
        g = 0.0
        b = 0.0

    # intensity fall-off at edges of visible spectrum
    if wavelength_nm < 420:
        factor = 0.3 + 0.7 * (wavelength_nm - 380) / (420 - 380)
    elif wavelength_nm > 700:
        factor = 0.3 + 0.7 * (780 - wavelength_nm) / (780 - 700)
    else:
        factor = 1.0

    return (r * factor, g * factor, b * factor)


def visualize(atom: BohrAtom, animate: bool = True, interval_ms: int = 30):
    """Render an animated Bohr model diagram of the atom.

    Parameters
    ----------
    atom : BohrAtom
        The atom to visualize.
    animate : bool
        If True, electrons orbit the nucleus. If False, show a static frame.
    interval_ms : int
        Milliseconds between animation frames.
    """
    n_shells = atom.num_shells
    if n_shells == 0:
        print("No electrons to visualize.")
        return

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), facecolor="black")
    ax.set_facecolor("black")
    ax.set_aspect("equal")
    ax.set_xlim(-n_shells - 1, n_shells + 1)
    ax.set_ylim(-n_shells - 1, n_shells + 1)
    ax.axis("off")

    # Title
    charge_label = ""
    if atom.charge > 0:
        charge_label = f"$^{{+{atom.charge}}}$"
    elif atom.charge < 0:
        charge_label = f"$^{{{atom.charge}}}$"
    ax.set_title(
        f"Bohr Model - {atom.name} ({atom.symbol}{charge_label})  "
        f"Z={atom.atomic_number}  e={atom.num_electrons}",
        color="white", fontsize=14, pad=12,
    )

    # ----- Nucleus -----
    nucleus_radius = 0.25 + 0.03 * atom.atomic_number**0.33
    nucleus = Circle((0, 0), nucleus_radius, color="#ff6633", zorder=10)
    ax.add_patch(nucleus)
    ax.text(0, 0, f"{atom.protons}p\n{atom.neutrons}n",
            ha="center", va="center", fontsize=7, color="white",
            fontweight="bold", zorder=11)

    # ----- Orbital rings -----
    shell_radii = []
    for n in range(1, n_shells + 1):
        r = n  # use integer spacing for visual clarity
        shell_radii.append(r)
        orbit = Circle((0, 0), r, fill=False, edgecolor="#334466",
                        linewidth=0.8, linestyle="--", zorder=1)
        ax.add_patch(orbit)
        ax.text(r + 0.15, 0.15, f"n={n}", fontsize=7, color="#5588aa", zorder=2)

    # ----- Electron dots (initial positions) -----
    electron_artists = []
    electron_positions = []  # (shell_index, angle_offset, shell_radius)

    for shell_idx, count in enumerate(atom.shell_config):
        r = shell_radii[shell_idx]
        for e in range(count):
            angle = 2 * math.pi * e / count
            x = r * math.cos(angle)
            y = r * math.sin(angle)
            dot = ax.plot(x, y, 'o', color="#44ccff", markersize=5, zorder=5)[0]
            electron_artists.append(dot)
            electron_positions.append((shell_idx, angle, r, count))

    # ----- Shell electron count labels -----
    for shell_idx, count in enumerate(atom.shell_config):
        r = shell_radii[shell_idx]
        ax.text(-r - 0.15, -0.25, str(count), fontsize=8, color="#88bbdd",
                ha="right", zorder=2)

    # ----- Energy level sidebar -----
    sidebar_x = n_shells + 0.6
    e_min = atom.energy_level(1)
    e_max = atom.energy_level(n_shells) if n_shells > 1 else e_min * 0.1
    e_range = abs(e_max - e_min) if abs(e_max - e_min) > 0 else 1.0
    bar_bottom = -n_shells
    bar_height = 2 * n_shells

    for n in range(1, n_shells + 1):
        e = atom.energy_level(n)
        if n_shells > 1:
            y_pos = bar_bottom + bar_height * (e - e_min) / e_range
        else:
            y_pos = 0
        ax.plot([sidebar_x, sidebar_x + 0.5], [y_pos, y_pos],
                color="#ffaa33", linewidth=1.5, zorder=3)
        ax.text(sidebar_x + 0.6, y_pos, f"{e:.2f} eV",
                fontsize=6, color="#ffcc66", va="center", zorder=3)

    # ----- Animation -----
    def update(frame):
        for i, (shell_idx, angle0, r, count) in enumerate(electron_positions):
            n = shell_idx + 1
            angular_speed = 0.06 / n  # outer shells orbit slower
            angle = angle0 + angular_speed * frame
            x = r * math.cos(angle)
            y = r * math.sin(angle)
            electron_artists[i].set_data([x], [y])
        return electron_artists

    if animate:
        anim = animation.FuncAnimation(
            fig, update, frames=None, interval=interval_ms, blit=True,
        )

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Convenience constructors
# ---------------------------------------------------------------------------
def from_symbol(symbol: str, mass_number: int = None, charge: int = 0) -> BohrAtom:
    """Create a BohrAtom from an element symbol string."""
    symbol = symbol.capitalize()
    if len(symbol) > 1:
        symbol = symbol[0].upper() + symbol[1:].lower()
    if symbol not in SYMBOL_TO_Z:
        raise ValueError(f"Unknown element symbol: {symbol}")
    return BohrAtom(SYMBOL_TO_Z[symbol], mass_number, charge)


def from_name(name: str, mass_number: int = None, charge: int = 0) -> BohrAtom:
    """Create a BohrAtom from an element name string."""
    name_lower = name.strip().lower()
    for z, (sym, elem_name, _) in ELEMENTS.items():
        if elem_name.lower() == name_lower:
            return BohrAtom(z, mass_number, charge)
    raise ValueError(f"Unknown element name: {name}")


# ---------------------------------------------------------------------------
# Main -demo
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Bohr Model of the Atom -Interactive Demo\n")

    # Build a few example atoms
    hydrogen = BohrAtom(1)
    carbon = BohrAtom(6)
    sodium = BohrAtom(11)
    iron = from_symbol("Fe")

    for atom in [hydrogen, carbon, sodium, iron]:
        print(atom.summary())
        print()

    # Hydrogen spectral series (Balmer: transitions to n=2)
    print("Hydrogen Balmer series (visible spectrum):")
    h = BohrAtom(1)
    for ni in range(3, 7):
        wl = h.transition_wavelength_nm(ni, 2)
        rgb = wavelength_to_rgb(wl)
        print(f"  n={ni} -> n=2 : {wl:.1f} nm  (RGB ~ {tuple(round(c,2) for c in rgb)})")
    print()

    # Visualize -pick an element
    print("Launching animated visualisation for Carbon (Z=6)...")
    visualize(carbon)
