"""Molecule report generator.

Produces a comprehensive ASCII report about a molecule's composition,
bonding, functional groups, and connectivity.  All output is cp1252-safe.
"""

from __future__ import annotations

from collections import Counter

from typing import TYPE_CHECKING

from molbuilder.reports.text_formatter import (
    section_header,
    subsection_header,
    ascii_table,
    bullet_list,
    key_value_block,
    format_percent,
)

if TYPE_CHECKING:
    from molbuilder.reports import MoleculeLike

# Attempt to import functional-group detection.  If the reactions module
# is unavailable the report simply omits that section.
try:
    from molbuilder.reactions.functional_group_detect import detect_functional_groups as _detect_fg
    _HAS_FG_DETECT = True
except Exception:
    _HAS_FG_DETECT = False


# Standard atomic weights for the most common organic elements (g/mol).
_ATOMIC_WEIGHTS: dict[str, float] = {
    "H":   1.008,
    "He":  4.003,
    "Li":  6.941,
    "Be":  9.012,
    "B":  10.811,
    "C":  12.011,
    "N":  14.007,
    "O":  15.999,
    "F":  18.998,
    "Ne": 20.180,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.982,
    "Si": 28.086,
    "P":  30.974,
    "S":  32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K":  39.098,
    "Ca": 40.078,
    "Ti": 47.867,
    "Cr": 51.996,
    "Mn": 54.938,
    "Fe": 55.845,
    "Co": 58.933,
    "Ni": 58.693,
    "Cu": 63.546,
    "Zn": 65.380,
    "Br": 79.904,
    "I": 126.904,
    "Pd": 106.42,
    "Sn": 118.71,
    "Pt": 195.08,
}


def _hill_formula(counts: dict[str, int]) -> str:
    """Return molecular formula in Hill system order.

    Hill system: C first, then H, then everything else alphabetically.
    If no carbon, all elements alphabetically.
    """
    parts: list[str] = []
    remaining = dict(counts)

    if "C" in remaining:
        n = remaining.pop("C")
        parts.append("C" if n == 1 else f"C{n}")
        if "H" in remaining:
            n = remaining.pop("H")
            parts.append("H" if n == 1 else f"H{n}")

    for elem in sorted(remaining.keys()):
        n = remaining[elem]
        parts.append(elem if n == 1 else f"{elem}{n}")

    return "".join(parts)


def _molecular_weight(counts: dict[str, int]) -> float:
    """Compute molecular weight from element counts."""
    total = 0.0
    for elem, n in counts.items():
        total += _ATOMIC_WEIGHTS.get(elem, 0.0) * n
    return total


def _bond_order_label(order) -> str:
    """Human-readable bond order label."""
    try:
        order_val = float(order)
    except (TypeError, ValueError):
        return "unknown"
    if order_val == 1.0:
        return "single"
    if order_val == 1.5:
        return "aromatic"
    if order_val == 2.0:
        return "double"
    if order_val == 3.0:
        return "triple"
    return f"order {order_val}"


# =====================================================================
#  Public API
# =====================================================================

def generate_molecule_report(mol: MoleculeLike) -> str:
    """Generate a comprehensive ASCII report about a molecule.

    Uses duck typing -- *mol* should expose:

    * ``.name`` -- molecule name (str)
    * ``.atoms`` -- list of atom objects with ``.symbol``
    * ``.bonds`` -- list of bond objects with ``.atom_i``, ``.atom_j``, ``.order``
    * ``.neighbors(idx)`` -- list of bonded-atom indices
    * ``.get_bond(i, j)`` -- Bond object or None

    Returns a single multi-line string.
    """
    if not hasattr(mol, 'atoms'):
        raise TypeError(
            f"mol must have an 'atoms' attribute, "
            f"got {type(mol).__name__}"
        )

    lines: list[str] = []

    # ------------------------------------------------------------------
    # 1. Header
    # ------------------------------------------------------------------
    mol_name = getattr(mol, "name", "Unknown") or "Unknown"
    lines.append(section_header(f"Molecule Report: {mol_name}"))
    lines.append("")

    atoms = getattr(mol, "atoms", [])
    bonds = getattr(mol, "bonds", [])

    # ------------------------------------------------------------------
    # 2. Basic Properties
    # ------------------------------------------------------------------
    elem_counts: Counter[str] = Counter()
    for atom in atoms:
        elem_counts[atom.symbol] += 1

    formula = _hill_formula(dict(elem_counts))
    mw = _molecular_weight(dict(elem_counts))

    lines.append(subsection_header("Basic Properties"))
    props = [
        ("Molecular Formula", formula),
        ("Molecular Weight",  f"{mw:.3f} g/mol"),
        ("Total Atoms",       str(len(atoms))),
        ("Total Bonds",       str(len(bonds))),
    ]
    lines.append(key_value_block(props))
    lines.append("")

    # ------------------------------------------------------------------
    # 3. Atom Composition
    # ------------------------------------------------------------------
    lines.append(subsection_header("Atom Composition"))
    comp_headers = ["Element", "Count", "Mass (g/mol)", "Mass %"]
    comp_rows: list[list[str]] = []
    for elem in sorted(elem_counts.keys()):
        count = elem_counts[elem]
        mass = _ATOMIC_WEIGHTS.get(elem, 0.0) * count
        pct = (mass / mw * 100.0) if mw > 0 else 0.0
        comp_rows.append([
            elem,
            str(count),
            f"{mass:.3f}",
            format_percent(pct),
        ])
    lines.append(ascii_table(
        comp_headers, comp_rows,
        alignments=["l", "r", "r", "r"],
        min_widths=[8, 6, 12, 8],
    ))
    lines.append("")

    # ------------------------------------------------------------------
    # 4. Bond Summary
    # ------------------------------------------------------------------
    lines.append(subsection_header("Bond Summary"))
    order_counts: Counter[str] = Counter()
    for bond in bonds:
        label = _bond_order_label(bond.order)
        order_counts[label] += 1

    bond_headers = ["Bond Type", "Count"]
    bond_rows = [[btype, str(cnt)]
                 for btype, cnt in sorted(order_counts.items())]
    lines.append(ascii_table(
        bond_headers, bond_rows,
        alignments=["l", "r"],
        min_widths=[12, 6],
    ))
    lines.append("")

    # ------------------------------------------------------------------
    # 5. Functional Groups
    # ------------------------------------------------------------------
    if _HAS_FG_DETECT:
        lines.append(subsection_header("Functional Groups"))
        try:
            groups = _detect_fg(mol)
            if groups:
                fg_names = [g.name for g in groups]
                fg_counts: Counter[str] = Counter(fg_names)
                fg_items = [
                    f"{name} (x{cnt})" if cnt > 1 else name
                    for name, cnt in sorted(fg_counts.items())
                ]
                lines.append(bullet_list(fg_items))
            else:
                lines.append("  No common functional groups detected.")
        except Exception:
            lines.append("  Functional group detection unavailable.")
        lines.append("")

    # ------------------------------------------------------------------
    # 6. Connectivity
    # ------------------------------------------------------------------
    lines.append(subsection_header("Connectivity"))

    # Degree distribution
    degree_counts: Counter[int] = Counter()
    for idx in range(len(atoms)):
        try:
            nbrs = mol.neighbors(idx)
            degree_counts[len(nbrs)] += 1
        except Exception:
            pass

    if degree_counts:
        deg_headers = ["Degree", "Atom Count"]
        deg_rows = [[str(d), str(c)]
                    for d, c in sorted(degree_counts.items())]
        lines.append(ascii_table(
            deg_headers, deg_rows,
            alignments=["r", "r"],
            min_widths=[8, 12],
        ))
    lines.append("")

    # Ring detection heuristic (Euler formula: rings = bonds - atoms + 1
    # for connected graph).  This gives the cyclomatic number.
    n_atoms = len(atoms)
    n_bonds = len(bonds)
    ring_count = n_bonds - n_atoms + 1
    if ring_count > 0:
        lines.append(f"  Ring structures detected (cyclomatic number: {ring_count})")
    else:
        lines.append("  No ring structures detected (acyclic molecule)")
    lines.append("")

    # Footer
    lines.append("=" * 70)
    lines.append("  End of Molecule Report")
    lines.append("=" * 70)

    return "\n".join(lines)
