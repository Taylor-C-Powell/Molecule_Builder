"""Functional group detection via subgraph pattern matching.

This module walks the atoms and bonds of a ``Molecule`` object and
identifies common organic functional groups by inspecting the local
neighbourhood of each atom.  The detection is heuristic (based on
connectivity and bond order) rather than relying on SMARTS matching
so that it works with the graph representation already available in
the ``molbuilder.molecule.graph`` module.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from molbuilder.molecule.graph import Molecule, Hybridization


# =====================================================================
#  Data class
# =====================================================================

@dataclass
class FunctionalGroup:
    """A detected functional group occurrence in a molecule.

    Attributes
    ----------
    name : str
        Human-readable name (e.g. ``"alcohol"``, ``"ketone"``).
    smarts_like : str
        Simplified pattern description for display.
    atoms : list[int]
        All atom indices participating in the group.
    center : int
        Primary (most characteristic) atom index.
    """
    name: str
    smarts_like: str
    atoms: list[int] = field(default_factory=list)
    center: int = -1

    def __repr__(self) -> str:
        return f"FunctionalGroup({self.name!r}, center={self.center})"


# =====================================================================
#  Helper utilities
# =====================================================================


# Standard valences for implicit hydrogen inference (same data as
# smiles.tokenizer.DEFAULT_VALENCE but upper-cased and used here so
# that FG detection does not depend on the SMILES subpackage).
_STANDARD_VALENCE: dict[str, list[int]] = {
    "B": [3], "C": [4], "N": [3, 5], "O": [2], "P": [3, 5],
    "S": [2, 4, 6], "F": [1], "Cl": [1], "Br": [1], "I": [1],
}


def _element(mol: Molecule, idx: int) -> str:
    """Return the element symbol of atom *idx* (upper-cased first letter)."""
    return mol.atoms[idx].symbol


def _neighbors(mol: Molecule, idx: int) -> list[int]:
    """Return indices of atoms bonded to *idx*."""
    return mol.neighbors(idx)


def _neighbor_elements(mol: Molecule, idx: int) -> list[str]:
    """Return element symbols of all neighbours of atom *idx*."""
    return [_element(mol, n) for n in _neighbors(mol, idx)]


def _bond_order(mol: Molecule, i: int, j: int) -> float:
    """Return the bond order between atoms *i* and *j*, or 0.0 if no bond."""
    bond = mol.get_bond(i, j)
    return bond.order if bond is not None else 0.0


def _sum_bond_orders(mol: Molecule, idx: int) -> int:
    """Sum of bond orders for all bonds on atom *idx*."""
    total = 0
    for n in _neighbors(mol, idx):
        total += int(_bond_order(mol, idx, n))
    return total


def _h_count(mol: Molecule, idx: int) -> int:
    """Total hydrogen count on atom *idx* (explicit + implicit).

    Explicit H atoms are counted from the neighbour list.  Implicit H
    atoms are inferred from standard valence rules when the molecule
    has fewer explicit neighbours than expected.  This makes FG
    detection work regardless of whether H atoms are represented as
    explicit nodes in the graph (SMILES-built molecules) or are
    absent (PDB/XYZ imports without H).
    """
    explicit_h = sum(1 for e in _neighbor_elements(mol, idx) if e == "H")
    if explicit_h > 0:
        return explicit_h

    # No explicit H found -- infer from valence rules
    sym = _element(mol, idx)
    valences = _STANDARD_VALENCE.get(sym)
    if valences is None:
        return 0
    bond_order_sum = _sum_bond_orders(mol, idx)
    # Pick the smallest standard valence that accommodates current bonds
    for v in valences:
        implicit = v - bond_order_sum
        if implicit >= 0:
            return implicit
    return 0


def _has_h(mol: Molecule, idx: int) -> bool:
    """Return True if atom *idx* has at least one hydrogen (explicit or implicit)."""
    return _h_count(mol, idx) > 0


def _count_element_neighbors(mol: Molecule, idx: int, elem: str) -> int:
    """Count how many neighbours of *idx* have element *elem*."""
    return sum(1 for e in _neighbor_elements(mol, idx) if e == elem)


def _double_bonded_to(mol: Molecule, idx: int, elem: str) -> list[int]:
    """Return neighbour indices that are *elem* and double-bonded to *idx*."""
    result = []
    for n in _neighbors(mol, idx):
        if _element(mol, n) == elem and _bond_order(mol, idx, n) == 2.0:
            result.append(n)
    return result


def _single_bonded_to(mol: Molecule, idx: int, elem: str) -> list[int]:
    """Return neighbour indices that are *elem* and single-bonded to *idx*."""
    result = []
    for n in _neighbors(mol, idx):
        if _element(mol, n) == elem and _bond_order(mol, idx, n) == 1.0:
            result.append(n)
    return result


def _triple_bonded_to(mol: Molecule, idx: int, elem: str) -> list[int]:
    """Return neighbour indices that are *elem* and triple-bonded to *idx*."""
    result = []
    for n in _neighbors(mol, idx):
        if _element(mol, n) == elem and _bond_order(mol, idx, n) == 3.0:
            result.append(n)
    return result


# =====================================================================
#  Master dispatcher
# =====================================================================

def detect_functional_groups(mol: Molecule) -> list[FunctionalGroup]:
    """Detect all recognisable functional groups in *mol*.

    Returns a list of ``FunctionalGroup`` instances, one per occurrence.
    The same atom may appear in more than one group (e.g. an ester
    contains both a C=O and a C-O-C linkage).
    """
    groups: list[FunctionalGroup] = []
    groups.extend(_detect_carboxylic_acids(mol))
    groups.extend(_detect_esters(mol))
    groups.extend(_detect_amides(mol))
    groups.extend(_detect_aldehydes(mol))
    groups.extend(_detect_ketones(mol))
    groups.extend(_detect_alcohols(mol))
    groups.extend(_detect_amines(mol))
    groups.extend(_detect_alkyl_halides(mol))
    groups.extend(_detect_alkenes(mol))
    groups.extend(_detect_alkynes(mol))
    groups.extend(_detect_ethers(mol))
    groups.extend(_detect_thiols(mol))
    groups.extend(_detect_nitriles(mol))
    groups.extend(_detect_nitro(mol))
    groups.extend(_detect_aromatic_rings(mol))
    groups.extend(_detect_epoxides(mol))
    groups.extend(_detect_acid_chlorides(mol))
    groups.extend(_detect_anhydrides(mol))
    groups.extend(_detect_sulfoxides(mol))
    groups.extend(_detect_sulfones(mol))
    groups.extend(_detect_imines(mol))
    groups.extend(_detect_boronic_acids(mol))
    groups.extend(_detect_phosphonates(mol))
    groups.extend(_detect_sulfonamides(mol))
    return groups


# =====================================================================
#  Individual detectors
# =====================================================================

def _detect_alcohols(mol: Molecule) -> list[FunctionalGroup]:
    """Alcohol: O bonded to C with an H (explicit or implicit).

    The O must be single-bonded to C and not part of a C=O or ester
    linkage.  Works with both explicit H in the graph and implicit H
    inferred from valence rules.
    """
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "O":
            continue
        nbrs = _neighbors(mol, idx)
        elems = [_element(mol, n) for n in nbrs]

        # Need at least one C neighbour, single-bonded
        c_indices = [nbrs[i] for i, e in enumerate(elems) if e == "C"]
        if not c_indices:
            continue

        for c_idx in c_indices:
            if _bond_order(mol, idx, c_idx) != 1.0:
                continue
            # Check for H: explicit neighbour OR implicit from valence
            if "H" in elems:
                h_idx = nbrs[elems.index("H")]
                found.append(FunctionalGroup(
                    name="alcohol", smarts_like="[C]-[OH]",
                    atoms=[c_idx, idx, h_idx], center=idx,
                ))
                break
            elif _h_count(mol, idx) >= 1:
                # Implicit H -- no explicit H atom index to record
                found.append(FunctionalGroup(
                    name="alcohol", smarts_like="[C]-[OH]",
                    atoms=[c_idx, idx], center=idx,
                ))
                break
    return found


def _detect_aldehydes(mol: Molecule) -> list[FunctionalGroup]:
    """Aldehyde: C=O where C also has an H (terminal carbonyl)."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        if not dbl_o:
            continue
        if _has_h(mol, idx):
            o_idx = dbl_o[0]
            found.append(FunctionalGroup(
                name="aldehyde", smarts_like="[CX3H1](=O)",
                atoms=[idx, o_idx], center=idx,
            ))
    return found


def _detect_ketones(mol: Molecule) -> list[FunctionalGroup]:
    """Ketone: C=O where C is bonded to two other carbons (no H, no O-single)."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        if not dbl_o:
            continue
        # Must not have H on carbonyl C (that would be aldehyde)
        if _has_h(mol, idx):
            continue
        # The other two neighbours should both be C
        c_nbrs = _single_bonded_to(mol, idx, "C")
        if len(c_nbrs) >= 2:
            o_idx = dbl_o[0]
            found.append(FunctionalGroup(
                name="ketone", smarts_like="[CX3](=O)([C])[C]",
                atoms=[idx, o_idx] + c_nbrs[:2], center=idx,
            ))
    return found


def _detect_carboxylic_acids(mol: Molecule) -> list[FunctionalGroup]:
    """Carboxylic acid: C with C=O and C-OH."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        sgl_o = _single_bonded_to(mol, idx, "O")
        if dbl_o and sgl_o:
            # Check that the single-bonded O has an H
            for o_idx in sgl_o:
                if _has_h(mol, o_idx):
                    found.append(FunctionalGroup(
                        name="carboxylic_acid",
                        smarts_like="[CX3](=O)[OH]",
                        atoms=[idx, dbl_o[0], o_idx], center=idx,
                    ))
                    break
    return found


def _detect_esters(mol: Molecule) -> list[FunctionalGroup]:
    """Ester: C(=O)-O-C where the single-bonded O has no H."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        sgl_o = _single_bonded_to(mol, idx, "O")
        if dbl_o and sgl_o:
            for o_idx in sgl_o:
                if not _has_h(mol, o_idx):
                    # Check that O is also bonded to a C (ester, not anhydride check)
                    o_c_nbrs = _single_bonded_to(mol, o_idx, "C")
                    other_c = [c for c in o_c_nbrs if c != idx]
                    if other_c:
                        found.append(FunctionalGroup(
                            name="ester",
                            smarts_like="[CX3](=O)[O][C]",
                            atoms=[idx, dbl_o[0], o_idx, other_c[0]],
                            center=idx,
                        ))
                        break
    return found


def _detect_amides(mol: Molecule) -> list[FunctionalGroup]:
    """Amide: C(=O)-N."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        sgl_n = _single_bonded_to(mol, idx, "N")
        if dbl_o and sgl_n:
            found.append(FunctionalGroup(
                name="amide", smarts_like="[CX3](=O)[NX3]",
                atoms=[idx, dbl_o[0], sgl_n[0]], center=idx,
            ))
    return found


def _detect_amines(mol: Molecule) -> list[FunctionalGroup]:
    """Primary, secondary, and tertiary amines (not amides)."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "N":
            continue
        # Skip if N is part of an amide (bonded to a carbonyl C)
        is_amide = False
        for c_idx in _single_bonded_to(mol, idx, "C"):
            if _double_bonded_to(mol, c_idx, "O"):
                is_amide = True
                break
        if is_amide:
            continue

        c_count = _count_element_neighbors(mol, idx, "C")
        h_count = _h_count(mol, idx)

        if c_count == 1 and h_count == 2:
            found.append(FunctionalGroup(
                name="primary_amine", smarts_like="[NX3H2][C]",
                atoms=[idx] + _single_bonded_to(mol, idx, "C"),
                center=idx,
            ))
        elif c_count == 2 and h_count == 1:
            found.append(FunctionalGroup(
                name="secondary_amine", smarts_like="[NX3H1]([C])[C]",
                atoms=[idx] + _single_bonded_to(mol, idx, "C"),
                center=idx,
            ))
        elif c_count == 3 and h_count == 0:
            found.append(FunctionalGroup(
                name="tertiary_amine", smarts_like="[NX3]([C])([C])[C]",
                atoms=[idx] + _single_bonded_to(mol, idx, "C"),
                center=idx,
            ))
    return found


def _detect_alkyl_halides(mol: Molecule) -> list[FunctionalGroup]:
    """Alkyl halide: C bonded to F, Cl, Br, or I."""
    halogens = {"F", "Cl", "Br", "I"}
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        for n in _neighbors(mol, idx):
            if _element(mol, n) in halogens and _bond_order(mol, idx, n) == 1.0:
                hal = _element(mol, n)
                found.append(FunctionalGroup(
                    name=f"alkyl_halide_{hal.lower()}",
                    smarts_like=f"[C][{hal}]",
                    atoms=[idx, n], center=idx,
                ))
    return found


def _detect_alkenes(mol: Molecule) -> list[FunctionalGroup]:
    """Alkene: C=C double bond."""
    found: list[FunctionalGroup] = []
    seen: set[tuple[int, int]] = set()
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        for n in _double_bonded_to(mol, idx, "C"):
            pair = (min(idx, n), max(idx, n))
            if pair not in seen:
                seen.add(pair)
                found.append(FunctionalGroup(
                    name="alkene", smarts_like="[C]=[C]",
                    atoms=list(pair), center=pair[0],
                ))
    return found


def _detect_alkynes(mol: Molecule) -> list[FunctionalGroup]:
    """Alkyne: C#C triple bond."""
    found: list[FunctionalGroup] = []
    seen: set[tuple[int, int]] = set()
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        for n in _triple_bonded_to(mol, idx, "C"):
            pair = (min(idx, n), max(idx, n))
            if pair not in seen:
                seen.add(pair)
                found.append(FunctionalGroup(
                    name="alkyne", smarts_like="[C]#[C]",
                    atoms=list(pair), center=pair[0],
                ))
    return found


def _detect_ethers(mol: Molecule) -> list[FunctionalGroup]:
    """Ether: C-O-C (oxygen single-bonded to two carbons, no C=O on either)."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "O":
            continue
        c_nbrs = _single_bonded_to(mol, idx, "C")
        if len(c_nbrs) != 2:
            continue
        # Exclude if either C has a C=O (that would be ester)
        is_ester = any(_double_bonded_to(mol, c, "O") for c in c_nbrs)
        if is_ester:
            continue
        found.append(FunctionalGroup(
            name="ether", smarts_like="[C]-[O]-[C]",
            atoms=[c_nbrs[0], idx, c_nbrs[1]], center=idx,
        ))
    return found


def _detect_thiols(mol: Molecule) -> list[FunctionalGroup]:
    """Thiol: S bonded to C with an H (explicit or implicit)."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "S":
            continue
        nbrs = _neighbors(mol, idx)
        elems = [_element(mol, n) for n in nbrs]
        c_indices = [nbrs[i] for i, e in enumerate(elems) if e == "C"]
        if not c_indices:
            continue
        if _has_h(mol, idx):
            c_idx = c_indices[0]
            found.append(FunctionalGroup(
                name="thiol", smarts_like="[C]-[SH]",
                atoms=[c_idx, idx], center=idx,
            ))
    return found


def _detect_nitriles(mol: Molecule) -> list[FunctionalGroup]:
    """Nitrile: C#N triple bond."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        for n in _triple_bonded_to(mol, idx, "N"):
            found.append(FunctionalGroup(
                name="nitrile", smarts_like="[C]#[N]",
                atoms=[idx, n], center=idx,
            ))
    return found


def _detect_nitro(mol: Molecule) -> list[FunctionalGroup]:
    """Nitro group: N bonded to two O atoms with at least one N=O."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "N":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        sgl_o = _single_bonded_to(mol, idx, "O")
        total_o = len(dbl_o) + len(sgl_o)
        if total_o >= 2 and len(dbl_o) >= 1:
            found.append(FunctionalGroup(
                name="nitro", smarts_like="[N](=O)[O]",
                atoms=[idx] + dbl_o + sgl_o, center=idx,
            ))
    return found


def _detect_aromatic_rings(mol: Molecule) -> list[FunctionalGroup]:
    """Simplified aromatic ring detection.

    Looks for six-membered rings composed entirely of carbons (or with
    one nitrogen for pyridine) where all ring bonds have order >= 1.5
    (aromatic) *or* alternating single/double bonds that form a
    conjugated cycle.

    This uses a breadth-first ring-finding approach limited to length 6.
    """
    found: list[FunctionalGroup] = []
    n_atoms = len(mol.atoms)
    seen_rings: set[tuple[int, ...]] = set()

    for start in range(n_atoms):
        if _element(mol, start) not in ("C", "N", "O", "S"):
            continue
        # BFS / DFS for 6-membered rings from start
        rings = _find_rings_of_size(mol, start, 6)
        # Also search for 5-membered rings (furan, thiophene, pyrrole, etc.)
        rings += _find_rings_of_size(mol, start, 5)
        for ring in rings:
            canon = _canonicalise_ring(ring)
            if canon in seen_rings:
                continue
            seen_rings.add(canon)
            # Check that ring is plausibly aromatic
            if _ring_is_aromatic(mol, ring):
                found.append(FunctionalGroup(
                    name="aromatic_ring",
                    smarts_like="c1ccccc1",
                    atoms=list(ring), center=ring[0],
                ))
    return found


def _find_rings_of_size(mol: Molecule, start: int, size: int) -> list[tuple[int, ...]]:
    """Return all simple rings of exactly *size* atoms that include *start*.

    Uses iterative depth-limited DFS.  To keep cost manageable the
    search only proceeds through C and N atoms.
    """
    results: list[tuple[int, ...]] = []
    allowed = {"C", "N", "O", "S"}
    # stack entries: (current_atom, path_so_far)
    stack: list[tuple[int, list[int]]] = [(start, [start])]
    while stack:
        current, path = stack.pop()
        if len(path) == size:
            # Check if we can close the ring back to start
            if start in [n for n in _neighbors(mol, current)]:
                results.append(tuple(path))
            continue
        for nbr in _neighbors(mol, current):
            if nbr == start and len(path) >= 3:
                # Early closure -- ring smaller than *size*; skip
                continue
            if nbr in path:
                continue
            if _element(mol, nbr) not in allowed:
                continue
            stack.append((nbr, path + [nbr]))
    return results


def _canonicalise_ring(ring: tuple[int, ...]) -> tuple[int, ...]:
    """Return a canonical form for a ring so that rotations/reflections match."""
    min_val = min(ring)
    min_idx = ring.index(min_val)
    forward = ring[min_idx:] + ring[:min_idx]
    backward = (ring[min_idx],) + tuple(reversed(ring[:min_idx])) + tuple(reversed(ring[min_idx + 1:]))
    return min(forward, backward)


def _ring_is_aromatic(mol: Molecule, ring: tuple[int, ...]) -> bool:
    """Heuristically decide if a ring is aromatic.

    A ring is considered aromatic if any of:
    - All bond orders are >= 1.5 (explicit aromatic annotation), **or**
    - The ring consists of alternating single (1.0) and double (2.0)
      bonds forming a fully conjugated system, **or**
    - All ring atoms have SP2 hybridization (aromatic SMILES atoms are
      assigned SP2 by the parser even though bonds are stored as order 1;
      this catches furan, thiophene, pyrrole and other heteroaromatics).
    """
    n = len(ring)
    orders = []
    for i in range(n):
        a, b = ring[i], ring[(i + 1) % n]
        orders.append(_bond_order(mol, a, b))

    # All aromatic-annotated bonds
    if all(o >= 1.5 for o in orders):
        return True

    # Alternating single/double
    if all(o in (1.0, 2.0) for o in orders):
        alternating = all(orders[i] != orders[(i + 1) % n] for i in range(n))
        if alternating:
            return True

    # All atoms SP2-hybridized (aromatic SMILES atoms, or conjugated rings)
    if all(mol.atoms[idx].hybridization == Hybridization.SP2 for idx in ring):
        return True

    return False


def _detect_epoxides(mol: Molecule) -> list[FunctionalGroup]:
    """Epoxide: 3-membered ring containing one O and two C atoms."""
    found: list[FunctionalGroup] = []
    seen_rings: set[tuple[int, ...]] = set()

    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "O":
            continue
        c_nbrs = _single_bonded_to(mol, idx, "C")
        if len(c_nbrs) < 2:
            continue
        # Check each pair of C neighbours for a bond between them
        for i in range(len(c_nbrs)):
            for j in range(i + 1, len(c_nbrs)):
                c1, c2 = c_nbrs[i], c_nbrs[j]
                if _bond_order(mol, c1, c2) > 0:
                    canon = _canonicalise_ring((idx, c1, c2))
                    if canon not in seen_rings:
                        seen_rings.add(canon)
                        found.append(FunctionalGroup(
                            name="epoxide",
                            smarts_like="C1OC1",
                            atoms=[c1, idx, c2], center=idx,
                        ))
    return found


def _detect_acid_chlorides(mol: Molecule) -> list[FunctionalGroup]:
    """Acid chloride (acyl chloride): C(=O)Cl."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        sgl_cl = _single_bonded_to(mol, idx, "Cl")
        if dbl_o and sgl_cl:
            found.append(FunctionalGroup(
                name="acid_chloride", smarts_like="[CX3](=O)[Cl]",
                atoms=[idx, dbl_o[0], sgl_cl[0]], center=idx,
            ))
    return found


def _detect_anhydrides(mol: Molecule) -> list[FunctionalGroup]:
    """Acid anhydride: C(=O)-O-C(=O)."""
    found: list[FunctionalGroup] = []
    seen: set[int] = set()
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "O":
            continue
        if idx in seen:
            continue
        c_nbrs = _single_bonded_to(mol, idx, "C")
        if len(c_nbrs) < 2:
            continue
        # Both C neighbours must have a C=O
        carbonyl_cs = [c for c in c_nbrs if _double_bonded_to(mol, c, "O")]
        if len(carbonyl_cs) >= 2:
            c1, c2 = carbonyl_cs[0], carbonyl_cs[1]
            o1 = _double_bonded_to(mol, c1, "O")[0]
            o2 = _double_bonded_to(mol, c2, "O")[0]
            seen.add(idx)
            found.append(FunctionalGroup(
                name="anhydride", smarts_like="[CX3](=O)[O][CX3](=O)",
                atoms=[c1, o1, idx, c2, o2], center=idx,
            ))
    return found


def _detect_sulfoxides(mol: Molecule) -> list[FunctionalGroup]:
    """Sulfoxide: S(=O) bonded to two carbons (no second O=S)."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "S":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        c_nbrs = _single_bonded_to(mol, idx, "C")
        if len(dbl_o) == 1 and len(c_nbrs) >= 2:
            found.append(FunctionalGroup(
                name="sulfoxide", smarts_like="[SX3](=O)([C])[C]",
                atoms=[idx, dbl_o[0]] + c_nbrs[:2], center=idx,
            ))
    return found


def _detect_sulfones(mol: Molecule) -> list[FunctionalGroup]:
    """Sulfone: S(=O)(=O) bonded to two carbons."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "S":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        c_nbrs = _single_bonded_to(mol, idx, "C")
        if len(dbl_o) >= 2 and len(c_nbrs) >= 2:
            found.append(FunctionalGroup(
                name="sulfone", smarts_like="[SX4](=O)(=O)([C])[C]",
                atoms=[idx] + dbl_o[:2] + c_nbrs[:2], center=idx,
            ))
    return found


def _detect_imines(mol: Molecule) -> list[FunctionalGroup]:
    """Imine: C=N (not part of nitrile C#N)."""
    found: list[FunctionalGroup] = []
    seen: set[tuple[int, int]] = set()
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "C":
            continue
        for n in _neighbors(mol, idx):
            if _element(mol, n) != "N":
                continue
            if _bond_order(mol, idx, n) != 2.0:
                continue
            pair = (min(idx, n), max(idx, n))
            if pair in seen:
                continue
            seen.add(pair)
            found.append(FunctionalGroup(
                name="imine", smarts_like="[C]=[N]",
                atoms=list(pair), center=idx,
            ))
    return found


def _detect_boronic_acids(mol: Molecule) -> list[FunctionalGroup]:
    """Boronic acid: B bonded to two O atoms (and one C)."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "B":
            continue
        o_nbrs = _single_bonded_to(mol, idx, "O")
        c_nbrs = _single_bonded_to(mol, idx, "C")
        if len(o_nbrs) >= 2 and len(c_nbrs) >= 1:
            found.append(FunctionalGroup(
                name="boronic_acid", smarts_like="[B]([OH])([OH])[C]",
                atoms=[idx] + o_nbrs[:2] + c_nbrs[:1], center=idx,
            ))
    return found


def _detect_phosphonates(mol: Molecule) -> list[FunctionalGroup]:
    """Phosphonate ester: P with one P=O and two P-O-C linkages."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "P":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        sgl_o = _single_bonded_to(mol, idx, "O")
        if len(dbl_o) >= 1 and len(sgl_o) >= 2:
            found.append(FunctionalGroup(
                name="phosphonate", smarts_like="[P](=O)([O])([O])",
                atoms=[idx] + dbl_o[:1] + sgl_o[:2], center=idx,
            ))
    return found


def _detect_sulfonamides(mol: Molecule) -> list[FunctionalGroup]:
    """Sulfonamide: S(=O)(=O)-N."""
    found: list[FunctionalGroup] = []
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol != "S":
            continue
        dbl_o = _double_bonded_to(mol, idx, "O")
        sgl_n = _single_bonded_to(mol, idx, "N")
        if len(dbl_o) >= 2 and len(sgl_n) >= 1:
            found.append(FunctionalGroup(
                name="sulfonamide", smarts_like="[S](=O)(=O)[N]",
                atoms=[idx] + dbl_o[:2] + sgl_n[:1], center=idx,
            ))
    return found
