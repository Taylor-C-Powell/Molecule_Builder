"""SMARTS string tokenizer.

Converts a SMARTS string into a sequence of SmartsToken objects for
the parser.  SMARTS extends SMILES with additional atom primitives
(atomic number, degree, valence, H-count, ring membership, ring size,
aromatic/aliphatic), logical operators (NOT, AND, OR), bond primitives
(any-bond, ring-bond), wildcards, and recursive SMARTS.

Handles:
- Bracket atoms with SMARTS primitives: [#6], [CD3], [!#7], [#6,#7]
- Wildcard atom: *
- Organic subset atoms: B, C, N, O, P, S, F, Cl, Br, I
- Aromatic atoms: b, c, n, o, p, s
- Bond types: - (single), = (double), # (triple), : (aromatic),
  ~ (any bond), @ (ring bond)
- Branch notation: ( and )
- Ring closure digits: 0-9, and %nn
- Dot disconnection: .
- Recursive SMARTS: $(...) patterns
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto


# ===================================================================
# Token types and data class
# ===================================================================

class SmartsTokenType(Enum):
    """Token categories for SMARTS lexing."""
    ATOM = auto()           # atom specification (organic, bracket, or wildcard)
    BOND = auto()           # -, =, #, :, ~, @
    BRANCH_OPEN = auto()    # (
    BRANCH_CLOSE = auto()   # )
    RING_DIGIT = auto()     # ring closure number
    DOT = auto()            # . (disconnection)
    WILDCARD = auto()       # * (any atom)


@dataclass
class SmartsToken:
    """A single lexical token from a SMARTS string.

    Attributes
    ----------
    type : SmartsTokenType
        The category of this token.
    value : str
        Raw text of the token or the atom symbol.
    atomic_number : int | None
        Atomic number constraint from ``#N`` primitive.
    degree : int | None
        Degree (heavy-atom connectivity) constraint from ``D`` or ``Dn``.
    valence : int | None
        Total valence constraint from ``v`` or ``vn``.
    hcount : int | None
        Hydrogen count constraint from ``H`` or ``Hn``.
    charge : int
        Formal charge constraint.
    aromatic : bool | None
        True if aromatic, False if aliphatic, None if unspecified.
    ring_membership : int | None
        Ring membership count from ``R`` or ``Rn``.
    ring_size : int | None
        Smallest ring size from ``r`` or ``rn``.
    is_not : bool
        True if the atom specification is negated with ``!``.
    symbol : str | None
        Element symbol extracted from the atom specification.
    recursive : str | None
        Recursive SMARTS string from ``$(...)`` pattern.
    or_specs : list | None
        List of SmartsToken for OR (,) combinations within brackets.
    and_specs : list | None
        List of SmartsToken for AND (&/;) combinations within brackets.
    """
    type: SmartsTokenType
    value: str
    atomic_number: int | None = None
    degree: int | None = None
    valence: int | None = None
    hcount: int | None = None
    charge: int = 0
    aromatic: bool | None = None
    ring_membership: int | None = None
    ring_size: int | None = None
    is_not: bool = False
    symbol: str | None = None
    recursive: str | None = None
    or_specs: list | None = None
    and_specs: list | None = None


# ===================================================================
# Constants
# ===================================================================

ORGANIC_SUBSET = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"}
AROMATIC_ATOMS = {"b", "c", "n", "o", "p", "s"}
TWO_LETTER_ORGANIC = {"Cl", "Br"}

SMARTS_BOND_CHARS = {"-", "=", "#", ":", "~", "@"}


# ===================================================================
# Bracket-atom parser for SMARTS primitives
# ===================================================================

def _parse_smarts_bracket(smarts: str, start: int) -> tuple[SmartsToken, int]:
    """Parse a SMARTS bracket atom ``[...]`` starting at *start* (the ``[``).

    Handles SMARTS-specific primitives in the following order:
    1. ``!`` (NOT prefix)
    2. Element symbol, ``*`` (wildcard), ``a``/``A`` (aromatic/aliphatic)
    3. ``#N`` (atomic number)
    4. ``D``, ``Dn`` (degree)
    5. ``v``, ``vn`` (valence)
    6. ``H``, ``Hn`` (H-count)
    7. ``R``, ``Rn`` (ring membership)
    8. ``r``, ``rn`` (ring size)
    9. ``+``, ``-``, ``+n``, ``-n`` (charge)
    10. ``,`` (OR), ``&`` (AND high), ``;`` (AND low)
    11. ``$(...)`` (recursive SMARTS)

    Returns the SmartsToken and the index after the closing ``]``.
    """
    pos = start + 1  # skip '['
    end = _find_bracket_close(smarts, pos)
    inner = smarts[pos:end]

    # Check for OR (,) or AND (&, ;) at the top level of the bracket
    or_parts = _split_top_level(inner, ",")
    if len(or_parts) > 1:
        # Parse each OR alternative as a separate spec
        or_specs = []
        for part in or_parts:
            sub_tok = _parse_single_spec("[" + part + "]", 0)
            or_specs.append(sub_tok)
        return SmartsToken(
            type=SmartsTokenType.ATOM,
            value=inner,
            or_specs=or_specs,
        ), end + 1

    and_parts = _split_top_level(inner, ";")
    if len(and_parts) > 1:
        and_specs = []
        for part in and_parts:
            sub_tok = _parse_single_spec("[" + part + "]", 0)
            and_specs.append(sub_tok)
        return SmartsToken(
            type=SmartsTokenType.ATOM,
            value=inner,
            and_specs=and_specs,
        ), end + 1

    # Also handle & as high-precedence AND within a single bracket
    and_hi_parts = _split_top_level(inner, "&")
    if len(and_hi_parts) > 1:
        and_specs = []
        for part in and_hi_parts:
            sub_tok = _parse_single_spec("[" + part + "]", 0)
            and_specs.append(sub_tok)
        return SmartsToken(
            type=SmartsTokenType.ATOM,
            value=inner,
            and_specs=and_specs,
        ), end + 1

    tok = _parse_single_spec(smarts, start)
    return tok, end + 1


def _find_bracket_close(smarts: str, start: int) -> int:
    """Find the matching `]` for a bracket atom, handling nested $(...).

    Parameters
    ----------
    start : int
        Index just after the opening ``[``.
    """
    depth = 0
    i = start
    n = len(smarts)
    while i < n:
        ch = smarts[i]
        if ch == "$" and i + 1 < n and smarts[i + 1] == "(":
            # Start of recursive SMARTS -- skip balanced parens
            depth += 1
            i += 2
            continue
        if depth > 0:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            i += 1
            continue
        if ch == "]":
            return i
        i += 1
    raise ValueError(f"Unclosed bracket atom starting near position {start - 1}")


def _split_top_level(inner: str, sep: str) -> list[str]:
    """Split *inner* on *sep* but only at the top level (not inside $(...))."""
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    i = 0
    n = len(inner)
    while i < n:
        ch = inner[i]
        if ch == "$" and i + 1 < n and inner[i + 1] == "(":
            depth += 1
            current.append(ch)
            i += 1
            current.append(inner[i])
            i += 1
            continue
        if depth > 0:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            current.append(ch)
            i += 1
            continue
        if ch == sep:
            parts.append("".join(current))
            current = []
            i += 1
            continue
        current.append(ch)
        i += 1
    parts.append("".join(current))
    return parts


def _parse_single_spec(smarts: str, start: int) -> SmartsToken:
    """Parse a single (non-composite) SMARTS bracket atom spec."""
    pos = start + 1  # skip '['
    end = _find_bracket_close(smarts, pos)
    inner = smarts[pos:end]

    is_not = False
    symbol: str | None = None
    atomic_number: int | None = None
    degree: int | None = None
    valence_val: int | None = None
    hcount: int | None = None
    charge = 0
    aromatic: bool | None = None
    ring_membership: int | None = None
    ring_size: int | None = None
    recursive: str | None = None

    i = 0
    n = len(inner)

    # 1. NOT prefix
    if i < n and inner[i] == "!":
        is_not = True
        i += 1

    # 2. Element symbol, *, a/A
    if i < n and inner[i] == "*":
        symbol = None  # wildcard
        i += 1
    elif i < n and inner[i] == "a":
        # Lowercase 'a' alone = aromatic any atom
        # But could also be start of element like "as" (arsenic aromatic)
        if i + 1 < n and inner[i + 1].islower() and inner[i + 1] not in "0123456789":
            # Check for aromatic element symbols: as, se
            if i + 1 < n and inner[i:i + 2] == "as":
                symbol = "As"
                aromatic = True
                i += 2
            elif i + 1 < n and inner[i:i + 2] == "se":
                symbol = "Se"
                aromatic = True
                i += 2
            else:
                # Just aromatic 'a' -- generic aromatic
                aromatic = True
                i += 1
        else:
            # Just 'a' -- generic aromatic atom
            aromatic = True
            i += 1
    elif i < n and inner[i] == "A":
        # Could be 'A' (aliphatic any) or start of element like Al, Ar
        if i + 1 < n and inner[i + 1].islower():
            # Try two-letter element
            maybe_sym = inner[i:i + 2]
            if maybe_sym in ("Al", "Ar", "As", "Au", "Ag", "At"):
                symbol = maybe_sym
                aromatic = False
                i += 2
            else:
                # Generic aliphatic
                aromatic = False
                i += 1
        else:
            aromatic = False
            i += 1
    elif i < n and inner[i] == "#":
        # Atomic number -- symbol stays None
        pass
    elif i < n and inner[i].isupper():
        symbol = inner[i]
        i += 1
        # Check for two-letter elements, but NOT if the second letter
        # is a SMARTS primitive (D, v, H, R, r, h, x) that should be
        # parsed as a constraint on the single-letter element.
        _smarts_primitives = {"D", "v", "H", "R", "r", "h", "x"}
        if i < n and inner[i].islower() and inner[i] not in _smarts_primitives:
            candidate = symbol + inner[i]
            # Common two-letter elements
            _two_letter = {
                "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Cl",
                "Ar", "Ca", "Sc", "Ti", "Mn", "Fe", "Co", "Ni",
                "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
                "Sr", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
                "Cd", "In", "Sn", "Sb", "Te", "Xe", "Cs", "Ba", "La",
                "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
                "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "Re", "Os",
                "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
                "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "Np", "Pu", "Am",
                "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
            }
            if candidate in _two_letter:
                symbol = candidate
                i += 1
        aromatic = False
    elif i < n and inner[i].islower():
        # Aromatic atom: c, n, o, s, p, etc.
        symbol = inner[i].upper()
        aromatic = True
        i += 1

    # 3. Atomic number #N
    if i < n and inner[i] == "#":
        i += 1
        num_start = i
        while i < n and inner[i].isdigit():
            i += 1
        if i > num_start:
            atomic_number = int(inner[num_start:i])

    # 4. Degree D, Dn
    if i < n and inner[i] == "D":
        i += 1
        num_start = i
        while i < n and inner[i].isdigit():
            i += 1
        if i > num_start:
            degree = int(inner[num_start:i])
        else:
            degree = 1  # bare D means degree 1

    # 5. Valence v, vn
    if i < n and inner[i] == "v":
        i += 1
        num_start = i
        while i < n and inner[i].isdigit():
            i += 1
        if i > num_start:
            valence_val = int(inner[num_start:i])
        else:
            valence_val = 1  # bare v means valence 1

    # 6. H-count H, Hn
    if i < n and inner[i] == "H":
        i += 1
        num_start = i
        while i < n and inner[i].isdigit():
            i += 1
        if i > num_start:
            hcount = int(inner[num_start:i])
        else:
            hcount = 1  # bare H means 1 hydrogen

    # 7. Ring membership R, Rn
    if i < n and inner[i] == "R":
        i += 1
        num_start = i
        while i < n and inner[i].isdigit():
            i += 1
        if i > num_start:
            ring_membership = int(inner[num_start:i])
        else:
            ring_membership = -1  # bare R means "in at least one ring"

    # 8. Ring size r, rn
    if i < n and inner[i] == "r":
        i += 1
        num_start = i
        while i < n and inner[i].isdigit():
            i += 1
        if i > num_start:
            ring_size = int(inner[num_start:i])
        else:
            ring_size = -1  # bare r means "in at least one ring"

    # 9. Charge +, -, +n, -n
    if i < n and inner[i] in ("+", "-"):
        sign = 1 if inner[i] == "+" else -1
        i += 1
        ch_start = i
        while i < n and inner[i].isdigit():
            i += 1
        if i > ch_start:
            charge = sign * int(inner[ch_start:i])
        else:
            # Count consecutive +/- signs
            extra = 0
            ch_char = "+" if sign == 1 else "-"
            while i < n and inner[i] == ch_char:
                extra += 1
                i += 1
            charge = sign * (1 + extra)

    # 10. Recursive SMARTS $(...)
    if i < n and inner[i] == "$" and i + 1 < n and inner[i + 1] == "(":
        i += 2  # skip $(
        depth = 1
        rec_start = i
        while i < n and depth > 0:
            if inner[i] == "(":
                depth += 1
            elif inner[i] == ")":
                depth -= 1
            if depth > 0:
                i += 1
        recursive = inner[rec_start:i]
        i += 1  # skip closing )

    return SmartsToken(
        type=SmartsTokenType.ATOM,
        value=inner,
        atomic_number=atomic_number,
        degree=degree,
        valence=valence_val,
        hcount=hcount,
        charge=charge,
        aromatic=aromatic,
        ring_membership=ring_membership,
        ring_size=ring_size,
        is_not=is_not,
        symbol=symbol,
        recursive=recursive,
    )


# ===================================================================
# Main tokenizer
# ===================================================================

def tokenize_smarts(smarts: str) -> list[SmartsToken]:
    """Convert a SMARTS string into a list of SmartsToken objects.

    Parameters
    ----------
    smarts : str
        A valid SMARTS string, e.g. ``"[#6]"``, ``"[OH]"``,
        ``"c1ccccc1"``, ``"[CD3](=O)[OD1]"``.

    Returns
    -------
    list[SmartsToken]
        Ordered sequence of tokens ready for the SMARTS parser.

    Raises
    ------
    ValueError
        If the string contains unexpected characters or malformed
        bracket atoms.
    """
    tokens: list[SmartsToken] = []
    i = 0
    n = len(smarts)

    while i < n:
        ch = smarts[i]

        # --- bracket atom ---
        if ch == "[":
            token, i = _parse_smarts_bracket(smarts, i)
            tokens.append(token)
            continue

        # --- wildcard ---
        if ch == "*":
            tokens.append(SmartsToken(
                type=SmartsTokenType.WILDCARD,
                value="*",
                symbol=None,
            ))
            i += 1
            continue

        # --- two-letter organic atoms (Cl, Br) ---
        if i + 1 < n and smarts[i:i + 2] in TWO_LETTER_ORGANIC:
            sym = smarts[i:i + 2]
            tokens.append(SmartsToken(
                type=SmartsTokenType.ATOM,
                value=sym,
                symbol=sym,
                aromatic=False,
            ))
            i += 2
            continue

        # --- single-letter organic subset atoms ---
        if ch in {"B", "C", "N", "O", "P", "S", "F", "I"}:
            tokens.append(SmartsToken(
                type=SmartsTokenType.ATOM,
                value=ch,
                symbol=ch,
                aromatic=False,
            ))
            i += 1
            continue

        # --- aromatic atoms ---
        if ch in {"b", "c", "n", "o", "p", "s"}:
            tokens.append(SmartsToken(
                type=SmartsTokenType.ATOM,
                value=ch,
                symbol=ch.upper(),
                aromatic=True,
            ))
            i += 1
            continue

        # --- bond symbols (SMARTS extends SMILES with ~ and @) ---
        if ch in SMARTS_BOND_CHARS:
            tokens.append(SmartsToken(
                type=SmartsTokenType.BOND,
                value=ch,
            ))
            i += 1
            continue

        # --- branch open / close ---
        if ch == "(":
            # Check for recursive SMARTS $(...) -- should already be
            # consumed by bracket parser; this handles bare parentheses
            tokens.append(SmartsToken(
                type=SmartsTokenType.BRANCH_OPEN,
                value="(",
            ))
            i += 1
            continue

        if ch == ")":
            tokens.append(SmartsToken(
                type=SmartsTokenType.BRANCH_CLOSE,
                value=")",
            ))
            i += 1
            continue

        # --- ring closure digit ---
        if ch.isdigit():
            tokens.append(SmartsToken(
                type=SmartsTokenType.RING_DIGIT,
                value=ch,
            ))
            i += 1
            continue

        # --- two-digit ring closure: %nn ---
        if ch == "%":
            if i + 2 >= n or not smarts[i + 1].isdigit() or not smarts[i + 2].isdigit():
                raise ValueError(
                    f"Expected two digits after '%' at position {i}")
            tokens.append(SmartsToken(
                type=SmartsTokenType.RING_DIGIT,
                value=smarts[i + 1:i + 3],
            ))
            i += 3
            continue

        # --- dot disconnection ---
        if ch == ".":
            tokens.append(SmartsToken(
                type=SmartsTokenType.DOT,
                value=".",
            ))
            i += 1
            continue

        # --- unexpected character ---
        raise ValueError(
            f"Unexpected character {ch!r} at position {i} in SMARTS "
            f"string {smarts!r}")

    return tokens
