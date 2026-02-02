"""SMILES string tokenizer.

Converts a SMILES string into a sequence of Token objects for the parser.

Handles:
- Organic subset atoms: B, C, N, O, P, S, F, Cl, Br, I (uppercase)
- Bracket atoms: [NH3+], [Fe], [13C], etc.
- Aromatic atoms: b, c, n, o, p, s (lowercase)
- Bond types: - (single), = (double), # (triple), : (aromatic)
- Branch notation: ( and )
- Ring closure digits: 0-9, and %nn for two-digit ring numbers
- Dot disconnection: .
- Hydrogen counts in brackets: [NH2], [CH3]
- Charges in brackets: +, -, +2, -1
- Isotopes in brackets: [13C], [2H]
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto


# ===================================================================
# Token types and data class
# ===================================================================

class TokenType(Enum):
    ATOM = auto()           # organic subset atom or bracket atom
    BOND = auto()           # -, =, #, :
    BRANCH_OPEN = auto()    # (
    BRANCH_CLOSE = auto()   # )
    RING_DIGIT = auto()     # ring closure number
    DOT = auto()            # . (disconnection)


@dataclass
class Token:
    """A single lexical token from a SMILES string.

    Attributes
    ----------
    type : TokenType
        The category of this token.
    value : str
        The atom symbol for ATOM tokens, the bond character for BOND
        tokens, or the digit string for RING_DIGIT tokens.
    isotope : int | None
        Mass number from a bracket atom, e.g. 13 in ``[13C]``.
    hcount : int | None
        Explicit hydrogen count from a bracket atom, e.g. 2 in ``[NH2]``.
        ``None`` means no explicit H specification (implicit semantics
        apply).  ``0`` means explicitly zero hydrogens.
    charge : int
        Formal charge from a bracket atom, e.g. +1 in ``[NH4+]``.
    aromatic : bool
        True if the atom was given in lowercase (aromatic notation).
    chirality : str | None
        Chirality marker from a bracket atom: ``"@"`` (anticlockwise),
        ``"@@"`` (clockwise), or ``None`` (no chirality specified).
    """
    type: TokenType
    value: str
    isotope: int | None = None
    hcount: int | None = None
    charge: int = 0
    aromatic: bool = False
    chirality: str | None = None


# ===================================================================
# Constants
# ===================================================================

# Organic subset: atoms that don't need brackets
ORGANIC_SUBSET = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"}
AROMATIC_ATOMS = {"b", "c", "n", "o", "p", "s"}

# Two-letter organic subset atoms (checked before single-letter)
TWO_LETTER_ORGANIC = {"Cl", "Br"}

# Default valence for implicit hydrogen calculation
DEFAULT_VALENCE: dict[str, list[int]] = {
    "B": [3], "C": [4], "N": [3, 5], "O": [2], "P": [3, 5],
    "S": [2, 4, 6], "F": [1], "Cl": [1], "Br": [1], "I": [1],
    "b": [3], "c": [4], "n": [3], "o": [2], "p": [3], "s": [2],
}

BOND_CHARS = {"-", "=", "#", ":", "/", "\\"}


# ===================================================================
# Bracket-atom parser helper
# ===================================================================

def _parse_bracket(smiles: str, start: int) -> tuple[Token, int]:
    """Parse a bracket atom ``[...]`` starting at *start* (the ``[``).

    Returns the Token and the index of the character after the closing
    ``]``.
    """
    pos = start + 1  # skip '['
    end = smiles.index("]", pos)
    inner = smiles[pos:end]

    isotope: int | None = None
    symbol: str = ""
    hcount: int | None = None
    charge: int = 0
    aromatic: bool = False
    chirality: str | None = None

    i = 0
    n = len(inner)

    # --- isotope (leading digits) ---
    iso_start = i
    while i < n and inner[i].isdigit():
        i += 1
    if i > iso_start:
        isotope = int(inner[iso_start:i])

    # --- element symbol ---
    # Symbol starts with an uppercase letter followed by optional lowercase,
    # OR a single lowercase letter for aromatic atoms.
    if i < n and inner[i].isupper():
        symbol = inner[i]
        i += 1
        while i < n and inner[i].islower():
            symbol += inner[i]
            i += 1
    elif i < n and inner[i].islower():
        # aromatic bracket atom
        symbol = inner[i]
        aromatic = True
        i += 1
    else:
        raise ValueError(
            f"Expected element symbol in bracket atom: [{inner}]")

    # --- chirality (@, @@) ---
    if i < n and inner[i] == "@":
        i += 1
        if i < n and inner[i] == "@":
            chirality = "@@"
            i += 1
        else:
            chirality = "@"

    # --- hydrogen count ---
    if i < n and inner[i] == "H":
        i += 1
        h_start = i
        while i < n and inner[i].isdigit():
            i += 1
        if i > h_start:
            hcount = int(inner[h_start:i])
        else:
            hcount = 1  # bare H means 1 hydrogen

    # --- charge ---
    if i < n and inner[i] in ("+", "-"):
        sign = 1 if inner[i] == "+" else -1
        i += 1
        ch_start = i
        while i < n and inner[i].isdigit():
            i += 1
        if i > ch_start:
            charge = sign * int(inner[ch_start:i])
        else:
            # Count consecutive +/- signs (e.g. ++ means +2)
            extra = 0
            ch = "+" if sign == 1 else "-"
            while i < n and inner[i] == ch:
                extra += 1
                i += 1
            charge = sign * (1 + extra)

    return Token(
        type=TokenType.ATOM,
        value=symbol,
        isotope=isotope,
        hcount=hcount,
        charge=charge,
        aromatic=aromatic,
        chirality=chirality,
    ), end + 1


# ===================================================================
# Main tokenizer
# ===================================================================

def tokenize(smiles: str) -> list[Token]:
    """Convert a SMILES string into a list of Token objects.

    Parameters
    ----------
    smiles : str
        A valid SMILES string, e.g. ``"CCO"``, ``"c1ccccc1"``,
        ``"[NH4+]"``.

    Returns
    -------
    list[Token]
        Ordered sequence of tokens ready for the parser.

    Raises
    ------
    ValueError
        If the string contains unexpected characters or malformed
        bracket atoms.
    """
    tokens: list[Token] = []
    i = 0
    n = len(smiles)

    while i < n:
        ch = smiles[i]

        # --- bracket atom ---
        if ch == "[":
            token, i = _parse_bracket(smiles, i)
            tokens.append(token)
            continue

        # --- two-letter organic atoms (Cl, Br) ---
        if i + 1 < n and smiles[i:i + 2] in TWO_LETTER_ORGANIC:
            tokens.append(Token(type=TokenType.ATOM, value=smiles[i:i + 2]))
            i += 2
            continue

        # --- single-letter organic subset atoms ---
        if ch in {"B", "C", "N", "O", "P", "S", "F", "I"}:
            tokens.append(Token(type=TokenType.ATOM, value=ch))
            i += 1
            continue

        # --- aromatic atoms ---
        if ch in {"b", "c", "n", "o", "p", "s"}:
            tokens.append(Token(
                type=TokenType.ATOM, value=ch, aromatic=True))
            i += 1
            continue

        # --- bond symbols ---
        if ch in BOND_CHARS:
            tokens.append(Token(type=TokenType.BOND, value=ch))
            i += 1
            continue

        # --- branch open / close ---
        if ch == "(":
            tokens.append(Token(type=TokenType.BRANCH_OPEN, value="("))
            i += 1
            continue

        if ch == ")":
            tokens.append(Token(type=TokenType.BRANCH_CLOSE, value=")"))
            i += 1
            continue

        # --- ring closure digit ---
        if ch.isdigit():
            tokens.append(Token(type=TokenType.RING_DIGIT, value=ch))
            i += 1
            continue

        # --- two-digit ring closure: %nn ---
        if ch == "%":
            if i + 2 >= n or not smiles[i + 1].isdigit() or not smiles[i + 2].isdigit():
                raise ValueError(
                    f"Expected two digits after '%' at position {i}")
            tokens.append(Token(
                type=TokenType.RING_DIGIT,
                value=smiles[i + 1:i + 3]))
            i += 3
            continue

        # --- dot disconnection ---
        if ch == ".":
            tokens.append(Token(type=TokenType.DOT, value="."))
            i += 1
            continue

        # --- unexpected character ---
        raise ValueError(
            f"Unexpected character {ch!r} at position {i} in SMILES "
            f"string {smiles!r}")

    return tokens
