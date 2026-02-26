"""SMARTS pattern matching engine for substructure search.

Provides parsing of SMARTS strings into pattern graphs and VF2-style
subgraph isomorphism matching against Molecule objects.

Public API
----------
parse_smarts(smarts) -> SmartsPattern
    Parse a SMARTS string into a pattern graph.
find_matches(mol, pattern) -> list[dict[int, int]]
    Find all subgraph matches of *pattern* in *mol*.
has_match(mol, pattern) -> bool
    Return True if *pattern* matches anywhere in *mol*.
detect_by_smarts(mol, smarts_str, group_name) -> list[FunctionalGroup]
    Convenience wrapper that returns FunctionalGroup objects.
"""

from molbuilder.smarts.parser import parse_smarts
from molbuilder.smarts.matcher import find_matches, has_match, detect_by_smarts
