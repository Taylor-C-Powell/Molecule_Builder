"""RetroCast adapter -- convert MolBuilder RetrosynthesisTree to RetroCast Route schema.

The RetroCast canonical Route schema is a dict-based JSON-serializable format
used for benchmarking synthesis planners.  This module walks MolBuilder's
retrosynthesis tree and emits one Route dict per root-level disconnection,
with steps ordered from leaves (starting materials) to root (target) --
i.e. forward synthesis order.

Public API
----------
tree_to_retrocast_routes(tree) -> list[dict]
    Convert a RetrosynthesisTree to a list of RetroCast Route dicts.

export_retrocast_json(tree, path) -> None
    Write the routes as gzipped JSON to *path*.
"""

from __future__ import annotations

import gzip
import json
from typing import Any

import molbuilder
from molbuilder.reactions.retrosynthesis import (
    RetrosynthesisTree,
    RetroNode,
    Disconnection,
)


# =====================================================================
#  Internal helpers
# =====================================================================

def _make_molecule_dict(smiles: str, synthesis_step: int | None) -> dict[str, Any]:
    """Build a RetroCast Molecule dict."""
    return {
        "smiles": smiles,
        "synthesis_step": synthesis_step,
    }


def _collect_steps_dfs(
    node: RetroNode,
    steps: list[dict[str, Any]],
    mol_step_map: dict[str, int | None],
) -> None:
    """Walk the tree depth-first along best_disconnection, collecting steps.

    Recurses into children first so that leaf steps appear before their
    parents -- producing forward synthesis order (leaves first, target last).

    Parameters
    ----------
    node : RetroNode
        Current node being processed.
    steps : list[dict]
        Accumulator for ReactionStep dicts.  Modified in place.
    mol_step_map : dict[str, int | None]
        Maps SMILES -> step index that produces this molecule, or None if
        the molecule is a purchasable leaf.  Modified in place.
    """
    if node.is_purchasable or node.best_disconnection is None:
        # Leaf node -- not produced by any step
        mol_step_map.setdefault(node.smiles, None)
        return

    # Recurse into children first (depth-first, forward order)
    for child in node.children:
        _collect_steps_dfs(child, steps, mol_step_map)

    disc = node.best_disconnection
    step_index = len(steps)

    # Build reactant Molecule dicts
    reactants = []
    for precursor in disc.precursors:
        syn_step = mol_step_map.get(precursor.smiles)
        reactants.append(_make_molecule_dict(precursor.smiles, syn_step))

    # Record the product
    mol_step_map[node.smiles] = step_index

    product = _make_molecule_dict(node.smiles, step_index)

    step = {
        "index": step_index,
        "reaction_name": disc.template.name,
        "reactants": reactants,
        "product": product,
    }
    steps.append(step)


def _extract_route_from_root_disconnection(
    root: RetroNode,
    disconnection: Disconnection,
    original_children: list[RetroNode],
) -> dict[str, Any]:
    """Build one RetroCast Route dict from a specific root disconnection.

    Temporarily patches the root node so that its ``best_disconnection``
    points at *disconnection* and its ``children`` are matched accordingly,
    then walks the tree with ``_collect_steps_dfs``.

    Parameters
    ----------
    root : RetroNode
        Root (target) node of the retrosynthesis tree.
    disconnection : Disconnection
        The specific root-level disconnection to follow.
    original_children : list[RetroNode]
        The children that correspond to *disconnection*.  For the
        tree's actual best_disconnection these are ``root.children``;
        for alternate disconnections, we build leaf-only children.

    Returns
    -------
    dict
        A partial Route dict (without ``rank``; caller sets that).
    """
    steps: list[dict[str, Any]] = []
    mol_step_map: dict[str, int | None] = {}

    # Temporarily swap the root's best_disconnection and children
    saved_best = root.best_disconnection
    saved_children = root.children

    root.best_disconnection = disconnection
    root.children = original_children

    try:
        _collect_steps_dfs(root, steps, mol_step_map)
    finally:
        root.best_disconnection = saved_best
        root.children = saved_children

    target_mol = _make_molecule_dict(
        root.smiles,
        mol_step_map.get(root.smiles),
    )

    route: dict[str, Any] = {
        "target": target_mol,
        "steps": steps,
        "metadata": {
            "source": "MolBuilder",
            "version": molbuilder.__version__,
            "score": disconnection.score,
        },
    }
    return route


def _build_leaf_children(disconnection: Disconnection, depth: int) -> list[RetroNode]:
    """Build stub RetroNode children for an alternate disconnection.

    Since the tree only expands the *best* disconnection, alternate
    disconnections at the root level do not have pre-built child nodes.
    We create purchasable-leaf stubs for each precursor so that
    ``_collect_steps_dfs`` can still produce a valid (shallow) route.
    """
    from molbuilder.molecule.graph import Molecule

    children: list[RetroNode] = []
    for precursor in disconnection.precursors:
        mol = precursor.molecule
        if mol is None:
            # Create a minimal placeholder molecule
            mol = Molecule(name=precursor.name or precursor.smiles)
        child = RetroNode(
            smiles=precursor.smiles,
            molecule=mol,
            is_purchasable=True,
            depth=depth + 1,
        )
        children.append(child)
    return children


# =====================================================================
#  Public API
# =====================================================================

def tree_to_retrocast_routes(tree: RetrosynthesisTree) -> list[dict]:
    """Convert a MolBuilder RetrosynthesisTree to RetroCast Route dicts.

    Each disconnection at the root level produces one route.  The
    primary (best) disconnection follows the full tree expansion; alternate
    disconnections produce shallow single-step routes since the tree only
    fully expands along the best path.

    Steps are ordered from leaves (starting materials) to root (target) --
    i.e. forward synthesis order.

    Parameters
    ----------
    tree : RetrosynthesisTree
        The retrosynthesis tree produced by ``retrosynthesis()``.

    Returns
    -------
    list[dict]
        A list of Route dicts sorted by score (highest first).
        Each dict follows the RetroCast canonical schema::

            {
                "target": {"smiles": str, "synthesis_step": int | None},
                "steps": [{"index": int, "reaction_name": str,
                           "reactants": [...], "product": {...}}, ...],
                "rank": int,
                "metadata": {"source": "MolBuilder", "version": str,
                             "score": float},
            }
    """
    root = tree.target

    # If the target is purchasable (leaf), there are no synthesis steps
    if root.is_purchasable or not root.disconnections:
        return []

    routes: list[dict[str, Any]] = []

    for disc in root.disconnections:
        # The best disconnection uses the tree's fully expanded children
        if disc is root.best_disconnection:
            children = root.children
        else:
            # Alternate disconnections: build leaf-only stubs
            children = _build_leaf_children(disc, root.depth)

        route = _extract_route_from_root_disconnection(
            root, disc, children,
        )
        routes.append(route)

    # Sort by score descending (highest score = rank 1)
    routes.sort(key=lambda r: r["metadata"]["score"], reverse=True)

    # Assign ranks (1-indexed)
    for i, route in enumerate(routes):
        route["rank"] = i + 1

    return routes


def export_retrocast_json(tree: RetrosynthesisTree, path: str) -> None:
    """Export routes as gzipped JSON in the RetroCast format.

    Parameters
    ----------
    tree : RetrosynthesisTree
        The retrosynthesis tree to export.
    path : str
        Filesystem path for the output file (typically ending in ``.json.gz``).
    """
    routes = tree_to_retrocast_routes(tree)
    data = json.dumps(routes, indent=2, ensure_ascii=True)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(data)
