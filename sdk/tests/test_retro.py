"""Tests for the retrosynthesis endpoint."""

from __future__ import annotations

import httpx
import respx

from molbuilder_client import (
    BestRoute,
    Disconnection,
    MolBuilder,
    Precursor,
    RetroNode,
    RetrosynthesisPlan,
    RouteStep,
)


RETRO_RESPONSE = {
    "tree": {
        "smiles": "CCO",
        "is_purchasable": False,
        "functional_groups": ["alcohol"],
        "best_disconnection": {
            "reaction_name": "hydration",
            "named_reaction": "Markovnikov hydration",
            "category": "addition",
            "score": 0.85,
            "precursors": [
                {"smiles": "C=C", "name": "ethylene", "cost_per_kg": 0.5},
            ],
        },
        "children": [
            {
                "smiles": "C=C",
                "is_purchasable": True,
                "functional_groups": ["alkene"],
                "best_disconnection": None,
                "children": [],
            },
        ],
    },
    "routes_found": 1,
    "max_depth": 5,
    "beam_width": 5,
    "best_route": {
        "target_smiles": "CCO",
        "target_name": "ethanol",
        "total_steps": 1,
        "overall_yield": 0.85,
        "longest_linear_sequence": 1,
        "starting_materials": [
            {"smiles": "C=C", "name": "ethylene", "cost_per_kg": 0.5},
        ],
        "steps": [
            {
                "step_number": 1,
                "reaction_name": "hydration",
                "named_reaction": "Markovnikov hydration",
                "category": "addition",
                "precursor_smiles": ["C=C"],
                "product_smiles": "CCO",
                "product_name": "ethanol",
                "conditions": "H2SO4, H2O, heat",
                "expected_yield": 0.85,
                "notes": "",
            },
        ],
    },
}


def test_retrosynthesis(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/retrosynthesis/plan").mock(
        return_value=httpx.Response(200, json=RETRO_RESPONSE)
    )

    result = client.retrosynthesis("CCO")

    assert isinstance(result, RetrosynthesisPlan)
    assert result.routes_found == 1
    assert result.max_depth == 5

    # Tree structure
    assert isinstance(result.tree, RetroNode)
    assert result.tree.smiles == "CCO"
    assert result.tree.is_purchasable is False
    assert result.tree.functional_groups == ["alcohol"]

    # Disconnection
    disc = result.tree.best_disconnection
    assert isinstance(disc, Disconnection)
    assert disc.reaction_name == "hydration"
    assert disc.score == 0.85
    assert len(disc.precursors) == 1
    assert isinstance(disc.precursors[0], Precursor)
    assert disc.precursors[0].smiles == "C=C"

    # Children (recursive)
    assert len(result.tree.children) == 1
    child = result.tree.children[0]
    assert isinstance(child, RetroNode)
    assert child.is_purchasable is True
    assert child.best_disconnection is None
    assert child.children == []

    # Best route
    route = result.best_route
    assert isinstance(route, BestRoute)
    assert route.total_steps == 1
    assert route.overall_yield == 0.85
    assert len(route.steps) == 1
    assert isinstance(route.steps[0], RouteStep)
    assert route.steps[0].conditions == "H2SO4, H2O, heat"


def test_retrosynthesis_custom_params(client: MolBuilder, mock_api: respx.Router) -> None:
    route = mock_api.post("/retrosynthesis/plan").mock(
        return_value=httpx.Response(200, json=RETRO_RESPONSE)
    )

    client.retrosynthesis("CCO", max_depth=3, beam_width=2)

    import json
    body = json.loads(route.calls.last.request.content)
    assert body["max_depth"] == 3
    assert body["beam_width"] == 2


def test_retrosynthesis_null_best_route(client: MolBuilder, mock_api: respx.Router) -> None:
    response = {**RETRO_RESPONSE, "best_route": None}
    mock_api.post("/retrosynthesis/plan").mock(
        return_value=httpx.Response(200, json=response)
    )

    result = client.retrosynthesis("CCO")
    assert result.best_route is None
