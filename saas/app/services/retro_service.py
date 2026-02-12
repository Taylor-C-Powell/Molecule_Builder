"""Retrosynthesis service wrapping molbuilder's retro engine."""

from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeoutError

from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis
from molbuilder.reactions.synthesis_route import extract_best_route
from app.models.retrosynthesis import (
    RetroNodeResponse,
    DisconnectionResponse,
    PrecursorResponse,
    SynthesisStepResponse,
    SynthesisRouteResponse,
    RetroResponse,
)


def _serialize_node(node) -> RetroNodeResponse:
    best = None
    if node.best_disconnection:
        d = node.best_disconnection
        best = DisconnectionResponse(
            reaction_name=d.template.name,
            named_reaction=d.template.named_reaction,
            category=d.template.category.name,
            score=d.score,
            precursors=[
                PrecursorResponse(
                    smiles=p.smiles, name=p.name, cost_per_kg=p.cost_per_kg
                )
                for p in d.precursors
            ],
        )

    return RetroNodeResponse(
        smiles=node.smiles,
        is_purchasable=node.is_purchasable,
        functional_groups=[fg.name for fg in node.functional_groups],
        best_disconnection=best,
        children=[_serialize_node(c) for c in node.children],
    )


def _serialize_route(route) -> SynthesisRouteResponse:
    return SynthesisRouteResponse(
        target_smiles=route.target_smiles,
        target_name=route.target_name,
        total_steps=route.total_steps,
        overall_yield=route.overall_yield,
        longest_linear_sequence=route.longest_linear_sequence,
        starting_materials=[
            PrecursorResponse(
                smiles=p.smiles, name=p.name, cost_per_kg=p.cost_per_kg
            )
            for p in route.starting_materials
        ],
        steps=[
            SynthesisStepResponse(
                step_number=s.step_number,
                reaction_name=s.template.name,
                named_reaction=s.template.named_reaction,
                category=s.template.category.name,
                precursor_smiles=[p.smiles for p in s.precursors],
                product_smiles=s.product_smiles,
                product_name=s.product_name,
                conditions=s.conditions,
                expected_yield=s.expected_yield,
                notes=s.notes,
            )
            for s in route.steps
        ],
    )


RETRO_TIMEOUT_SECONDS = 30

_executor = ThreadPoolExecutor(max_workers=2)


def _run_retro_sync(
    smiles: str, max_depth: int, beam_width: int
) -> RetroResponse:
    """Core retrosynthesis logic (runs in thread pool)."""
    mol = parse(smiles)
    tree = retrosynthesis(mol, max_depth=max_depth, beam_width=beam_width)

    best_route = None
    if tree.target.best_disconnection:
        route = extract_best_route(tree)
        best_route = _serialize_route(route)

    return RetroResponse(
        tree=_serialize_node(tree.target),
        routes_found=tree.routes_found,
        max_depth=tree.max_depth,
        beam_width=tree.beam_width,
        best_route=best_route,
    )


def run_retrosynthesis(
    smiles: str, max_depth: int = 5, beam_width: int = 5
) -> RetroResponse:
    """Run retrosynthesis with timeout protection."""
    future = _executor.submit(_run_retro_sync, smiles, max_depth, beam_width)
    try:
        return future.result(timeout=RETRO_TIMEOUT_SECONDS)
    except FuturesTimeoutError:
        future.cancel()
        raise ValueError(
            f"Retrosynthesis timed out after {RETRO_TIMEOUT_SECONDS}s "
            f"for SMILES: {smiles}"
        )
