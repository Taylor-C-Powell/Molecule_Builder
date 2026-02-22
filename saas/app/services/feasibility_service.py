"""Service layer for feasibility scoring.

Orchestrates the feasibility engine with timeout protection and
serialization to API response models.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, TimeoutError

from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis
from molbuilder.reactions.synthesis_route import extract_best_route
from molbuilder.reactions.feasibility import FeasibilityEngine, FeasibilityResult

from app.models.feasibility import (
    FeasibilityResponse,
    FeasibilityCompareResponse,
    DimensionScores,
    RouteSummary,
)

_executor = ThreadPoolExecutor(max_workers=2)


def _serialize_result(result: FeasibilityResult) -> FeasibilityResponse:
    """Convert a FeasibilityResult to API response model."""
    d = result.to_dict()
    return FeasibilityResponse(
        target_smiles=result.target_smiles,
        composite_score=result.composite_score,
        grade=result.grade,
        is_feasible=result.is_feasible,
        dimensions=DimensionScores(**d["dimensions"]),
        route_summary=RouteSummary(**d["route_summary"]),
    )


def _run_feasibility_sync(
    smiles: str, max_depth: int, beam_width: int
) -> FeasibilityResponse:
    """Run feasibility scoring synchronously (for thread pool)."""
    mol = parse(smiles)
    tree = retrosynthesis(mol, max_depth, beam_width)
    route = extract_best_route(tree)

    engine = FeasibilityEngine(max_depth=max_depth, beam_width=beam_width)
    result = engine.score(smiles, tree=tree, route=route)

    return _serialize_result(result)


def run_feasibility(
    smiles: str,
    max_depth: int = 5,
    beam_width: int = 5,
    timeout: int = 30,
) -> FeasibilityResponse:
    """Run feasibility scoring with timeout protection."""
    future = _executor.submit(_run_feasibility_sync, smiles, max_depth, beam_width)
    try:
        return future.result(timeout=timeout)
    except TimeoutError:
        future.cancel()
        raise ValueError(f"Feasibility scoring timed out after {timeout}s")


def _run_compare_sync(
    smiles_list: list[str], max_depth: int, beam_width: int
) -> FeasibilityCompareResponse:
    """Run comparison synchronously."""
    engine = FeasibilityEngine(max_depth=max_depth, beam_width=beam_width)
    results = engine.compare(smiles_list)
    responses = [_serialize_result(r) for r in results]

    best = responses[0] if responses else None
    return FeasibilityCompareResponse(
        results=responses,
        best_target=best.target_smiles if best else "",
        best_score=best.composite_score if best else 0.0,
    )


def run_feasibility_compare(
    smiles_list: list[str],
    max_depth: int = 5,
    beam_width: int = 5,
    timeout: int = 60,
) -> FeasibilityCompareResponse:
    """Compare feasibility of multiple targets with timeout."""
    future = _executor.submit(_run_compare_sync, smiles_list, max_depth, beam_width)
    try:
        return future.result(timeout=timeout)
    except TimeoutError:
        future.cancel()
        raise ValueError(f"Feasibility comparison timed out after {timeout}s")
