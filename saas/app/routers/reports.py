"""Report endpoints: PDF generation for process engineering reports."""

from fastapi import APIRouter, Depends, Query
from fastapi.responses import Response

from app.dependencies import UserContext, check_expensive_rate_limit
from app.exceptions import InvalidSMILES

router = APIRouter(prefix="/api/v1/reports", tags=["reports"])


@router.post("/process-pdf")
def process_pdf(
    smiles: str = Query(..., min_length=1, max_length=2000, description="Target SMILES"),
    scale_kg: float = Query(1.0, gt=0, description="Production scale in kg"),
    user: UserContext = Depends(check_expensive_rate_limit),
):
    """Generate a PDF process engineering report for a molecule.

    Runs retrosynthesis, process evaluation, and safety assessment,
    then renders all results into a professional PDF document.
    """
    from molbuilder.smiles.parser import parse
    from molbuilder.reactions.retrosynthesis import retrosynthesis
    from molbuilder.reactions.synthesis_route import extract_best_route
    from molbuilder.process.costing import estimate_cost, CostParameters
    from molbuilder.process.safety import assess_safety
    from molbuilder.reports.pdf_report import generate_process_pdf

    try:
        mol = parse(smiles)
    except (ValueError, KeyError) as e:
        raise InvalidSMILES(smiles, str(e))

    try:
        tree = retrosynthesis(mol)
        route = extract_best_route(tree)
    except (ValueError, KeyError) as e:
        raise InvalidSMILES(smiles, str(e))

    if route is None:
        raise InvalidSMILES(smiles, "No synthesis route found")

    # Cost and safety
    cost = None
    safety = None
    try:
        cost = estimate_cost(route, CostParameters(scale_kg=scale_kg))
    except Exception:
        pass
    try:
        safety = assess_safety(route)
    except Exception:
        pass

    pdf_bytes = generate_process_pdf(
        mol_name=mol.name or smiles,
        route=route,
        cost=cost,
        safety=safety,
    )

    filename = f"process_report_{smiles[:20].replace('/', '_')}.pdf"
    return Response(
        content=pdf_bytes,
        media_type="application/pdf",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )
