"""File I/O endpoints: import and export molecule files."""

from fastapi import APIRouter, Depends, UploadFile
from fastapi.responses import Response

from app.dependencies import UserContext, check_rate_limit
from app.exceptions import MoleculeNotFound, InvalidFileFormat
from app.models.file_io import FileImportResponse
from app.models.molecule import MoleculeResponse
from app.services.file_io_service import (
    detect_format,
    parse_file,
    export_molecule,
    CONTENT_TYPES,
)
from app.services.molecule_store import molecule_store
from molbuilder.smiles.writer import to_smiles

router = APIRouter(prefix="/api/v1/molecule", tags=["file-io"])


@router.post("/import-file", response_model=FileImportResponse)
async def import_file(
    file: UploadFile,
    user: UserContext = Depends(check_rate_limit),
):
    """Import a molecule file (XYZ, MOL, SDF, PDB, or JSON)."""
    raw = await file.read()
    content = raw.decode("utf-8", errors="replace")
    filename = file.filename or "unknown"

    fmt = detect_format(filename, content)
    molecules = parse_file(content, fmt)

    results = []
    for mol in molecules:
        try:
            smiles = to_smiles(mol)
        except Exception:
            smiles = mol.name or "unknown"
        mol_id = molecule_store.put(mol, smiles)
        results.append(MoleculeResponse(
            id=mol_id,
            name=mol.name,
            smiles=smiles,
            num_atoms=len(mol.atoms),
            num_bonds=len(mol.bonds),
        ))

    return FileImportResponse(
        molecules=results,
        format=fmt,
        count=len(results),
    )


@router.get("/{mol_id}/export/{fmt}")
def export_file(
    mol_id: str,
    fmt: str,
    user: UserContext = Depends(check_rate_limit),
):
    """Export a molecule in the specified format (xyz, mol, pdb, json)."""
    item = molecule_store.get(mol_id)
    if item is None:
        raise MoleculeNotFound(mol_id)
    mol, _ = item

    if fmt not in CONTENT_TYPES:
        raise InvalidFileFormat(fmt)

    content = export_molecule(mol, fmt)
    media_type = CONTENT_TYPES[fmt]
    ext = fmt
    filename = f"{mol.name or mol_id}.{ext}"

    return Response(
        content=content,
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )
