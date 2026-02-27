"""File import/export service for molecule files."""

from __future__ import annotations

from molbuilder.molecule.graph import Molecule
from molbuilder.io.xyz import from_xyz_string, to_xyz_string
from molbuilder.io.mol_sdf import from_mol_string, to_mol_string
from molbuilder.io.pdb import from_pdb_string, to_pdb_string
from molbuilder.io.json_io import from_json_string, to_json_string

from app.exceptions import InvalidFileFormat


# Supported formats and their content types
CONTENT_TYPES: dict[str, str] = {
    "xyz": "chemical/x-xyz",
    "mol": "chemical/x-mdl-molfile",
    "sdf": "chemical/x-mdl-sdfile",
    "pdb": "chemical/x-pdb",
    "json": "application/json",
}

# Maximum molecules per SDF upload
MAX_SDF_MOLECULES = 100


def detect_format(filename: str, content: str) -> str:
    """Detect file format from extension or content sniffing.

    Returns one of: 'xyz', 'mol', 'sdf', 'pdb', 'json'.
    Raises InvalidFileFormat if unrecognised.
    """
    ext = filename.rsplit(".", 1)[-1].lower() if "." in filename else ""

    ext_map = {
        "xyz": "xyz",
        "mol": "mol",
        "sdf": "sdf",
        "pdb": "pdb",
        "json": "json",
    }
    if ext in ext_map:
        return ext_map[ext]

    # Content sniffing as fallback
    stripped = content.strip()
    if stripped.startswith("{"):
        return "json"
    if "$$$$" in content:
        return "sdf"
    lines = stripped.splitlines()
    if len(lines) >= 2 and lines[0].strip().isdigit():
        return "xyz"
    if any(ln.startswith(("ATOM", "HETATM")) for ln in lines[:20]):
        return "pdb"
    # Check for V2000 counts line pattern
    if len(lines) >= 4:
        counts_line = lines[3].strip()
        if "V2000" in counts_line or (
            len(counts_line.split()) >= 2
            and counts_line.split()[0].isdigit()
            and counts_line.split()[1].isdigit()
        ):
            return "mol"

    raise InvalidFileFormat(filename)


def parse_file(content: str, fmt: str) -> list[Molecule]:
    """Parse file content into Molecule objects.

    For SDF, splits on $$$$ and parses each MOL block.
    Returns a list (single-element for non-SDF formats).
    """
    try:
        if fmt == "xyz":
            return [from_xyz_string(content)]
        elif fmt == "mol":
            return [from_mol_string(content)]
        elif fmt == "sdf":
            blocks = content.split("$$$$")
            molecules = []
            for block in blocks:
                block = block.strip()
                if not block:
                    continue
                if len(molecules) >= MAX_SDF_MOLECULES:
                    break
                molecules.append(from_mol_string(block))
            return molecules
        elif fmt == "pdb":
            return [from_pdb_string(content)]
        elif fmt == "json":
            return [from_json_string(content)]
        else:
            raise InvalidFileFormat(fmt)
    except InvalidFileFormat:
        raise
    except Exception as e:
        raise InvalidFileFormat(f"{fmt}: {e}")


def export_molecule(mol: Molecule, fmt: str) -> str:
    """Export a Molecule to the specified format string."""
    if fmt == "xyz":
        return to_xyz_string(mol)
    elif fmt == "mol":
        return to_mol_string(mol)
    elif fmt == "pdb":
        return to_pdb_string(mol)
    elif fmt == "json":
        return to_json_string(mol)
    else:
        raise InvalidFileFormat(fmt)
