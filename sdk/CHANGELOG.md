# Changelog

## v0.2.0 (2026-02-26)

### Added

- **File I/O**: `import_file()` and `export_file()` methods for uploading/downloading molecule files (XYZ, MOL, SDF, PDB, JSON)
- **PDF Reports**: `download_report()` method for generating process engineering PDF reports
- **SA Score**: `MoleculeProperties.sa_score` field for synthetic accessibility scoring (1-10 scale)
- **FileImportResult** dataclass for file import responses

### Changed

- Bumped `USER_AGENT` to `molbuilder-client/0.2.0`

## v0.1.1 (2026-02-26)

### Added

- **Library**: `library_save()`, `library_get()`, `library_list()`, `library_update()`, `library_delete()`, `library_import()` methods
- **Batch**: `batch_submit()`, `batch_status()`, `batch_list()`, `batch_cancel()` methods
- Corresponding async methods on `AsyncMolBuilder`
- `LibraryMolecule`, `LibraryList`, `LibraryImport`, `BatchSubmit`, `BatchStatus`, `BatchList`, `BatchJobSummary` dataclasses

## v0.1.0 (2026-02-25)

### Added

- Initial release
- `MolBuilder` sync client and `AsyncMolBuilder` async client
- Authentication: `register()`, `get_token()`
- Molecule: `from_smiles()`, `get_properties()`, `get_3d()`
- Elements: `get_element()`
- Retrosynthesis: `retrosynthesis()`
- Process evaluation: `process_evaluate()`
- Billing: `create_checkout()`, `billing_status()`
- Full typed dataclass models for all API responses
- Typed exception hierarchy for API errors
