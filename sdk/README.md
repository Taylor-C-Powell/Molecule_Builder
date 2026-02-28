# molbuilder-client

Python SDK for the [MolBuilder](https://www.molbuilder.io) REST API.

## Install

```bash
pip install molbuilder-client
```

Requires Python 3.11+. The only runtime dependency is [httpx](https://www.python-httpx.org/).

## Quickstart

```python
from molbuilder_client import MolBuilder

client = MolBuilder(api_key="mb_...")

# Parse a SMILES string
mol = client.from_smiles("CCO", name="ethanol")
print(mol.id, mol.num_atoms)  # mol_xxxx 9

# Get molecular properties
props = client.get_properties(mol.id)
print(props.formula, props.molecular_weight)  # C2H6O 46.07

# Get 3D coordinates
structure = client.get_3d(mol.id)
for atom in structure.atoms:
    print(atom.symbol, atom.position)

# Look up an element
fe = client.get_element("Fe")
print(fe.name, fe.atomic_weight)  # Iron 55.845

# Retrosynthesis
plan = client.retrosynthesis("c1ccccc1O")  # phenol
print(plan.routes_found, plan.best_route.total_steps)

# Process evaluation
proc = client.process_evaluate("CCO", scale_kg=100.0)
print(proc.cost.total_usd, proc.cost.per_kg_usd)
```

## Async usage

```python
import asyncio
from molbuilder_client import AsyncMolBuilder

async def main():
    async with AsyncMolBuilder(api_key="mb_...") as client:
        mol = await client.from_smiles("CCO")
        props = await client.get_properties(mol.id)
        print(props.formula)

asyncio.run(main())
```

## Error handling

All API errors raise typed exceptions that inherit from `MolBuilderError`:

```python
from molbuilder_client import MolBuilder, ValidationError, RateLimitError

client = MolBuilder(api_key="mb_...")

try:
    client.from_smiles("not_valid")
except ValidationError as e:
    print(f"Bad SMILES: {e.message}")
except RateLimitError as e:
    print(f"Slow down! Retry after {e.retry_after}s")
```

| HTTP Status | Exception                  |
|-------------|----------------------------|
| 401         | `AuthenticationError`      |
| 403         | `ForbiddenError`           |
| 404         | `NotFoundError`            |
| 422         | `ValidationError`          |
| 429         | `RateLimitError`           |
| 500         | `ServerError`              |
| 501/503     | `ServiceUnavailableError`  |

## API reference

### Auth

| Method                    | Returns        | Description                         |
|---------------------------|----------------|-------------------------------------|
| `register(email)`         | `APIKeyInfo`   | Register a new free-tier API key    |
| `get_token()`             | `Token`        | Exchange API key for a JWT          |

### Molecule

| Method                        | Returns              | Description                  |
|-------------------------------|----------------------|------------------------------|
| `from_smiles(smiles, name="")` | `MoleculeInfo`     | Parse SMILES into a molecule |
| `get_properties(mol_id)`      | `MoleculeProperties` | Computed molecular properties |
| `get_3d(mol_id)`              | `Molecule3D`        | 3D coordinates and bonds     |

### Elements

| Method                | Returns   | Description                  |
|-----------------------|-----------|------------------------------|
| `get_element(symbol)` | `Element` | Look up element by symbol    |

### Retrosynthesis

| Method                                              | Returns              | Description              |
|-----------------------------------------------------|----------------------|--------------------------|
| `retrosynthesis(smiles, max_depth=5, beam_width=5)` | `RetrosynthesisPlan` | Plan retrosynthetic routes |

### Process

| Method                                                            | Returns             | Description                    |
|-------------------------------------------------------------------|---------------------|--------------------------------|
| `process_evaluate(smiles, scale_kg=1.0, max_depth=5, beam_width=5)` | `ProcessEvaluation` | Full process engineering report |

### Batch Processing

| Method                                          | Returns       | Description                       |
|-------------------------------------------------|---------------|-----------------------------------|
| `submit_batch(smiles_list, job_type, **params)` | `BatchJob`    | Submit async batch job            |
| `get_batch_status(job_id)`                      | `BatchJob`    | Poll job status and results       |
| `list_batches(limit=20, offset=0)`              | `BatchList`   | List your batch jobs              |
| `cancel_batch(job_id)`                          | `BatchJob`    | Cancel a running batch job        |
| `wait_for_batch(job_id, timeout=300)`           | `BatchJob`    | Poll until complete or timeout    |

### Molecule Library

| Method                                   | Returns          | Description                   |
|------------------------------------------|------------------|-------------------------------|
| `save_molecule(mol_id, **kwargs)`        | `LibraryEntry`   | Save molecule to library      |
| `get_molecule(entry_id)`                 | `LibraryEntry`   | Get library entry by ID       |
| `list_molecules(limit=20, offset=0)`     | `LibraryList`    | List saved molecules          |
| `update_molecule(entry_id, **kwargs)`    | `LibraryEntry`   | Update library entry metadata |
| `delete_molecule(entry_id)`              | `None`           | Delete library entry          |

### File I/O

| Method                                          | Returns           | Description                    |
|-------------------------------------------------|-------------------|--------------------------------|
| `import_file(file_path)`                        | `FileImportResult`| Upload and parse a molecule file |
| `export_file(mol_id, format, save_to=None)`     | `bytes`           | Export molecule in given format |

### Reports

| Method                                             | Returns | Description                      |
|----------------------------------------------------|---------|----------------------------------|
| `download_report(smiles, scale_kg=1.0, save_to=None)` | `bytes` | Generate and download PDF report |

### Billing

| Method                  | Returns           | Description                        |
|-------------------------|-------------------|------------------------------------|
| `create_checkout(plan)` | `CheckoutSession` | Create Stripe checkout session     |
| `billing_status()`      | `BillingStatus`   | Current subscription status        |

## License

Apache 2.0
