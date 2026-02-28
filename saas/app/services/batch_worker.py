"""Batch processing worker using ThreadPoolExecutor."""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, Future, TimeoutError as FuturesTimeoutError

from app.services.job_db import get_job_db

logger = logging.getLogger("molbuilder.batch")

# Maximum time (seconds) a single batch job is allowed to run before timeout
BATCH_JOB_TIMEOUT_SECONDS = 300


def _process_molecule(smiles: str, job_type: str, params: dict) -> dict:
    """Process a single molecule. Returns result dict or error dict."""
    try:
        from molbuilder.smiles.parser import parse
        from molbuilder.smiles.writer import to_smiles

        mol = parse(smiles)
        canonical = to_smiles(mol)

        if job_type == "properties":
            from molbuilder.molecule.properties import lipinski_properties
            from molbuilder.reactions.functional_group_detect import detect_functional_groups
            from molbuilder.io.json_io import _molecular_formula
            from molbuilder.core.elements import atomic_weight
            fgs = detect_functional_groups(mol)
            props = lipinski_properties(mol)
            mw = sum(atomic_weight(a.symbol) for a in mol.atoms)
            return {
                "smiles": canonical,
                "molecular_weight": round(mw, 4),
                "formula": _molecular_formula(mol),
                "logp": round(props.logp, 4),
                "heavy_atom_count": props.heavy_atom_count,
                "functional_groups": [fg.name for fg in fgs],
                "status": "success",
            }
        elif job_type == "retrosynthesis":
            from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree
            tree = retrosynthesis(mol, max_depth=params.get("max_depth", 3))
            return {
                "smiles": canonical,
                "tree": format_tree(tree),
                "status": "success",
            }
        elif job_type == "conditions":
            from molbuilder.process.condition_prediction import predict_conditions
            pred = predict_conditions(
                smiles,
                reaction_name=params.get("reaction_name"),
                scale_kg=params.get("scale_kg", 1.0),
            )
            best = None
            if pred.best_match:
                best = {
                    "template": pred.best_match.template_name,
                    "score": pred.best_match.match_score,
                    "conditions": {
                        "temperature_c": pred.best_match.conditions.temperature_celsius,
                        "solvent": pred.best_match.conditions.solvent,
                    },
                }
            return {
                "smiles": canonical,
                "confidence": pred.overall_confidence,
                "best_match": best,
                "num_candidates": len(pred.candidates),
                "status": "success",
            }
        elif job_type == "evaluate":
            return {
                "smiles": canonical,
                "name": mol.name,
                "num_atoms": len(mol.atoms),
                "num_bonds": len(mol.bonds),
                "status": "success",
            }
        else:
            return {"smiles": smiles, "status": "error", "error": f"Unknown job type: {job_type}"}
    except Exception as e:
        return {"smiles": smiles, "status": "error", "error": str(e)}


def _run_batch(job_id: str, smiles_list: list[str], job_type: str, params: dict) -> None:
    """Execute a batch job, updating progress in the DB."""
    db = get_job_db()
    db.update_status(job_id, "running", 0.0)

    results = []
    total = len(smiles_list)
    errors = 0

    for i, smiles in enumerate(smiles_list):
        # Check if job was cancelled or timed out
        job = db.get_job(job_id)
        if job and job["status"] in ("cancelled", "failed"):
            return

        result = _process_molecule(smiles, job_type, params)
        results.append(result)
        if result.get("status") == "error":
            errors += 1

        progress = ((i + 1) / total) * 100.0
        db.update_status(job_id, "running", progress)

    # Final status check: do not overwrite if job was cancelled or timed out
    job = db.get_job(job_id)
    if job and job["status"] in ("cancelled", "failed"):
        return

    db.set_result(job_id, {
        "results": results,
        "total": total,
        "succeeded": total - errors,
        "failed": errors,
    })


class BatchWorker:
    """Manages async batch processing with a thread pool."""

    def __init__(self, max_workers: int = 4):
        self._executor = ThreadPoolExecutor(max_workers=max_workers)
        self._futures: dict[str, Future] = {}

    def submit_batch(self, job_id: str, smiles_list: list[str],
                     job_type: str, params: dict | None = None) -> None:
        future = self._executor.submit(
            self._run_with_timeout, job_id, smiles_list, job_type, params or {}
        )
        self._futures[job_id] = future

    @staticmethod
    def _run_with_timeout(job_id: str, smiles_list: list[str],
                          job_type: str, params: dict) -> None:
        """Run a batch job with a timeout guard."""
        import threading

        result_holder: list[Exception | None] = [None]
        done_event = threading.Event()

        def _target():
            try:
                _run_batch(job_id, smiles_list, job_type, params)
            except Exception as exc:
                result_holder[0] = exc
            finally:
                done_event.set()

        thread = threading.Thread(target=_target, daemon=True)
        thread.start()

        if not done_event.wait(timeout=BATCH_JOB_TIMEOUT_SECONDS):
            logger.error("Batch job %s timed out after %d seconds",
                         job_id, BATCH_JOB_TIMEOUT_SECONDS)
            db = get_job_db()
            db.update_status(job_id, "failed", 0.0)
            db.set_error(job_id, f"Job timed out after {BATCH_JOB_TIMEOUT_SECONDS} seconds")
            return

        if result_holder[0] is not None:
            raise result_holder[0]

    def cancel(self, job_id: str) -> bool:
        future = self._futures.get(job_id)
        if future and not future.done():
            return future.cancel()
        return False

    def shutdown(self) -> None:
        self._executor.shutdown(wait=False)


_batch_worker: BatchWorker | None = None


def get_batch_worker() -> BatchWorker:
    global _batch_worker
    if _batch_worker is None:
        _batch_worker = BatchWorker()
    return _batch_worker


def set_batch_worker(worker: BatchWorker | None) -> None:
    global _batch_worker
    _batch_worker = worker
