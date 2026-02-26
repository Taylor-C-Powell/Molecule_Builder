"""MolBuilder Python SDK -- typed client for the MolBuilder REST API.

Usage::

    from molbuilder_client import MolBuilder

    client = MolBuilder(api_key="mb_...")
    mol = client.from_smiles("CCO")
    print(mol.id, mol.num_atoms)
"""

from molbuilder_client._async_client import AsyncMolBuilder
from molbuilder_client._client import MolBuilder
from molbuilder_client._exceptions import (
    AuthenticationError,
    ForbiddenError,
    MolBuilderError,
    NotFoundError,
    RateLimitError,
    ServerError,
    ServiceUnavailableError,
    ValidationError,
)
from molbuilder_client._models import (
    APIKeyInfo,
    Atom3D,
    BatchJobSummary,
    BatchList,
    BatchStatus,
    BatchSubmit,
    BestRoute,
    BillingStatus,
    Bond3D,
    CheckoutSession,
    Conditions,
    CostBreakdown,
    CostEstimate,
    Disconnection,
    Element,
    Hazard,
    LibraryImport,
    LibraryList,
    LibraryMolecule,
    Molecule3D,
    MoleculeInfo,
    MoleculeProperties,
    Precursor,
    ProcessEvaluation,
    PurificationMethod,
    Reactor,
    RetroNode,
    RetrosynthesisPlan,
    RouteStep,
    SafetyAssessment,
    ScaleUp,
    StepDetail,
    Token,
)

__all__ = [
    # Clients
    "MolBuilder",
    "AsyncMolBuilder",
    # Exceptions
    "MolBuilderError",
    "AuthenticationError",
    "ForbiddenError",
    "NotFoundError",
    "ValidationError",
    "RateLimitError",
    "ServerError",
    "ServiceUnavailableError",
    # Models
    "APIKeyInfo",
    "Token",
    "MoleculeInfo",
    "MoleculeProperties",
    "Atom3D",
    "Bond3D",
    "Molecule3D",
    "Element",
    "Precursor",
    "Disconnection",
    "RetroNode",
    "RouteStep",
    "BestRoute",
    "RetrosynthesisPlan",
    "Reactor",
    "Conditions",
    "PurificationMethod",
    "StepDetail",
    "Hazard",
    "SafetyAssessment",
    "CostBreakdown",
    "CostEstimate",
    "ScaleUp",
    "ProcessEvaluation",
    "CheckoutSession",
    "BillingStatus",
    "BatchJobSummary",
    "BatchList",
    "BatchStatus",
    "BatchSubmit",
    "LibraryImport",
    "LibraryList",
    "LibraryMolecule",
]

__version__ = "0.1.0"
