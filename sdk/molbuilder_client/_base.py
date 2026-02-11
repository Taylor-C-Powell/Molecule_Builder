"""Shared utilities for sync and async clients.

Centralises URL construction, header building, error mapping, and
JSON-to-dataclass deserialization so both client variants stay DRY.
"""

from __future__ import annotations

import dataclasses
import types
import typing
from typing import Any, TypeVar, get_type_hints

import httpx

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

T = TypeVar("T")

DEFAULT_BASE_URL = "https://molbuilder-api-production.up.railway.app"
API_PREFIX = "/api/v1"
USER_AGENT = "molbuilder-client/0.1.0"


# ---------------------------------------------------------------------------
# Header / URL helpers
# ---------------------------------------------------------------------------

def build_headers(api_key: str) -> dict[str, str]:
    return {
        "X-API-Key": api_key,
        "User-Agent": USER_AGENT,
        "Accept": "application/json",
    }


def build_url(base_url: str, path: str) -> str:
    """Join base URL and path, avoiding double slashes."""
    return f"{base_url.rstrip('/')}{API_PREFIX}/{path.lstrip('/')}"


# ---------------------------------------------------------------------------
# Error mapping
# ---------------------------------------------------------------------------

def raise_for_status(response: httpx.Response) -> None:
    """Map HTTP error status codes to typed SDK exceptions."""
    if response.is_success:
        return

    code = response.status_code
    try:
        body = response.json()
        message = body.get("detail") or body.get("error") or response.text
    except Exception:
        message = response.text or f"HTTP {code}"

    if code == 401:
        raise AuthenticationError(str(message))
    if code == 403:
        raise ForbiddenError(str(message))
    if code == 404:
        raise NotFoundError(str(message))
    if code == 422:
        raise ValidationError(str(message))
    if code == 429:
        retry = response.headers.get("Retry-After")
        retry_after = int(retry) if retry else None
        raise RateLimitError(str(message), retry_after=retry_after)
    if code in (501, 503):
        raise ServiceUnavailableError(str(message))
    if code >= 500:
        raise ServerError(str(message))

    raise MolBuilderError(str(message), status_code=code)


# ---------------------------------------------------------------------------
# JSON → dataclass deserializer
# ---------------------------------------------------------------------------

def _get_origin(tp: Any) -> Any:
    return getattr(tp, "__origin__", None)


def _get_args(tp: Any) -> tuple[Any, ...]:
    return getattr(tp, "__args__", ())


def _is_dataclass_type(tp: Any) -> bool:
    return isinstance(tp, type) and dataclasses.is_dataclass(tp)


def _resolve_optional(tp: Any) -> tuple[Any, bool]:
    """If *tp* is ``X | None``, return ``(X, True)``.  Else ``(tp, False)``."""
    # types.UnionType (X | Y syntax) has no __origin__; check isinstance first
    if isinstance(tp, types.UnionType):
        args = [a for a in tp.__args__ if a is not type(None)]
        if len(args) == 1:
            return args[0], True
    origin = _get_origin(tp)
    if origin is typing.Union:
        args = [a for a in _get_args(tp) if a is not type(None)]
        if len(args) == 1:
            return args[0], True
    return tp, False


def _coerce_value(value: Any, target_type: Any) -> Any:
    """Recursively coerce *value* into *target_type*."""
    if value is None:
        return None

    inner, is_optional = _resolve_optional(target_type)
    if is_optional:
        return _coerce_value(value, inner)

    # Dataclass → recursive descent
    if _is_dataclass_type(inner):
        return from_dict(inner, value)

    origin = _get_origin(inner)

    # list[X]
    if origin is list:
        (item_type,) = _get_args(inner)
        return [_coerce_value(item, item_type) for item in value]

    # tuple[X, Y, Z] (fixed-length, e.g. position)
    if origin is tuple:
        args = _get_args(inner)
        if args and args[-1] is not Ellipsis:
            return tuple(_coerce_value(v, t) for v, t in zip(value, args))
        return tuple(value)

    return value


def from_dict(cls: type[T], data: dict[str, Any]) -> T:
    """Construct a frozen dataclass *cls* from a JSON-derived *data* dict.

    Handles nested dataclasses, ``list[X]``, ``X | None``, and
    ``tuple[float, ...]`` fields.  Unknown keys in *data* are silently
    ignored so the SDK stays forward-compatible with new server fields.
    """
    hints = get_type_hints(cls)
    fields = {f.name for f in dataclasses.fields(cls)}
    kwargs: dict[str, Any] = {}

    for name in fields:
        if name not in data:
            continue
        kwargs[name] = _coerce_value(data[name], hints[name])

    return cls(**kwargs)
