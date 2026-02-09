"""In-memory usage tracking for API analytics."""

import threading
import time
from collections import defaultdict
from dataclasses import dataclass, field


@dataclass
class EndpointStats:
    """Aggregate stats for a single endpoint."""
    total_calls: int = 0
    total_errors: int = 0
    total_latency_ms: float = 0.0
    min_latency_ms: float = float("inf")
    max_latency_ms: float = 0.0

    @property
    def avg_latency_ms(self) -> float:
        if self.total_calls == 0:
            return 0.0
        return self.total_latency_ms / self.total_calls


@dataclass
class RequestRecord:
    """A single API request log entry."""
    timestamp: float
    method: str
    path: str
    status_code: int
    latency_ms: float
    user_email: str
    user_tier: str
    query_params: dict = field(default_factory=dict)
    body_summary: dict = field(default_factory=dict)


class UsageTracker:
    """Thread-safe in-memory usage tracker.

    Stores per-endpoint aggregate stats and a rolling window of
    recent individual request records for inspection.
    """

    def __init__(self, max_records: int = 10_000) -> None:
        self._lock = threading.Lock()
        self._max_records = max_records
        self._records: list[RequestRecord] = []
        self._endpoint_stats: dict[str, EndpointStats] = defaultdict(EndpointStats)
        self._tier_calls: dict[str, int] = defaultdict(int)
        self._user_calls: dict[str, int] = defaultdict(int)
        self._smiles_queries: dict[str, int] = defaultdict(int)
        self._start_time: float = time.time()

    def record(self, rec: RequestRecord) -> None:
        with self._lock:
            # Rolling window
            self._records.append(rec)
            if len(self._records) > self._max_records:
                self._records = self._records[-self._max_records:]

            # Endpoint aggregates
            key = f"{rec.method} {rec.path}"
            stats = self._endpoint_stats[key]
            stats.total_calls += 1
            if rec.status_code >= 400:
                stats.total_errors += 1
            stats.total_latency_ms += rec.latency_ms
            if rec.latency_ms < stats.min_latency_ms:
                stats.min_latency_ms = rec.latency_ms
            if rec.latency_ms > stats.max_latency_ms:
                stats.max_latency_ms = rec.latency_ms

            # Per-tier and per-user counts
            if rec.user_tier:
                self._tier_calls[rec.user_tier] += 1
            if rec.user_email:
                self._user_calls[rec.user_email] += 1

            # Track SMILES queries
            smiles = rec.body_summary.get("smiles")
            if smiles:
                self._smiles_queries[smiles] += 1

    def get_summary(self) -> dict:
        with self._lock:
            uptime_s = time.time() - self._start_time
            total = sum(s.total_calls for s in self._endpoint_stats.values())
            total_errors = sum(s.total_errors for s in self._endpoint_stats.values())

            endpoints = {}
            for path, stats in sorted(
                self._endpoint_stats.items(),
                key=lambda x: x[1].total_calls,
                reverse=True,
            ):
                endpoints[path] = {
                    "calls": stats.total_calls,
                    "errors": stats.total_errors,
                    "avg_latency_ms": round(stats.avg_latency_ms, 2),
                    "min_latency_ms": round(stats.min_latency_ms, 2)
                    if stats.min_latency_ms != float("inf")
                    else 0,
                    "max_latency_ms": round(stats.max_latency_ms, 2),
                }

            top_smiles = sorted(
                self._smiles_queries.items(), key=lambda x: x[1], reverse=True
            )[:20]

            top_users = sorted(
                self._user_calls.items(), key=lambda x: x[1], reverse=True
            )[:20]

            return {
                "uptime_seconds": round(uptime_s, 1),
                "total_requests": total,
                "total_errors": total_errors,
                "error_rate": round(total_errors / total, 4) if total > 0 else 0,
                "requests_per_minute": round(total / (uptime_s / 60), 2)
                if uptime_s > 0
                else 0,
                "by_endpoint": endpoints,
                "by_tier": dict(self._tier_calls),
                "by_user": dict(top_users),
                "top_smiles": dict(top_smiles),
                "recent_records": len(self._records),
            }

    def get_recent(self, limit: int = 50) -> list[dict]:
        with self._lock:
            records = self._records[-limit:]
            return [
                {
                    "timestamp": r.timestamp,
                    "method": r.method,
                    "path": r.path,
                    "status_code": r.status_code,
                    "latency_ms": round(r.latency_ms, 2),
                    "user_email": r.user_email,
                    "user_tier": r.user_tier,
                    "body_summary": r.body_summary,
                }
                for r in reversed(records)
            ]

    def get_user_usage(self, email: str) -> dict:
        with self._lock:
            user_records = [r for r in self._records if r.user_email == email]
            endpoints_used = defaultdict(int)
            smiles_used = defaultdict(int)
            total_latency = 0.0
            errors = 0
            for r in user_records:
                key = f"{r.method} {r.path}"
                endpoints_used[key] += 1
                total_latency += r.latency_ms
                if r.status_code >= 400:
                    errors += 1
                smi = r.body_summary.get("smiles")
                if smi:
                    smiles_used[smi] += 1

            return {
                "email": email,
                "total_requests": len(user_records),
                "total_errors": errors,
                "avg_latency_ms": round(total_latency / len(user_records), 2)
                if user_records
                else 0,
                "endpoints": dict(
                    sorted(endpoints_used.items(), key=lambda x: x[1], reverse=True)
                ),
                "smiles_queries": dict(
                    sorted(smiles_used.items(), key=lambda x: x[1], reverse=True)[:10]
                ),
            }


usage_tracker = UsageTracker()
