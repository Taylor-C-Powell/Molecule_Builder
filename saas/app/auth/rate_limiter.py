"""Sliding-window rate limiter (per-user and per-IP, in-memory)."""

import threading
import time
from app.config import Tier, settings


class SlidingWindowLimiter:
    """Thread-safe sliding-window rate limiter."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        # email -> list of request timestamps
        self._windows: dict[str, list[float]] = {}
        self._expensive_windows: dict[str, list[float]] = {}
        # IP -> list of timestamps (for unauthenticated endpoints)
        self._ip_windows: dict[str, dict[str, list[float]]] = {}

    def _prune(self, timestamps: list[float], window_seconds: float) -> list[float]:
        cutoff = time.monotonic() - window_seconds
        return [t for t in timestamps if t > cutoff]

    def check(self, email: str, tier: Tier) -> bool:
        """Return True if request is allowed, False if rate limited."""
        limit = settings.rpm_for_tier(tier)
        now = time.monotonic()
        with self._lock:
            ts = self._windows.get(email, [])
            ts = self._prune(ts, 60.0)
            if len(ts) >= limit:
                self._windows[email] = ts
                return False
            ts.append(now)
            self._windows[email] = ts
            return True

    def check_expensive(self, email: str, tier: Tier) -> bool:
        """Check hourly limit for expensive endpoints."""
        limit = settings.expensive_hourly_for_tier(tier)
        now = time.monotonic()
        with self._lock:
            ts = self._expensive_windows.get(email, [])
            ts = self._prune(ts, 3600.0)
            if len(ts) >= limit:
                self._expensive_windows[email] = ts
                return False
            ts.append(now)
            self._expensive_windows[email] = ts
            return True

    def check_ip(self, ip: str, endpoint: str, limit: int) -> bool:
        """Per-IP rate limit for unauthenticated endpoints (per minute)."""
        now = time.monotonic()
        with self._lock:
            ep_windows = self._ip_windows.setdefault(endpoint, {})
            ts = ep_windows.get(ip, [])
            ts = self._prune(ts, 60.0)
            if len(ts) >= limit:
                ep_windows[ip] = ts
                return False
            ts.append(now)
            ep_windows[ip] = ts
            return True


rate_limiter = SlidingWindowLimiter()
