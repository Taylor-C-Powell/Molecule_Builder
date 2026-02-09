"""Application settings loaded from environment variables."""

from enum import Enum
from pydantic_settings import BaseSettings


class Tier(str, Enum):
    FREE = "free"
    PRO = "pro"
    TEAM = "team"
    ACADEMIC = "academic"
    ENTERPRISE = "enterprise"


class Settings(BaseSettings):
    jwt_secret_key: str = "dev-secret-change-in-production"
    jwt_algorithm: str = "HS256"
    jwt_expiry_minutes: int = 60

    rate_limit_free: int = 10
    rate_limit_pro: int = 60
    rate_limit_team: int = 120
    rate_limit_academic: int = 30
    rate_limit_enterprise: int = 300

    expensive_hourly_free: int = 5
    expensive_hourly_pro: int = 30
    expensive_hourly_team: int = 60
    expensive_hourly_academic: int = 15
    expensive_hourly_enterprise: int = 200

    molecule_store_max: int = 10_000

    audit_db_path: str = "molbuilder_audit.db"

    model_config = {"env_file": ".env", "extra": "ignore"}

    def rpm_for_tier(self, tier: Tier) -> int:
        return {
            Tier.FREE: self.rate_limit_free,
            Tier.PRO: self.rate_limit_pro,
            Tier.TEAM: self.rate_limit_team,
            Tier.ACADEMIC: self.rate_limit_academic,
            Tier.ENTERPRISE: self.rate_limit_enterprise,
        }[tier]

    def expensive_hourly_for_tier(self, tier: Tier) -> int:
        return {
            Tier.FREE: self.expensive_hourly_free,
            Tier.PRO: self.expensive_hourly_pro,
            Tier.TEAM: self.expensive_hourly_team,
            Tier.ACADEMIC: self.expensive_hourly_academic,
            Tier.ENTERPRISE: self.expensive_hourly_enterprise,
        }[tier]


settings = Settings()
