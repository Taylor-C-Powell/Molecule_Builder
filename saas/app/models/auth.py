"""Auth request/response models."""

from pydantic import BaseModel, Field
from app.config import Tier
from app.auth.roles import Role


class APIKeyCreate(BaseModel):
    """Public self-registration: always creates free/chemist keys."""
    email: str = Field(..., description="Email for the API key owner")


class AdminKeyCreate(BaseModel):
    """Admin-only provisioning: can set any tier and role."""
    email: str = Field(..., description="Email for the API key owner")
    tier: Tier = Tier.FREE
    role: Role = Role.CHEMIST


class APIKeyResponse(BaseModel):
    api_key: str = Field(..., description="Your API key (store securely, shown once)")
    email: str
    tier: Tier
    role: Role


class TokenRequest(BaseModel):
    api_key: str


class TokenResponse(BaseModel):
    access_token: str
    token_type: str = "bearer"
    expires_in: int = Field(..., description="Seconds until expiry")


class UserInfo(BaseModel):
    email: str
    tier: Tier
    role: Role
