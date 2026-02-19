"""Legal endpoints: Terms of Service and Privacy Policy."""

from fastapi import APIRouter
from fastapi.responses import JSONResponse

router = APIRouter(prefix="/api/v1/legal", tags=["legal"])

TERMS_OF_SERVICE = {
    "title": "MolBuilder Terms of Service",
    "effective_date": "2026-02-19",
    "sections": [
        {
            "heading": "1. Acceptance of Terms",
            "content": (
                "By accessing or using the MolBuilder API, web dashboard, studio editor, "
                "or any associated services (collectively, the 'Service'), you agree to be "
                "bound by these Terms of Service. If you are using the Service on behalf of "
                "an organization, you represent that you have authority to bind that organization."
            ),
        },
        {
            "heading": "2. Description of Service",
            "content": (
                "MolBuilder provides computational chemistry tools including SMILES parsing, "
                "molecular property calculation, retrosynthetic analysis, process engineering "
                "evaluation, and 3D molecular visualization. The Service is offered via a REST "
                "API with tiered access (Free, Pro, Team, Academic, Enterprise)."
            ),
        },
        {
            "heading": "3. Account Registration and API Keys",
            "content": (
                "You must provide a valid email address to register. You are responsible for "
                "maintaining the confidentiality of your API key. You must not share API keys "
                "or attempt to circumvent rate limits. We reserve the right to revoke keys that "
                "violate these terms."
            ),
        },
        {
            "heading": "4. Acceptable Use",
            "content": (
                "You may use the Service for lawful purposes including academic research, "
                "commercial drug discovery, process chemistry planning, and education. You "
                "must not: (a) use the Service to facilitate the synthesis of controlled "
                "substances, chemical weapons, or other prohibited materials; (b) attempt to "
                "reverse-engineer, decompile, or extract the source code of proprietary "
                "components; (c) exceed rate limits or abuse the API; (d) resell API access "
                "without written authorization."
            ),
        },
        {
            "heading": "5. Intellectual Property",
            "content": (
                "You retain all rights to molecular structures, SMILES strings, and data you "
                "submit to the Service. MolBuilder does not claim ownership of your input data "
                "or the results generated from it. The MolBuilder core library is open source "
                "under the MIT License. The API service, web dashboard, and studio editor are "
                "proprietary to Materia Foundation."
            ),
        },
        {
            "heading": "6. Data Handling and Scientific Accuracy",
            "content": (
                "MolBuilder uses rule-based algorithms and empirical data for molecular "
                "property predictions, retrosynthetic analysis, and process engineering "
                "estimates. Results are computational approximations and must not be used as "
                "the sole basis for safety-critical decisions. Always validate computational "
                "results with experimental data and consult qualified professionals before "
                "proceeding with synthesis, manufacturing, or clinical applications. "
                "MolBuilder is not a substitute for professional chemical engineering judgment."
            ),
        },
        {
            "heading": "7. Billing and Subscriptions",
            "content": (
                "Paid subscriptions are billed through Stripe. Pro plans may be billed monthly "
                "or annually. You may cancel at any time; access continues through the end of "
                "the current billing period. Refunds are handled on a case-by-case basis. "
                "Pricing is subject to change with 30 days notice."
            ),
        },
        {
            "heading": "8. Service Level and Availability",
            "content": (
                "We aim for high availability but do not guarantee uninterrupted service. "
                "Scheduled maintenance will be announced in advance when possible. Free-tier "
                "users have no SLA. Pro and Enterprise users should contact us for SLA terms."
            ),
        },
        {
            "heading": "9. Limitation of Liability",
            "content": (
                "THE SERVICE IS PROVIDED 'AS IS' WITHOUT WARRANTIES OF ANY KIND. MOLBUILDER "
                "AND MATERIA FOUNDATION SHALL NOT BE LIABLE FOR ANY INDIRECT, INCIDENTAL, "
                "SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING FROM USE OF THE SERVICE, INCLUDING "
                "BUT NOT LIMITED TO DAMAGES FROM RELIANCE ON COMPUTATIONAL RESULTS FOR "
                "SYNTHESIS PLANNING, SAFETY ASSESSMENTS, OR COST ESTIMATES. OUR TOTAL "
                "LIABILITY SHALL NOT EXCEED THE AMOUNT YOU PAID FOR THE SERVICE IN THE "
                "PRECEDING 12 MONTHS."
            ),
        },
        {
            "heading": "10. Compliance",
            "content": (
                "The Service includes 21 CFR Part 11 compliant audit trail functionality "
                "for users who require it. However, achieving full regulatory compliance is "
                "the responsibility of the user and their organization. MolBuilder provides "
                "tools to support compliance but does not guarantee it."
            ),
        },
        {
            "heading": "11. Termination",
            "content": (
                "We may suspend or terminate your access for violation of these terms. You may "
                "request deletion of your account and associated data at any time via the "
                "API (DELETE /api/v1/auth/me). Upon termination, your API keys are revoked "
                "and stored data is deleted within 30 days."
            ),
        },
        {
            "heading": "12. Governing Law",
            "content": (
                "These terms are governed by the laws of the State of North Carolina, "
                "United States, without regard to conflict of law provisions."
            ),
        },
        {
            "heading": "13. Contact",
            "content": (
                "For questions about these terms, contact: legal@molbuilder.io or open an "
                "issue at https://github.com/Taylor-C-Powell/Molecule_Builder/issues"
            ),
        },
    ],
}

PRIVACY_POLICY = {
    "title": "MolBuilder Privacy Policy",
    "effective_date": "2026-02-19",
    "sections": [
        {
            "heading": "1. Information We Collect",
            "content": (
                "We collect: (a) Email address, provided at registration; (b) API usage "
                "data including endpoints called, timestamps, IP addresses, and response "
                "codes; (c) Molecular data (SMILES strings, molecule names) submitted through "
                "the API during your session; (d) Payment information processed by Stripe "
                "(we do not store card numbers)."
            ),
        },
        {
            "heading": "2. How We Use Your Data",
            "content": (
                "We use your data to: (a) Authenticate API requests and manage your account; "
                "(b) Enforce rate limits and prevent abuse; (c) Generate usage analytics "
                "(visible to administrators only); (d) Maintain 21 CFR Part 11 audit trails "
                "when required; (e) Process billing through Stripe; (f) Improve the Service."
            ),
        },
        {
            "heading": "3. Molecular Data Handling",
            "content": (
                "SMILES strings and molecular structures you submit are processed in memory "
                "for computation and may be temporarily stored for session-based retrieval "
                "(e.g., molecule ID lookups). We do not use your molecular data to train "
                "models, share with third parties, or retain beyond your session unless you "
                "explicitly save molecules to your library. We recognize that molecular "
                "structures may constitute trade secrets or proprietary research data and "
                "treat them accordingly."
            ),
        },
        {
            "heading": "4. Data Retention",
            "content": (
                "Account data (email, API key hashes) is retained while your account is "
                "active. Audit trail records are retained for the minimum period required "
                "by applicable regulations (typically 2 years for 21 CFR Part 11). Usage "
                "analytics are retained for 12 months. You may request deletion of your "
                "account and all associated data at any time."
            ),
        },
        {
            "heading": "5. Data Security",
            "content": (
                "API keys are stored as HMAC-SHA256 hashes, not in plaintext. All API "
                "traffic is encrypted via TLS. JWT tokens expire after 60 minutes. "
                "Audit trail records are cryptographically signed to prevent tampering. "
                "Access to user data is restricted to administrators via role-based access "
                "control (RBAC)."
            ),
        },
        {
            "heading": "6. Third-Party Services",
            "content": (
                "We use: (a) Railway for API hosting; (b) Vercel for frontend hosting; "
                "(c) Stripe for payment processing; (d) GitHub for source code hosting. "
                "Each service has its own privacy policy. We do not sell your data to any "
                "third party."
            ),
        },
        {
            "heading": "7. Your Rights (GDPR and Global Privacy)",
            "content": (
                "Regardless of your location, you have the right to: (a) Access your data "
                "(GET /api/v1/auth/me/export); (b) Delete your account and data "
                "(DELETE /api/v1/auth/me); (c) Know what data we hold about you; "
                "(d) Object to processing of your data. For EU/EEA residents, processing "
                "is based on contractual necessity (Article 6(1)(b) GDPR) and legitimate "
                "interest (Article 6(1)(f) GDPR). For scientific research use, additional "
                "protections under GDPR Recital 159 may apply."
            ),
        },
        {
            "heading": "8. Cookies",
            "content": (
                "The MolBuilder web dashboard uses localStorage to persist authentication "
                "state (API key and JWT token). We do not use tracking cookies or third-party "
                "analytics cookies."
            ),
        },
        {
            "heading": "9. Children",
            "content": (
                "The Service is not directed at children under 16. We do not knowingly "
                "collect data from children."
            ),
        },
        {
            "heading": "10. Changes to This Policy",
            "content": (
                "We will notify registered users via email of material changes to this "
                "policy at least 14 days before they take effect."
            ),
        },
        {
            "heading": "11. Contact",
            "content": (
                "For privacy inquiries: privacy@molbuilder.io or open an issue at "
                "https://github.com/Taylor-C-Powell/Molecule_Builder/issues"
            ),
        },
    ],
}


@router.get("/terms")
def get_terms():
    """Returns the Terms of Service as structured JSON."""
    return JSONResponse(content=TERMS_OF_SERVICE)


@router.get("/privacy")
def get_privacy():
    """Returns the Privacy Policy as structured JSON."""
    return JSONResponse(content=PRIVACY_POLICY)
