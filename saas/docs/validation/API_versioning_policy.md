# MolBuilder API Versioning Policy

**Document ID**: MB-VP-001
**Effective Date**: _______________
**Version**: 1.0

## 1. Versioning Scheme

The MolBuilder API uses **Semantic Versioning** (SemVer):

```
MAJOR.MINOR.PATCH
```

| Component | Incremented When |
|-----------|-----------------|
| MAJOR | Breaking changes to existing endpoints |
| MINOR | New endpoints or backward-compatible features |
| PATCH | Bug fixes, documentation updates |

**Current version**: 1.0.0

## 2. Version Communication

Every API response includes version headers:

| Header | Description |
|--------|-------------|
| `X-API-Version` | Current API version |
| `X-Min-Supported-Version` | Oldest supported version |
| `X-Deprecation-Warning` | Present when client version is below minimum |

The version endpoint (`GET /api/v1/version`) returns full version details and this policy.

## 3. Deprecation Policy

### 3.1 Timeline

- **Minimum support window**: 12 months from the release of a new MAJOR version
- **Deprecation notice**: Communicated via `X-Deprecation-Warning` header and changelog
- **End of life**: After 12-month window, deprecated versions may return `410 Gone`

### 3.2 Deprecation Process

1. New MAJOR version released (e.g., v2.0.0)
2. Previous version (v1.x.x) enters deprecation immediately
3. All requests to deprecated version receive `X-Deprecation-Warning` header
4. 6-month reminder: email notification to all registered API key holders
5. 11-month final notice: last-chance email notification
6. 12 months: deprecated version endpoints may be removed

### 3.3 Breaking Changes

The following are considered breaking changes (require MAJOR version bump):

- Removing an endpoint
- Changing the response schema of an existing endpoint
- Changing required request parameters
- Changing authentication mechanisms
- Changing error response formats

The following are NOT breaking changes:

- Adding new endpoints
- Adding optional request parameters
- Adding new fields to response objects
- Performance improvements
- Bug fixes that correct behavior to match documentation

## 4. Migration Guide Template

When a new MAJOR version is released, a migration guide will be published covering:

1. **Summary of breaking changes**
2. **Before/after examples** for each changed endpoint
3. **Timeline** with key dates
4. **Support contacts** for migration assistance

### Example Migration Steps

```
# Update base URL
- OLD: https://api.molbuilder.com/api/v1/...
+ NEW: https://api.molbuilder.com/api/v2/...

# Update changed request formats (example)
- OLD: POST /molecule/from-smiles {"smiles": "CCO"}
+ NEW: POST /molecule/parse {"smiles": "CCO", "format": "smiles"}
```

## 5. Client Best Practices

1. **Always send** `X-API-Version` header with your requests
2. **Monitor** response headers for `X-Deprecation-Warning`
3. **Subscribe** to the MolBuilder changelog for version announcements
4. **Test** against new versions in staging before upgrading production
5. **Plan migrations** within the first 6 months of a deprecation notice

## 6. Changelog

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2026-02-09 | Initial release: auth, molecule, retrosynthesis, process, elements, analytics, audit, RBAC |

## 7. Contact

For versioning questions or migration assistance:
- API Documentation: `/api/v1/version`
- GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
