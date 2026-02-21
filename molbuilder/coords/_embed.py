"""Distance geometry embedding via eigendecomposition.

Samples a random distance matrix within the bounds, converts to a metric
matrix, and extracts the top 3 eigenvectors as initial 3D coordinates.
Follows with a quick violation minimization pass.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def embed(lower: np.ndarray, upper: np.ndarray, n: int,
          seed: int = 42) -> np.ndarray:
    """Embed atoms in 3D from distance bounds.

    Parameters
    ----------
    lower, upper : ndarray of shape (N, N)
        Lower and upper distance bounds.
    n : int
        Number of atoms.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    coords : ndarray of shape (N, 3)
        Initial 3D coordinates.
    """
    rng = np.random.RandomState(seed)

    # Sample random distances within bounds
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            lo = max(lower[i, j], 0.01)
            hi = max(upper[i, j], lo + 0.01)
            d = rng.uniform(lo, hi)
            D[i, j] = D[j, i] = d

    # Convert distance matrix to metric matrix (Gram matrix)
    # G = -0.5 * J * D^2 * J where J = I - 1/n * 11^T
    D2 = D ** 2
    row_mean = D2.mean(axis=1)
    col_mean = D2.mean(axis=0)
    total_mean = D2.mean()
    G = -0.5 * (D2 - row_mean[:, None] - col_mean[None, :] + total_mean)

    # Eigendecomposition: take top 3 positive eigenvalues
    eigenvalues, eigenvectors = np.linalg.eigh(G)

    # eigh returns eigenvalues in ascending order; take last 3
    idx = np.argsort(eigenvalues)[::-1][:3]
    vals = eigenvalues[idx]
    vecs = eigenvectors[:, idx]

    # Clamp negative eigenvalues to small positive
    vals = np.maximum(vals, 0.01)

    coords = vecs * np.sqrt(vals)[None, :]

    # Quick violation minimization
    coords = _minimize_violations(coords, lower, upper)

    return coords


def _minimize_violations(coords: np.ndarray, lower: np.ndarray,
                         upper: np.ndarray, max_iter: int = 100) -> np.ndarray:
    """Minimize distance bound violations with L-BFGS-B.

    Penalty function: sum of squared violations of lower and upper bounds.
    """
    n = coords.shape[0]

    # Precompute pairs
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((i, j, lower[i, j], upper[i, j]))

    def objective(x):
        pos = x.reshape(n, 3)
        energy = 0.0
        grad = np.zeros_like(pos)
        for i, j, lo, hi in pairs:
            diff = pos[j] - pos[i]
            d = np.linalg.norm(diff)
            if d < 1e-10:
                continue
            if d < lo:
                violation = lo - d
                energy += violation * violation
                f = 2.0 * violation / d * diff
                grad[i] += f
                grad[j] -= f
            elif d > hi:
                violation = d - hi
                energy += violation * violation
                f = 2.0 * violation / d * diff
                grad[i] -= f
                grad[j] += f
        return energy, grad.ravel()

    result = minimize(
        objective, coords.ravel(), method='L-BFGS-B', jac=True,
        options={'maxiter': max_iter, 'ftol': 1e-6},
    )
    return result.x.reshape(n, 3)
