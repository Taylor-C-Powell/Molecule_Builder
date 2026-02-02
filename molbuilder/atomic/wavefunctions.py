"""
Quantum Mechanical Wave Functions for Hydrogen-like Atoms

Provides exact analytical solutions to the time-independent Schrodinger
equation for hydrogen-like (single-electron) atoms.  Multi-electron atoms
can use these with an effective nuclear charge (Z_eff) from Slater's rules.

Key components:
    - Radial wave function  R_nl(r)        (associated Laguerre polynomials)
    - Angular wave function Y_l^m(theta, phi) (spherical harmonics)
    - Full wave function    psi_nlm = R_nl * Y_l^m
    - Probability densities and expectation values

Physics conventions:
    - theta : polar angle     [0, pi]
    - phi   : azimuthal angle [0, 2*pi]
    - Condon-Shortley phase is included (via scipy)
"""

import functools
import math
import numpy as np
from scipy import special
from molbuilder.core.constants import (
    BOHR_RADIUS_M as BOHR_RADIUS,
    HBAR,
    ELECTRON_MASS,
    ELECTRON_CHARGE,
    SPEED_OF_LIGHT,
    EV_TO_JOULES,
    RYDBERG_ENERGY_EV,
)
from molbuilder.atomic.quantum_numbers import SUBSHELL_LETTER, ORBITAL_NAMES


# ===================================================================
# Radial wave function
# ===================================================================

def radial_wavefunction(n: int, l: int, r, Z: int = 1):
    """Compute R_nl(r) for a hydrogen-like atom.

    R_nl(r) = N_nl * rho^l * exp(-rho/2) * L_{n-l-1}^{2l+1}(rho)

    where rho = 2*Z*r / (n * a_0) and N_nl is the normalization constant.

    Parameters
    ----------
    n : int   — principal quantum number (>= 1)
    l : int   — angular momentum quantum number (0 <= l < n)
    r : float or ndarray — radial distance in metres
    Z : int   — nuclear charge (atomic number for hydrogen-like)

    Returns
    -------
    R : same shape as r
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if l < 0 or l >= n:
        raise ValueError(f"l must satisfy 0 <= l < n, got l={l}, n={n}")

    r = np.asarray(r, dtype=float)
    a0 = BOHR_RADIUS
    rho = 2.0 * Z * r / (n * a0)

    # Normalization in log-space to avoid overflow for large n
    log_norm = (
        1.5 * math.log(2.0 * Z / (n * a0))
        + 0.5 * (math.lgamma(n - l) - math.log(2.0 * n) - math.lgamma(n + l + 1))
    )
    norm = math.exp(log_norm)

    laguerre = special.genlaguerre(n - l - 1, 2 * l + 1)

    R = norm * np.power(rho, l) * np.exp(-rho / 2.0) * laguerre(rho)
    return R


def radial_wavefunction_scaled(n: int, l: int, r_over_a0, Z: int = 1):
    """R_nl evaluated with r given in units of a_0 (dimensionless).

    Convenient for plotting — returns R in units of a_0^{-3/2}.
    """
    r_metres = np.asarray(r_over_a0, dtype=float) * BOHR_RADIUS
    return radial_wavefunction(n, l, r_metres, Z)


# ===================================================================
# Spherical harmonics
# ===================================================================

def spherical_harmonic(l: int, m: int, theta, phi):
    """Complex spherical harmonic Y_l^m(theta, phi).

    Parameters
    ----------
    l     : int          — degree (>= 0)
    m     : int          — order (-l <= m <= l)
    theta : array-like   — polar angle [0, pi]
    phi   : array-like   — azimuthal angle [0, 2*pi]

    Returns
    -------
    Y : complex ndarray — spherical harmonic values
    """
    if abs(m) > l:
        raise ValueError(f"|m| must be <= l, got l={l}, m={m}")
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    # scipy convention: sph_harm_y(n, m, theta_polar, theta_az)
    return special.sph_harm_y(l, m, theta, phi)


def real_spherical_harmonic(l: int, m: int, theta, phi):
    """Real spherical harmonic S_l^m(theta, phi).

    These produce the standard orbital shapes:
        m > 0 : proportional to cos(m*phi)  (e.g. p_x, d_xz, d_x2-y2)
        m = 0 : axially symmetric            (e.g. p_z, d_z2)
        m < 0 : proportional to sin(|m|*phi) (e.g. p_y, d_yz, d_xy)

    Convention:
        m > 0 : sqrt(2) * (-1)^m * Re(Y_l^m)
        m = 0 : Y_l^0
        m < 0 : sqrt(2) * (-1)^|m| * Im(Y_l^|m|)
    """
    if abs(m) > l:
        raise ValueError(f"|m| must be <= l, got l={l}, m={m}")

    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)

    if m > 0:
        Y = special.sph_harm_y(l, m, theta, phi)
        return np.sqrt(2.0) * (-1)**m * Y.real
    elif m < 0:
        Y = special.sph_harm_y(l, abs(m), theta, phi)
        return np.sqrt(2.0) * (-1)**abs(m) * Y.imag
    else:
        return special.sph_harm_y(l, 0, theta, phi).real


# ===================================================================
# Full wave function
# ===================================================================

def wavefunction(n: int, l: int, m: int, r, theta, phi, Z: int = 1):
    """Full hydrogen-like wave function psi_nlm(r, theta, phi).

    psi = R_nl(r) * Y_l^m(theta, phi)

    Returns complex values (due to complex spherical harmonics).
    """
    R = radial_wavefunction(n, l, r, Z)
    Y = spherical_harmonic(l, m, theta, phi)
    return R * Y


def wavefunction_real(n: int, l: int, m: int, r, theta, phi, Z: int = 1):
    """Wave function using real spherical harmonics.

    psi = R_nl(r) * S_l^m(theta, phi)

    Returns real values. These correspond to the standard orbital shapes
    (1s, 2px, 2py, 2pz, 3dxy, etc.).
    """
    R = radial_wavefunction(n, l, r, Z)
    S = real_spherical_harmonic(l, m, theta, phi)
    return R * S


# ===================================================================
# Probability densities
# ===================================================================

def probability_density(n, l, m, r, theta, phi, Z=1):
    """Volume probability density |psi|^2 in 1/m^3."""
    psi = wavefunction(n, l, m, r, theta, phi, Z)
    return np.abs(psi)**2


def probability_density_real(n, l, m, r, theta, phi, Z=1):
    """Volume probability density using real orbital wave functions."""
    psi = wavefunction_real(n, l, m, r, theta, phi, Z)
    return psi**2


def radial_probability_density(n: int, l: int, r, Z: int = 1):
    """Radial probability density P(r) = r^2 * |R_nl(r)|^2.

    Integral of P(r) dr from 0 to infinity equals 1.
    The peak of P(r) gives the most probable radial distance.
    """
    r = np.asarray(r, dtype=float)
    R = radial_wavefunction(n, l, r, Z)
    return r**2 * np.abs(R)**2


# ===================================================================
# Expectation values (analytical, hydrogen-like)
# ===================================================================

def expectation_r(n: int, l: int, Z: int = 1) -> float:
    """<r> — mean radial distance, in metres.

    <r> = (a_0 / Z) * (1/2) * [3*n^2 - l*(l+1)]
    """
    a0 = BOHR_RADIUS
    return (a0 / Z) * 0.5 * (3 * n**2 - l * (l + 1))


def expectation_r2(n: int, l: int, Z: int = 1) -> float:
    """<r^2> — mean square radial distance, in m^2.

    <r^2> = (a_0/Z)^2 * (n^2/2) * [5*n^2 + 1 - 3*l*(l+1)]
    """
    a0 = BOHR_RADIUS
    return (a0 / Z)**2 * (n**2 / 2.0) * (5 * n**2 + 1 - 3 * l * (l + 1))


def expectation_1_over_r(n: int, l: int, Z: int = 1) -> float:
    """<1/r> in 1/metres.

    <1/r> = Z / (n^2 * a_0)
    """
    return Z / (n**2 * BOHR_RADIUS)


@functools.lru_cache(maxsize=256)
def most_probable_radius(n: int, l: int, Z: int = 1) -> float:
    """Most probable radius (peak of radial probability density), in metres.

    For l = n-1 (circular orbits): r_mp = n^2 * a_0 / Z
    General case found numerically.
    """
    if l == n - 1:
        return n**2 * BOHR_RADIUS / Z

    # Numerical search for the peak of r^2 |R_nl|^2
    from scipy.optimize import minimize_scalar
    r_mean = expectation_r(n, l, Z)
    result = minimize_scalar(
        lambda r: -radial_probability_density(n, l, r, Z),
        bounds=(1e-15, 5 * r_mean),
        method="bounded",
    )
    return result.x


# ===================================================================
# Angular momentum
# ===================================================================

def orbital_angular_momentum(l: int) -> float:
    """Magnitude of orbital angular momentum L = hbar * sqrt(l*(l+1)), in J*s."""
    return HBAR * math.sqrt(l * (l + 1))


def angular_momentum_z(m: int) -> float:
    """z-component of orbital angular momentum L_z = m * hbar, in J*s."""
    return m * HBAR


# ===================================================================
# Energy
# ===================================================================

def energy_level_eV(n: int, Z: int = 1) -> float:
    """Energy of the n-th level: E_n = -13.6 eV * Z^2 / n^2."""
    return -RYDBERG_ENERGY_EV * Z**2 / n**2


def energy_level_joules(n: int, Z: int = 1) -> float:
    """Energy of the n-th level in joules."""
    return energy_level_eV(n, Z) * EV_TO_JOULES


def transition_wavelength_nm(n_i: int, n_f: int, Z: int = 1) -> float:
    """Wavelength (nm) of a photon from transition n_i -> n_f."""
    dE = abs(energy_level_eV(n_i, Z) - energy_level_eV(n_f, Z)) * EV_TO_JOULES
    if dE == 0:
        return float("inf")
    lam = (HBAR * 2 * math.pi * SPEED_OF_LIGHT) / dE
    return lam * 1e9


# ===================================================================
# Orbital name helper
# ===================================================================

def orbital_label(n: int, l: int, m: int = None) -> str:
    """Return a human-readable orbital label, e.g. '3d' or '2p_x'."""
    base = f"{n}{SUBSHELL_LETTER.get(l, '?')}"
    if m is not None and (l, m) in ORBITAL_NAMES:
        return f"{n}{ORBITAL_NAMES[(l, m)]}"
    return base
