"""
Unit tests for molbuilder.atomic modules:
    bohr, quantum_numbers, quantum_atom, wavefunctions
"""

import math
import unittest

import numpy as np


# ---------------------------------------------------------------------------
# bohr  --  BohrAtom
# ---------------------------------------------------------------------------
class TestBohrAtom(unittest.TestCase):
    """Verify the Bohr model calculations for hydrogen-like atoms."""

    # --- construction ---
    def test_hydrogen_creation(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        self.assertEqual(h.atomic_number, 1)
        self.assertEqual(h.symbol, "H")
        self.assertEqual(h.name, "Hydrogen")
        self.assertEqual(h.num_electrons, 1)

    def test_helium_creation(self):
        from molbuilder.atomic.bohr import BohrAtom
        he = BohrAtom(2)
        self.assertEqual(he.symbol, "He")
        self.assertEqual(he.num_electrons, 2)

    def test_invalid_atomic_number_raises(self):
        from molbuilder.atomic.bohr import BohrAtom
        with self.assertRaises(ValueError):
            BohrAtom(0)
        with self.assertRaises(ValueError):
            BohrAtom(119)

    def test_charge_reduces_electrons(self):
        from molbuilder.atomic.bohr import BohrAtom
        na_ion = BohrAtom(11, charge=1)
        self.assertEqual(na_ion.num_electrons, 10)

    def test_excessive_charge_raises(self):
        from molbuilder.atomic.bohr import BohrAtom
        with self.assertRaises(ValueError):
            BohrAtom(1, charge=2)

    # --- ionization energy ---
    def test_hydrogen_ionization_energy(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        # Bohr model: ionization of H is 13.6 eV
        self.assertAlmostEqual(h.ionization_energy(), 13.6, delta=0.1)

    def test_he_plus_ionization_energy(self):
        from molbuilder.atomic.bohr import BohrAtom
        # He+ is hydrogen-like with Z=2; ionization from n=1: 13.6*4 = 54.4 eV
        he_plus = BohrAtom(2, charge=1)
        self.assertAlmostEqual(he_plus.ionization_energy(), 54.4, delta=0.2)

    # --- orbital radii ---
    def test_hydrogen_n1_radius(self):
        from molbuilder.atomic.bohr import BohrAtom
        from molbuilder.core.constants import BOHR_RADIUS_M
        h = BohrAtom(1)
        self.assertAlmostEqual(h.orbital_radius(1), BOHR_RADIUS_M, places=15)

    def test_orbital_radii_scale_with_n_squared(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        r1 = h.orbital_radius(1)
        r2 = h.orbital_radius(2)
        r3 = h.orbital_radius(3)
        self.assertAlmostEqual(r2 / r1, 4.0, places=10)
        self.assertAlmostEqual(r3 / r1, 9.0, places=10)

    def test_orbital_radius_pm(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        r_pm = h.orbital_radius_pm(1)
        r_m = h.orbital_radius(1)
        self.assertAlmostEqual(r_pm, r_m * 1e12, places=3)

    def test_orbital_radius_invalid_n_raises(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        with self.assertRaises(ValueError):
            h.orbital_radius(0)

    # --- energy levels ---
    def test_hydrogen_energy_n1(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        self.assertAlmostEqual(h.energy_level(1), -13.6, delta=0.1)

    def test_hydrogen_energy_n2(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        self.assertAlmostEqual(h.energy_level(2), -3.4, delta=0.1)

    def test_energy_level_joules(self):
        from molbuilder.atomic.bohr import BohrAtom
        from molbuilder.core.constants import EV_TO_JOULES
        h = BohrAtom(1)
        e_eV = h.energy_level(1)
        e_J = h.energy_level_joules(1)
        self.assertAlmostEqual(e_J, e_eV * EV_TO_JOULES, places=25)

    # --- transitions ---
    def test_transition_energy_hydrogen_2_to_1(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        # E(2) - E(1) = -3.4 - (-13.6) = 10.2 eV (photon emitted)
        delta_e = h.transition_energy(2, 1)
        self.assertAlmostEqual(delta_e, 10.2, delta=0.1)

    def test_transition_energy_absorption_negative(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        # Absorption: n=1 -> n=2 gives negative transition_energy
        delta_e = h.transition_energy(1, 2)
        self.assertLess(delta_e, 0)

    # --- Balmer series wavelengths ---
    def test_balmer_h_alpha(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        # H-alpha: n=3 -> n=2, ~656 nm
        wl_nm = h.transition_wavelength_nm(3, 2)
        self.assertAlmostEqual(wl_nm, 656.0, delta=5.0)

    def test_balmer_h_beta(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        # H-beta: n=4 -> n=2, ~486 nm
        wl_nm = h.transition_wavelength_nm(4, 2)
        self.assertAlmostEqual(wl_nm, 486.0, delta=5.0)

    def test_balmer_h_gamma(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        # H-gamma: n=5 -> n=2, ~434 nm
        wl_nm = h.transition_wavelength_nm(5, 2)
        self.assertAlmostEqual(wl_nm, 434.0, delta=5.0)

    def test_transition_wavelength_same_level_inf(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        self.assertEqual(h.transition_wavelength(2, 2), float('inf'))

    # --- orbital velocity and period ---
    def test_orbital_velocity_positive(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        v = h.orbital_velocity(1)
        self.assertGreater(v, 0)
        # Should be on the order of ~2.2e6 m/s for hydrogen n=1
        self.assertAlmostEqual(v, 2.2e6, delta=0.1e6)

    def test_orbital_period_positive(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        T = h.orbital_period(1)
        self.assertGreater(T, 0)

    # --- coulomb force ---
    def test_coulomb_force_positive(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        f = h.coulomb_force_on_electron(1)
        self.assertGreater(f, 0)

    # --- shell configuration ---
    def test_hydrogen_shell_config(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        self.assertEqual(h.shell_config, [1])

    def test_neon_shell_config(self):
        from molbuilder.atomic.bohr import BohrAtom
        ne = BohrAtom(10)
        self.assertEqual(ne.shell_config, [2, 8])

    def test_sodium_shell_config(self):
        from molbuilder.atomic.bohr import BohrAtom
        na = BohrAtom(11)
        self.assertEqual(na.shell_config, [2, 8, 1])

    # --- convenience constructors ---
    def test_from_symbol(self):
        from molbuilder.atomic.bohr import from_symbol
        h = from_symbol("H")
        self.assertEqual(h.atomic_number, 1)

    def test_from_name(self):
        from molbuilder.atomic.bohr import from_name
        fe = from_name("Iron")
        self.assertEqual(fe.atomic_number, 26)

    # --- repr ---
    def test_repr(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        r = repr(h)
        self.assertIn("BohrAtom", r)
        self.assertIn("H", r)

    def test_summary_returns_string(self):
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        s = h.summary()
        self.assertIsInstance(s, str)
        self.assertIn("Hydrogen", s)


# ---------------------------------------------------------------------------
# quantum_numbers
# ---------------------------------------------------------------------------
class TestQuantumNumbers(unittest.TestCase):
    """Verify quantum number containers and aufbau ordering."""

    # --- aufbau_order ---
    def test_aufbau_first_five(self):
        from molbuilder.atomic.quantum_numbers import aufbau_order, SUBSHELL_LETTER
        order = aufbau_order()
        labels = [f"{n}{SUBSHELL_LETTER[l]}" for n, l in order]
        self.assertEqual(labels[0], "1s")
        self.assertEqual(labels[1], "2s")
        self.assertEqual(labels[2], "2p")
        self.assertEqual(labels[3], "3s")
        self.assertEqual(labels[4], "3p")

    def test_aufbau_includes_4s_before_3d(self):
        from molbuilder.atomic.quantum_numbers import aufbau_order
        order = aufbau_order()
        idx_4s = order.index((4, 0))
        idx_3d = order.index((3, 2))
        self.assertLess(idx_4s, idx_3d)

    def test_aufbau_length(self):
        from molbuilder.atomic.quantum_numbers import aufbau_order
        # For max_n=7, the sum is 1+2+3+4+5+6+7 = 28 subshells
        order = aufbau_order(7)
        self.assertEqual(len(order), 28)

    # --- QuantumState validation ---
    def test_valid_quantum_state(self):
        from molbuilder.atomic.quantum_numbers import QuantumState
        qs = QuantumState(2, 1, 0, 0.5)
        self.assertEqual(qs.n, 2)
        self.assertEqual(qs.l, 1)
        self.assertEqual(qs.ml, 0)
        self.assertAlmostEqual(qs.ms, 0.5)

    def test_quantum_state_invalid_n(self):
        from molbuilder.atomic.quantum_numbers import QuantumState
        with self.assertRaises(ValueError):
            QuantumState(0, 0, 0, 0.5)

    def test_quantum_state_invalid_l_too_large(self):
        from molbuilder.atomic.quantum_numbers import QuantumState
        with self.assertRaises(ValueError):
            QuantumState(1, 1, 0, 0.5)  # l must be < n

    def test_quantum_state_invalid_l_negative(self):
        from molbuilder.atomic.quantum_numbers import QuantumState
        with self.assertRaises(ValueError):
            QuantumState(2, -1, 0, 0.5)

    def test_quantum_state_invalid_ml(self):
        from molbuilder.atomic.quantum_numbers import QuantumState
        with self.assertRaises(ValueError):
            QuantumState(2, 1, 2, 0.5)  # ml must be in [-l, +l]

    def test_quantum_state_invalid_ms(self):
        from molbuilder.atomic.quantum_numbers import QuantumState
        with self.assertRaises(ValueError):
            QuantumState(1, 0, 0, 0.0)  # ms must be +/-0.5

    def test_quantum_state_subshell_label(self):
        from molbuilder.atomic.quantum_numbers import QuantumState
        qs = QuantumState(3, 2, 1, 0.5)
        self.assertEqual(qs.subshell_label, "3d")

    # --- Subshell ---
    def test_subshell_s_capacity(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        s = Subshell(1, 0)
        self.assertEqual(s.capacity, 2)

    def test_subshell_p_capacity(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        p = Subshell(2, 1)
        self.assertEqual(p.capacity, 6)

    def test_subshell_d_capacity(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        d = Subshell(3, 2)
        self.assertEqual(d.capacity, 10)

    def test_subshell_f_capacity(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        f = Subshell(4, 3)
        self.assertEqual(f.capacity, 14)

    def test_subshell_is_full(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        s = Subshell(1, 0, electron_count=2)
        self.assertTrue(s.is_full)

    def test_subshell_not_full(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        s = Subshell(1, 0, electron_count=1)
        self.assertFalse(s.is_full)

    def test_subshell_is_half_filled(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        d = Subshell(3, 2, electron_count=5)
        self.assertTrue(d.is_half_filled)

    def test_subshell_label(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        ss = Subshell(4, 3, electron_count=7)
        self.assertEqual(ss.label, "4f")

    def test_subshell_quantum_states_count(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        ss = Subshell(2, 1, electron_count=4)
        states = ss.quantum_states()
        self.assertEqual(len(states), 4)

    def test_subshell_quantum_states_hund(self):
        from molbuilder.atomic.quantum_numbers import Subshell
        # 3 electrons in 2p: should fill ml=-1,0,+1 all spin-up first
        ss = Subshell(2, 1, electron_count=3)
        states = ss.quantum_states()
        self.assertEqual(len(states), 3)
        for qs in states:
            self.assertAlmostEqual(qs.ms, 0.5)

    # --- AUFBAU_EXCEPTIONS ---
    def test_aufbau_exceptions_cr(self):
        from molbuilder.atomic.quantum_numbers import AUFBAU_EXCEPTIONS
        self.assertIn(24, AUFBAU_EXCEPTIONS)
        # Cr: should have 3d5 4s1
        cr_config = AUFBAU_EXCEPTIONS[24]
        d_count = sum(c for n, l, c in cr_config if n == 3 and l == 2)
        s_count = sum(c for n, l, c in cr_config if n == 4 and l == 0)
        self.assertEqual(d_count, 5)
        self.assertEqual(s_count, 1)

    def test_aufbau_exceptions_cu(self):
        from molbuilder.atomic.quantum_numbers import AUFBAU_EXCEPTIONS
        self.assertIn(29, AUFBAU_EXCEPTIONS)
        # Cu: should have 3d10 4s1
        cu_config = AUFBAU_EXCEPTIONS[29]
        d_count = sum(c for n, l, c in cu_config if n == 3 and l == 2)
        s_count = sum(c for n, l, c in cu_config if n == 4 and l == 0)
        self.assertEqual(d_count, 10)
        self.assertEqual(s_count, 1)


# ---------------------------------------------------------------------------
# quantum_atom  --  QuantumAtom
# ---------------------------------------------------------------------------
class TestQuantumAtom(unittest.TestCase):
    """Verify quantum mechanical atom model: configurations, valence, Slater."""

    # --- electron configurations ---
    def test_hydrogen_config(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        h = QuantumAtom(1)
        self.assertEqual(h.electron_configuration_string(), "1s1")

    def test_carbon_config(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        self.assertEqual(c.electron_configuration_string(), "1s2 2s2 2p2")

    def test_nitrogen_config(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        n = QuantumAtom(7)
        self.assertEqual(n.electron_configuration_string(), "1s2 2s2 2p3")

    def test_iron_config(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        fe = QuantumAtom(26)
        config = fe.electron_configuration_string()
        self.assertIn("3d6", config)
        self.assertIn("4s2", config)

    def test_chromium_config_exception(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        cr = QuantumAtom(24)
        config = cr.electron_configuration_string()
        # Cr is a known exception: 3d5 4s1 (half-filled stability)
        self.assertIn("3d5", config)
        self.assertIn("4s1", config)

    def test_copper_config_exception(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        cu = QuantumAtom(29)
        config = cu.electron_configuration_string()
        # Cu is a known exception: 3d10 4s1 (fully-filled stability)
        self.assertIn("3d10", config)
        self.assertIn("4s1", config)

    def test_total_electrons_sum(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        for z in [1, 6, 26, 29]:
            atom = QuantumAtom(z)
            total = sum(ss.electron_count for ss in atom.subshells)
            self.assertEqual(total, z, f"Electron count mismatch for Z={z}")

    # --- valence_electrons (property, NOT a method call) ---
    def test_hydrogen_valence_electrons(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        h = QuantumAtom(1)
        self.assertEqual(h.valence_electrons, 1)

    def test_carbon_valence_electrons(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        # Valence shell is n=2: 2s2 + 2p2 = 4
        self.assertEqual(c.valence_electrons, 4)

    def test_sodium_valence_electrons(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        na = QuantumAtom(11)
        # Valence shell is n=3: 3s1
        self.assertEqual(na.valence_electrons, 1)

    def test_valence_electrons_is_property(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        # Confirm it is a property (not callable)
        self.assertFalse(callable(c.valence_electrons))

    def test_valence_shell(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        self.assertEqual(c.valence_shell, 2)

    def test_core_electrons(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        self.assertEqual(c.core_electrons, 2)  # 1s2

    # --- noble gas notation ---
    def test_noble_gas_notation_carbon(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        notation = c.noble_gas_notation()
        self.assertIn("[He]", notation)
        self.assertIn("2s2", notation)
        self.assertIn("2p2", notation)

    def test_noble_gas_notation_iron(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        fe = QuantumAtom(26)
        notation = fe.noble_gas_notation()
        self.assertIn("[Ar]", notation)
        self.assertIn("3d6", notation)
        self.assertIn("4s2", notation)

    def test_noble_gas_notation_hydrogen_no_core(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        h = QuantumAtom(1)
        notation = h.noble_gas_notation()
        # No noble gas core for H, should be same as full config
        self.assertEqual(notation, "1s1")

    # --- Slater effective nuclear charge ---
    def test_slater_zeff_hydrogen(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        h = QuantumAtom(1)
        z_eff = h.effective_nuclear_charge(1, 0)
        # Hydrogen has 1 electron, so Z_eff = Z - 0 (no shielding, minus self correction 0.30)
        # Slater: same-group shielding = (1-1)*0.30 = 0 => Z_eff = 1.0
        self.assertAlmostEqual(z_eff, 1.0, places=5)

    def test_slater_zeff_helium_1s(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        he = QuantumAtom(2)
        z_eff = he.effective_nuclear_charge(1, 0)
        # He 1s2: shielding = 1 * 0.30 => Z_eff = 2 - 0.30 = 1.70
        self.assertAlmostEqual(z_eff, 1.70, places=2)

    def test_slater_zeff_positive(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        z_eff = c.effective_nuclear_charge(2, 1)
        self.assertGreater(z_eff, 0)

    # --- orbital energy ---
    def test_orbital_energy_negative(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        e = c.orbital_energy_eV(2, 1)
        self.assertLess(e, 0)

    # --- ionization energy ---
    def test_ionization_energy_positive(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        ie = c.ionization_energy_eV()
        self.assertGreater(ie, 0)

    # --- spin ---
    def test_nitrogen_spin_multiplicity(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        n = QuantumAtom(7)
        # N: 2p3 (3 unpaired electrons) => S = 3/2, 2S+1 = 4
        self.assertEqual(n.spin_multiplicity(), 4)

    def test_neon_spin_multiplicity(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        ne = QuantumAtom(10)
        # Ne: all paired, S = 0, 2S+1 = 1
        self.assertEqual(ne.spin_multiplicity(), 1)

    # --- charged atoms ---
    def test_sodium_cation(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        na_plus = QuantumAtom(11, charge=1)
        self.assertEqual(na_plus.num_electrons, 10)
        config = na_plus.electron_configuration_string()
        # Should look like neon
        self.assertIn("2p6", config)

    # --- config_tuples ---
    def test_config_tuples_carbon(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        tuples = c.config_tuples
        self.assertEqual(tuples, [(1, 0, 2), (2, 0, 2), (2, 1, 2)])

    # --- quantum states enumeration ---
    def test_quantum_states_count(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        self.assertEqual(len(c.quantum_states), 6)

    # --- convenience constructors ---
    def test_from_symbol(self):
        from molbuilder.atomic.quantum_atom import from_symbol
        fe = from_symbol("Fe")
        self.assertEqual(fe.atomic_number, 26)

    def test_from_name(self):
        from molbuilder.atomic.quantum_atom import from_name
        cu = from_name("Copper")
        self.assertEqual(cu.atomic_number, 29)

    # --- repr ---
    def test_repr(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        r = repr(c)
        self.assertIn("QuantumAtom", r)
        self.assertIn("C", r)

    def test_summary_returns_string(self):
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        s = c.summary()
        self.assertIsInstance(s, str)
        self.assertIn("Carbon", s)


# ---------------------------------------------------------------------------
# wavefunctions
# ---------------------------------------------------------------------------
class TestWavefunctions(unittest.TestCase):
    """Verify quantum wavefunctions: radial, spherical harmonics, energy."""

    # --- radial_wavefunction ---
    def test_radial_wavefunction_nonzero(self):
        from molbuilder.atomic.wavefunctions import radial_wavefunction
        from molbuilder.core.constants import BOHR_RADIUS_M
        r = BOHR_RADIUS_M  # 1 Bohr radius
        R_val = radial_wavefunction(1, 0, r, Z=1)
        self.assertNotAlmostEqual(float(R_val), 0.0)

    def test_radial_wavefunction_1s_at_origin_near_max(self):
        from molbuilder.atomic.wavefunctions import radial_wavefunction
        # R_10(0) should be maximum (for 1s), evaluate at very small r
        R_origin = radial_wavefunction(1, 0, 1e-15, Z=1)
        R_far = radial_wavefunction(1, 0, 1e-9, Z=1)
        self.assertGreater(abs(float(R_origin)), abs(float(R_far)))

    def test_radial_wavefunction_array_input(self):
        from molbuilder.atomic.wavefunctions import radial_wavefunction
        from molbuilder.core.constants import BOHR_RADIUS_M
        r_arr = np.linspace(0.1 * BOHR_RADIUS_M, 5 * BOHR_RADIUS_M, 50)
        R_arr = radial_wavefunction(1, 0, r_arr, Z=1)
        self.assertEqual(R_arr.shape, (50,))

    def test_radial_wavefunction_invalid_n(self):
        from molbuilder.atomic.wavefunctions import radial_wavefunction
        with self.assertRaises(ValueError):
            radial_wavefunction(0, 0, 1e-10, Z=1)

    def test_radial_wavefunction_invalid_l(self):
        from molbuilder.atomic.wavefunctions import radial_wavefunction
        with self.assertRaises(ValueError):
            radial_wavefunction(1, 1, 1e-10, Z=1)  # l must be < n

    def test_radial_wavefunction_2p_zero_at_origin(self):
        from molbuilder.atomic.wavefunctions import radial_wavefunction
        # For l > 0, R_nl(0) = 0 (r^l factor)
        R_at_zero = radial_wavefunction(2, 1, 0.0, Z=1)
        self.assertAlmostEqual(float(R_at_zero), 0.0, places=10)

    # --- energy_level_eV ---
    def test_energy_level_hydrogen_n1(self):
        from molbuilder.atomic.wavefunctions import energy_level_eV
        e = energy_level_eV(1, Z=1)
        self.assertAlmostEqual(e, -13.6, delta=0.1)

    def test_energy_level_hydrogen_n2(self):
        from molbuilder.atomic.wavefunctions import energy_level_eV
        e = energy_level_eV(2, Z=1)
        expected = -13.6 / 4.0
        self.assertAlmostEqual(e, expected, delta=0.1)

    def test_energy_level_hydrogen_n3(self):
        from molbuilder.atomic.wavefunctions import energy_level_eV
        e = energy_level_eV(3, Z=1)
        expected = -13.6 / 9.0
        self.assertAlmostEqual(e, expected, delta=0.1)

    def test_energy_level_scales_with_Z_squared(self):
        from molbuilder.atomic.wavefunctions import energy_level_eV
        e_h = energy_level_eV(1, Z=1)
        e_he = energy_level_eV(1, Z=2)
        self.assertAlmostEqual(e_he / e_h, 4.0, places=5)

    def test_energy_level_joules_conversion(self):
        from molbuilder.atomic.wavefunctions import energy_level_eV, energy_level_joules
        from molbuilder.core.constants import EV_TO_JOULES
        e_eV = energy_level_eV(1, Z=1)
        e_J = energy_level_joules(1, Z=1)
        self.assertAlmostEqual(e_J, e_eV * EV_TO_JOULES, places=25)

    # --- transition_wavelength_nm ---
    def test_transition_wavelength_lyman_alpha(self):
        from molbuilder.atomic.wavefunctions import transition_wavelength_nm
        # Lyman alpha: n=2 -> n=1, ~121.5 nm
        wl = transition_wavelength_nm(2, 1, Z=1)
        self.assertAlmostEqual(wl, 121.5, delta=2.0)

    def test_transition_wavelength_balmer_alpha(self):
        from molbuilder.atomic.wavefunctions import transition_wavelength_nm
        # Balmer alpha: n=3 -> n=2, ~656 nm
        wl = transition_wavelength_nm(3, 2, Z=1)
        self.assertAlmostEqual(wl, 656.0, delta=5.0)

    # --- spherical harmonics ---
    def test_spherical_harmonic_Y00_constant(self):
        from molbuilder.atomic.wavefunctions import spherical_harmonic
        # Y_0^0 = 1/(2*sqrt(pi)) -- constant for all angles
        Y = spherical_harmonic(0, 0, 0.5, 0.5)
        expected = 1.0 / (2.0 * math.sqrt(math.pi))
        self.assertAlmostEqual(abs(complex(Y)), expected, places=8)

    def test_spherical_harmonic_invalid_m_raises(self):
        from molbuilder.atomic.wavefunctions import spherical_harmonic
        with self.assertRaises(ValueError):
            spherical_harmonic(1, 2, 0.0, 0.0)

    # --- real_spherical_harmonic ---
    def test_real_spherical_harmonic_s_orbital(self):
        from molbuilder.atomic.wavefunctions import real_spherical_harmonic
        # m=0, l=0: should be same as Y_0^0
        S = real_spherical_harmonic(0, 0, 0.5, 0.5)
        expected = 1.0 / (2.0 * math.sqrt(math.pi))
        self.assertAlmostEqual(float(S), expected, places=8)

    # --- probability densities ---
    def test_probability_density_nonnegative(self):
        from molbuilder.atomic.wavefunctions import probability_density
        from molbuilder.core.constants import BOHR_RADIUS_M
        pd = probability_density(1, 0, 0, BOHR_RADIUS_M, 0.5, 0.5, Z=1)
        self.assertGreaterEqual(float(pd), 0.0)

    def test_radial_probability_density_nonnegative(self):
        from molbuilder.atomic.wavefunctions import radial_probability_density
        from molbuilder.core.constants import BOHR_RADIUS_M
        rpd = radial_probability_density(1, 0, BOHR_RADIUS_M, Z=1)
        self.assertGreaterEqual(float(rpd), 0.0)

    # --- angular momentum ---
    def test_orbital_angular_momentum_s(self):
        from molbuilder.atomic.wavefunctions import orbital_angular_momentum
        # l=0: L = 0
        L = orbital_angular_momentum(0)
        self.assertAlmostEqual(L, 0.0, places=40)

    def test_orbital_angular_momentum_p(self):
        from molbuilder.atomic.wavefunctions import orbital_angular_momentum
        from molbuilder.core.constants import HBAR
        # l=1: L = hbar * sqrt(2)
        L = orbital_angular_momentum(1)
        self.assertAlmostEqual(L, HBAR * math.sqrt(2), places=45)

    def test_angular_momentum_z(self):
        from molbuilder.atomic.wavefunctions import angular_momentum_z
        from molbuilder.core.constants import HBAR
        self.assertAlmostEqual(angular_momentum_z(0), 0.0, places=40)
        self.assertAlmostEqual(angular_momentum_z(1), HBAR, places=45)
        self.assertAlmostEqual(angular_momentum_z(-1), -HBAR, places=45)

    # --- expectation values ---
    def test_expectation_r_hydrogen_1s(self):
        from molbuilder.atomic.wavefunctions import expectation_r
        from molbuilder.core.constants import BOHR_RADIUS_M
        # <r> for 1s hydrogen = 1.5 * a_0
        r_mean = expectation_r(1, 0, Z=1)
        self.assertAlmostEqual(r_mean, 1.5 * BOHR_RADIUS_M, places=20)

    def test_most_probable_radius_1s(self):
        from molbuilder.atomic.wavefunctions import most_probable_radius
        from molbuilder.core.constants import BOHR_RADIUS_M
        # For 1s (n=1, l=0 => l=n-1), r_mp = a_0
        r_mp = most_probable_radius(1, 0, Z=1)
        self.assertAlmostEqual(r_mp, BOHR_RADIUS_M, places=15)

    # --- orbital_label ---
    def test_orbital_label_1s(self):
        from molbuilder.atomic.wavefunctions import orbital_label
        self.assertEqual(orbital_label(1, 0), "1s")

    def test_orbital_label_3d(self):
        from molbuilder.atomic.wavefunctions import orbital_label
        self.assertEqual(orbital_label(3, 2), "3d")

    def test_orbital_label_with_m(self):
        from molbuilder.atomic.wavefunctions import orbital_label
        label = orbital_label(2, 1, m=1)
        self.assertIn("p_x", label)

    # --- wavefunction ---
    def test_wavefunction_1s_complex(self):
        from molbuilder.atomic.wavefunctions import wavefunction
        from molbuilder.core.constants import BOHR_RADIUS_M
        psi = wavefunction(1, 0, 0, BOHR_RADIUS_M, 0.5, 0.5, Z=1)
        self.assertNotAlmostEqual(abs(complex(psi)), 0.0)

    def test_wavefunction_real_1s(self):
        from molbuilder.atomic.wavefunctions import wavefunction_real
        from molbuilder.core.constants import BOHR_RADIUS_M
        psi = wavefunction_real(1, 0, 0, BOHR_RADIUS_M, 0.5, 0.5, Z=1)
        self.assertNotAlmostEqual(float(psi), 0.0)


if __name__ == "__main__":
    unittest.main()
