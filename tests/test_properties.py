"""Tests for Lipinski Rule-of-5 molecular property calculations."""

import pytest

from molbuilder.smiles.parser import parse
from molbuilder.molecule.properties import (
    LipinskiProperties,
    lipinski_properties,
    hydrogen_bond_donors,
    hydrogen_bond_acceptors,
    rotatable_bond_count,
    crippen_logp,
    topological_polar_surface_area,
    heavy_atom_count,
)


# =====================================================================
#  Helper
# =====================================================================

def mol(smiles: str):
    """Parse a SMILES string into a Molecule."""
    return parse(smiles)


# =====================================================================
#  hydrogen_bond_donors
# =====================================================================

class TestHydrogenBondDonors:
    def test_ethanol(self):
        # CCO -> one OH
        assert hydrogen_bond_donors(mol("CCO")) == 1

    def test_ethylamine(self):
        # CCN -> one NH2 (counts as 1 donor atom)
        assert hydrogen_bond_donors(mol("CCN")) == 1

    def test_triethylamine(self):
        # N(CC)(CC)CC -> no N-H
        assert hydrogen_bond_donors(mol("N(CC)(CC)CC")) == 0

    def test_benzene(self):
        assert hydrogen_bond_donors(mol("c1ccccc1")) == 0

    def test_water(self):
        assert hydrogen_bond_donors(mol("O")) == 1

    def test_acetic_acid(self):
        # CC(=O)O -> one OH
        assert hydrogen_bond_donors(mol("CC(=O)O")) == 1

    def test_ethanediol(self):
        # OCCO -> two OH
        assert hydrogen_bond_donors(mol("OCCO")) == 2


# =====================================================================
#  hydrogen_bond_acceptors
# =====================================================================

class TestHydrogenBondAcceptors:
    def test_ethanol(self):
        # CCO -> 1 O
        assert hydrogen_bond_acceptors(mol("CCO")) == 1

    def test_acetic_acid(self):
        # CC(=O)O -> 2 O atoms
        assert hydrogen_bond_acceptors(mol("CC(=O)O")) == 2

    def test_aniline(self):
        # c1ccccc1N -> 1 N
        assert hydrogen_bond_acceptors(mol("c1ccccc1N")) == 1

    def test_benzene(self):
        assert hydrogen_bond_acceptors(mol("c1ccccc1")) == 0

    def test_dimethylether(self):
        # COC -> 1 O
        assert hydrogen_bond_acceptors(mol("COC")) == 1


# =====================================================================
#  rotatable_bond_count
# =====================================================================

class TestRotatableBonds:
    def test_ethane(self):
        # C-C only, both terminal -> 0 rotatable
        assert rotatable_bond_count(mol("CC")) == 0

    def test_butane(self):
        # C-C-C-C -> middle C-C is rotatable
        assert rotatable_bond_count(mol("CCCC")) == 1

    def test_pentane(self):
        # CCCCC -> 2 rotatable (middle two C-C)
        assert rotatable_bond_count(mol("CCCCC")) == 2

    def test_cyclohexane(self):
        # Ring bonds are not rotatable
        assert rotatable_bond_count(mol("C1CCCCC1")) == 0

    def test_ethene(self):
        # Double bond -> not rotatable
        assert rotatable_bond_count(mol("C=C")) == 0

    def test_methane(self):
        assert rotatable_bond_count(mol("C")) == 0


# =====================================================================
#  crippen_logp
# =====================================================================

class TestCrippenLogP:
    def test_methane(self):
        # CH4 ~ 1.1 (literature ~1.09)
        lp = crippen_logp(mol("C"))
        assert 0.4 < lp < 1.6, f"methane logP={lp}"

    def test_ethanol(self):
        # Ethanol ~ -0.31
        lp = crippen_logp(mol("CCO"))
        assert -1.0 < lp < 0.5, f"ethanol logP={lp}"

    def test_benzene(self):
        # Benzene ~ 2.13
        lp = crippen_logp(mol("c1ccccc1"))
        assert 1.0 < lp < 3.0, f"benzene logP={lp}"

    def test_octanol(self):
        # 1-octanol ~ 3.00
        lp = crippen_logp(mol("CCCCCCCCO"))
        assert 2.0 < lp < 4.0, f"octanol logP={lp}"

    def test_water(self):
        # Water logP < 0
        lp = crippen_logp(mol("O"))
        assert lp < 0.5, f"water logP={lp}"

    def test_hexane(self):
        # Hexane ~ 3.9
        lp = crippen_logp(mol("CCCCCC"))
        assert 2.5 < lp < 5.0, f"hexane logP={lp}"


# =====================================================================
#  topological_polar_surface_area
# =====================================================================

class TestTPSA:
    def test_benzene(self):
        # No polar atoms -> TPSA = 0
        assert topological_polar_surface_area(mol("c1ccccc1")) == 0.0

    def test_hexane(self):
        assert topological_polar_surface_area(mol("CCCCCC")) == 0.0

    def test_ethanol(self):
        # Ethanol: one -OH -> ~20.23
        tpsa = topological_polar_surface_area(mol("CCO"))
        assert 15.0 < tpsa < 25.0, f"ethanol TPSA={tpsa}"

    def test_acetic_acid(self):
        # CC(=O)O -> =O (~17.07) + -OH (~20.23) = ~37.3
        tpsa = topological_polar_surface_area(mol("CC(=O)O"))
        assert 30.0 < tpsa < 45.0, f"acetic acid TPSA={tpsa}"

    def test_methylamine(self):
        # CN -> -NH2 -> ~26.02
        tpsa = topological_polar_surface_area(mol("CN"))
        assert 20.0 < tpsa < 32.0, f"methylamine TPSA={tpsa}"


# =====================================================================
#  heavy_atom_count
# =====================================================================

class TestHeavyAtomCount:
    def test_methane(self):
        assert heavy_atom_count(mol("C")) == 1

    def test_ethanol(self):
        # CCO -> C, C, O = 3
        assert heavy_atom_count(mol("CCO")) == 3

    def test_benzene(self):
        assert heavy_atom_count(mol("c1ccccc1")) == 6

    def test_water(self):
        assert heavy_atom_count(mol("O")) == 1


# =====================================================================
#  lipinski_properties (integration)
# =====================================================================

class TestLipinskiProperties:
    def test_returns_dataclass(self):
        props = lipinski_properties(mol("CCO"))
        assert isinstance(props, LipinskiProperties)

    def test_ethanol_fields(self):
        props = lipinski_properties(mol("CCO"))
        assert props.heavy_atom_count == 3
        assert props.hba == 1
        assert props.hbd == 1
        assert props.molecular_weight > 40  # ~46
        assert props.lipinski_pass is True
        assert props.lipinski_violations == 0

    def test_aspirin(self):
        # Aspirin: CC(=O)Oc1ccccc1C(=O)O
        # MW ~180, logP ~1.2, HBD=1, HBA=4, should pass
        props = lipinski_properties(mol("CC(=O)Oc1ccccc1C(=O)O"))
        assert props.lipinski_violations == 0
        assert props.lipinski_pass is True
        assert props.hba == 4
        assert props.hbd == 1

    def test_benzene_no_violations(self):
        props = lipinski_properties(mol("c1ccccc1"))
        assert props.lipinski_violations == 0
        assert props.lipinski_pass is True
        assert props.hbd == 0
        assert props.hba == 0
        assert props.tpsa == 0.0

    def test_all_fields_present(self):
        props = lipinski_properties(mol("CCO"))
        assert hasattr(props, "molecular_weight")
        assert hasattr(props, "logp")
        assert hasattr(props, "hbd")
        assert hasattr(props, "hba")
        assert hasattr(props, "rotatable_bonds")
        assert hasattr(props, "tpsa")
        assert hasattr(props, "heavy_atom_count")
        assert hasattr(props, "lipinski_violations")
        assert hasattr(props, "lipinski_pass")

    def test_large_molecule_violations(self):
        # Build a large-ish molecule that should violate MW
        # 30-carbon chain: MW ~422, logP very high (>5), so 1-2 violations
        big = "C" * 30
        props = lipinski_properties(mol(big))
        assert props.heavy_atom_count == 30
        assert props.logp > 5  # should violate logP
        assert props.lipinski_violations >= 1

    def test_water(self):
        props = lipinski_properties(mol("O"))
        assert props.heavy_atom_count == 1
        assert props.hbd == 1
        assert props.hba == 1
        assert props.lipinski_pass is True

    def test_rotatable_bonds_in_properties(self):
        props = lipinski_properties(mol("CCCC"))
        assert props.rotatable_bonds == 1

    def test_tpsa_in_properties(self):
        props = lipinski_properties(mol("CCO"))
        assert props.tpsa > 0
