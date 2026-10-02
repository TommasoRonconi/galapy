
#------------------------------------------------------------------------------#
#           Test of galapy.internal.constants module.
#------------------------------------------------------------------------------#

import numpy as np
import pytest

from galapy.internal.constants import (
    Lsun, sunL, clight, hP, Mpc_to_cm, Ang_to_keV,
)


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

def test_lsun_value():
    assert Lsun == pytest.approx(3.828e+33)

def test_sunl_is_inverse_lsun():
    assert sunL * Lsun == pytest.approx(1.0)

def test_clight_keys():
    assert set(clight.keys()) == {'cm/s', 'm/s', 'A/s'}

def test_clight_consistency():
    # c in A/s must equal c in cm/s * 1e8
    assert clight['A/s'] == pytest.approx(clight['cm/s'] * 1.e8)

def test_clight_mps():
    assert clight['m/s'] == pytest.approx(2.99792458e8)

def test_mpc_to_cm():
    assert Mpc_to_cm == pytest.approx(3.086e+24)

def test_hp_keys():
    assert set(hP.keys()) == {'eV/Hz', 'erg*s'}

def test_hp_eVHz():
    assert hP['eV/Hz'] == pytest.approx(4.1357e-15)


# ---------------------------------------------------------------------------
# Ang_to_keV
# ---------------------------------------------------------------------------

def test_ang_to_keV_1ang():
    # 1 Angstrom ~ 12.4 keV
    result = Ang_to_keV(1.0)
    assert result == pytest.approx(12.398516685505998, rel=1e-6)

def test_ang_to_keV_10ang():
    assert Ang_to_keV(10.0) == pytest.approx(1.2398516685505998, rel=1e-6)

def test_ang_to_keV_inverse_scaling():
    # doubling wavelength halves the energy
    assert Ang_to_keV(2.0) == pytest.approx(Ang_to_keV(1.0) / 2.0, rel=1e-10)

def test_ang_to_keV_array_monotone():
    lam = np.array([1., 10., 100., 1000.])
    result = Ang_to_keV(lam)
    assert result.shape == lam.shape
    assert np.all(np.diff(result) < 0)

def test_ang_to_keV_array_values():
    lam = np.array([1., 10., 100.])
    expected = np.array([12.398516685505998, 1.2398516685505998, 0.12398516685505998])
    np.testing.assert_allclose(Ang_to_keV(lam), expected, rtol=1e-6)
