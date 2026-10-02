
# External imports
import numpy as np
import pytest


from galapy import StarFormationHistory as gpsfh
from galapy import CompositeStellarPopulation as gpcsp
from galapy import InterStellarMedium as gpism
from galapy import ActiveGalacticNucleus as gpagn
from galapy import XRayBinaries as gpxrb
from galapy import NebularFreeFree as gpnff
from galapy import PhotometricSystem as gppms
from galapy import Synchrotron as gpsyn
from galapy import Cosmology as gpcsm
from galapy.Galaxy import GXY
from galapy.internal.utils import trap_int

# ===========================================================================
# AGN fraction (GXY.get_fAGN)
# ===========================================================================

_GXY_KW = dict( age = 1.e+9, redshift = 0.5 )

def _stellar_luminosity ( gxy ) :
    """Unattenuated stellar luminosity in L_sun, integrated as in get_fAGN."""
    wn0 = gxy.components['stellar'] > 0
    return trap_int( gxy.wl()[wn0], gxy.components['stellar'][wn0] )

def test_gxy_fAGN_without_agn () :
    gxy = GXY( **_GXY_KW )
    assert gxy.get_fAGN() == 0.0

def test_gxy_fAGN_fritz2006_is_parameter () :
    gxy = GXY( do_AGN = True, agn = { 'model' : 'Fritz2006', 'fAGN' : 0.2 },
               **_GXY_KW )
    assert gxy.get_fAGN() == pytest.approx( 0.2 )

@pytest.mark.parametrize( 'ratio', [ 1.e-2, 1., 9. ] )
def test_gxy_fAGN_panchromatic ( ratio ) :
    # Lbol and the stellar emission are both in solar units: an AGN with
    # Lbol = ratio * L_stellar must give fAGN = ratio / (1 + ratio).
    # Guards against unit mismatches (e.g. a spurious erg/s -> L_sun factor
    # on Lbol, which would drive fAGN to ~0).
    probe = GXY( do_AGN = True, agn = { 'model' : 'Panchromatic' }, **_GXY_KW )
    Lstar = _stellar_luminosity( probe )
    gxy = GXY( do_AGN = True,
               agn = { 'model' : 'Panchromatic', 'Lbol' : ratio * Lstar },
               **_GXY_KW )
    assert gxy.get_fAGN() == pytest.approx( ratio / ( 1. + ratio ), rel = 1.e-6 )
