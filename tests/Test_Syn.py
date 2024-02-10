
#------------------------------------------------------------------------------#
#                     Test of galapy.Synchrotron module.
#------------------------------------------------------------------------------#

"""

Testing the Class SYN that implements a generic parameterized synchrotron 
emission.


    Parameters
    ----------
    ll: float or array-like
    Wavelength grid where the synchrotron emission is computed.


"""


# External imports
import numpy as np
import pytest

# Internal imports
import galapy
from galapy import Synchrotron as gpsyn



def test_syn_build_params () :

    """
    Check the function that sets the parameter values different from default 
    ones.

    """
    
    assert gpsyn.syn_build_params()  == pytest.approx({
        'alpha_syn' : 0.75, 'nu_self_syn' : 0.2
    })
    assert gpsyn.syn_build_params( alpha_syn = 0.5 , nu_self_syn = 0.1 ) == pytest.approx({
        'alpha_syn' : 0.5, 'nu_self_syn' : 0.1
    })
    
#------------------------------------------------------------------------------#

@pytest.fixture
def syn () :
        
    """
    SYN class is initialize with a pytest fixture to avoid repetitions.
    Required parameters: wavelength grid
  
    """

    ll = np.logspace( 6., 10., 5 )
    syn =  gpsyn.SYN( ll )
    return syn

#------------------------------------------------------------------------------#

def test_syn_init ( syn ):
   
    """
    Check the correct initialiation of the class.
    
    """
 
    assert isinstance( syn, galapy.Synchrotron.SYN )

#------------------------------------------------------------------------------#

def test_syn_set_params ( syn ):
    
    """
    Set the class parameters.
    
    """

    syn.set_parameters( alpha_syn = 0.5 )
    assert syn.params == pytest.approx({ 'alpha_syn' : 0.5, 'nu_self_syn' : 0.2 })

#------------------------------------------------------------------------------#

def test_syn_opt_depth (syn):
    
    """
    Check the value of the optical depth for different wavelengths.
    
    """

    assert np.all( syn.opt_depth() == pytest.approx( [ 2.68336778e-14,
                                                       4.77177768e-11,
                                                       8.48555400e-08,
                                                       1.50896860e-04,
                                                       2.68336778e-01 ] ) )

#------------------------------------------------------------------------------#

def test_syn_energy (syn):
    
    """
    Check the value of the normalized total energy radiated at different
    wavelengths.
    
    """

    assert np.all( syn.energy(1) == pytest.approx([1.86610704e-04,
                                                   2.84911159e-03,
                                                   3.50902904e-02,
                                                   3.16390758e-01,
                                                   1.92867026e+00] ))

#------------------------------------------------------------------------------#

def test_syn_emission (syn):
    
    """
    Check the value of the normalized synchrotron emission at different 
    wavelengths.
    
    """

    assert np.all( syn.emission(1) == pytest.approx([5.59444816e+02,
                                                     8.54142168e+01,
                                                     1.05198044e+01,
                                                     9.48515631e-01,
                                                     5.78200798e-02]))
    
#------------------------------------------------------------------------------#

def test_snsyn_init():

    ll = np.logspace(6,10,5)
    snsyn=gpsyn.SNSYN(ll)
    assert isinstance( snsyn, galapy.Synchrotron.SNSYN )

#------------------------------------------------------------------------------#

def test_snsyn_emission():

    ll = np.logspace(6,13,5)
    assert np.all( gpsyn.SNSYN(ll).emission() == pytest.approx(0))

