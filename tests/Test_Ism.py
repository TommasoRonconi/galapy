
#------------------------------------------------------------------------------#
#                  Test of galapy.InterStellarMedium module.
#------------------------------------------------------------------------------#

""" 
Testing the "galapy.InterStellarMedium" module which defines classes and 
functions for modelling the gas and dust components impact on the total 
spectrum of a galaxy .

In galapy, ISM is modelled as a medium constituted of two phases, 
a Molecular Cloud (MC) component, embedding star-forming regions, 
and a Diffuse Dust (DD) component which is roughly considered as a spherical 
region centred at the centre of the galaxy .

"""

# External imports
import numpy as np
import pytest

# Internal imports
import galapy
from galapy import InterStellarMedium as gpism



def test_ism_parameters ( ) :

    """
    Available phases for this module are only 'mc' and 'dd' .
    
    """

    phase = 'new'
    with pytest.raises( ValueError, match = f'Phase "{phase}" chosen not '
                        'available. Available phases are: "mc" "dd"' ) :
        gpism.ism_build_params( phase = phase )

#------------------------------------------------------------------------------#

@pytest.fixture
def ism ():
    
    ism = gpism.ISM()
    return ism
 

def test_ism_init ( ism ) :

    """
    Testing initialization of class "ISM" .
    
    """
    
    assert isinstance( ism, galapy.InterStellarMedium.ISM )
    
        
################################################################################
################################## DIFFUSE DUST ################################
################################################################################

def test_ism_dd_parameters ( ism ) :

    """
    Test the default parameters set for the DD component .
    
    """
    
    assert gpism.ism_build_params( phase = 'dd' ) == { 'f_MC': 0.5,
                                                      'norm_DD': 1.0,
                                                      'Mdust': 10000000.0,
                                                      'Rdust': 1000.0,
                                                       'f_PAH': 0.2,
                                                       'dDDlow': 0.7,
                                                       'dDDupp': 2.0 }
    
#------------------------------------------------------------------------------#

def test_ism_dd_temperature ( ism ) :

    """
    Test DD component temperature with a certain total emitted power .
    
    """
    
    E_tot = 4.e+43
    assert ism.dd.temperature( E_tot ) == pytest.approx( 37.1029777891508 )
    
#------------------------------------------------------------------------------#

def test_ism_dd_emission ( ism ) :

    """
    Test DD emission at T = 0 and T = 38 K .
    
    """

    ll = np.logspace( 3, 8, 5 )
    ism.dd.set_temperature( 0. )
    assert all ( [ ism.dd.emission( l ) == pytest.approx( 0. ) for l in ll ] )
    ism.dd.set_temperature( 38. )
    assert np.all ( ism.dd.emission( ll ) == pytest.approx( [ 0.,
                                                              0.,
                                                              1.67584039e+03,
                                                              2.20224278e+00,
                                                              9.75030567e-08 ] )
                    )

#------------------------------------------------------------------------------#

def test_ism_dd_attenuation ( ism ) :

    """
    Test the DD component attenuation factor at T = 38 K .
    
    """
    
    ll = np.logspace( 3, 8, 5 )
    ism.dd.set_temperature( 38. )
    assert np.all ( ism.dd.attenuation( ll ) == pytest.approx( [ 0.73803851,
                                                                 0.96030249,
                                                                 0.99461287,
                                                                 0.9999237,
                                                                 0.99999976 ] )
    )
    
#------------------------------------------------------------------------------#

def test_ism_dd_extinction ( ism ) :

    """
    Test the extinction for the DD component at T = 38 K .
    
    """

    ll = np.logspace( 3, 8, 5 )
    ism.dd.set_temperature( 38. )
    assert np.all ( ism.dd.extinction( ll ) == pytest.approx( [ 3.29802447e-01,
                                                                4.39798647e-02,
                                                                5.86480941e-03,
                                                                8.28426327e-05,
                                                                2.61971407e-07 ]
    ) )

#------------------------------------------------------------------------------#

def test_ism_dd_A_V ( ism ) :

    """
    Test the DD extiction value in the visible band at T = 38 K .

    """
    
    ism.dd.set_temperature( 38. )
    assert ism.dd.A_V() == 0.1


    
################################################################################
################################ MOLECULAR CLOUD ###############################
################################################################################


def test_ism_mc_parameters ( ism ) :

    """
    Test the default parameters set for the MC component .
    
    """
    
    assert gpism.ism_build_params( phase = 'mc' ) == { 'f_MC': 0.5,
                                                      'norm_MC': 100.0,
                                                      'N_MC': 1000.0,
                                                      'R_MC': 10.0,
                                                      'Zgas': 0.01,
                                                      'tau_esc': 10000000.0,
                                                       'Mgas': 1000000000.0,
                                                       'dMClow': 1.3,
                                                       'dMCupp': 1.6 }
    
#------------------------------------------------------------------------------#

def test_ism_mc_temperature ( ism ) :

    """
    Test MC component temperature with a certain total emitted power .
    
    """
    
    ism.set_parameters( Zgas = 0.0134 )
    E_tot = 4.e+43
    assert ism.mc.temperature( E_tot ) == pytest.approx( 27.381094065189018 )
    
#------------------------------------------------------------------------------#

def test_ism_mc_emission ( ism ) :

    """
    Test MC emission at T = 0 and T = 30 K .
    
    """

    ll = np.logspace( 3, 8, 5 )
    ism.mc.set_temperature( 0. )
    assert all ( [ ism.mc.emission( l ) == pytest.approx( 0. ) for l in ll ] )
    ism.mc.set_temperature( 30. )
    assert np.all ( ism.mc.emission( ll ) == pytest.approx( [ 0.,
                                                              0.,
                                                              4.65817437e+2,
                                                              1.31598377e+1,
                                                              2.03449423e-6 ] )
                    )

#------------------------------------------------------------------------------#

def test_ism_mc_attenuation ( ism ) :

    """
    Test the MC component attenuation factor at T = 30 K .
    
    """

    ll = np.logspace( 3, 8, 5 )
    ism.mc.set_temperature( 30. )
    assert np.all ( ism.mc.attenuation( ll ) == pytest.approx( [0.,
                                                                4.89009895e-09,
                                                                6.35217767e-01,
                                                                9.93610585e-01,
                                                                9.99935903e-01]
    ) )
#------------------------------------------------------------------------------#

def test_ism_mc_extinction (ism) :

    """
    Test the extinction for the MC component at T = 30 K .
    
    """
    
    ll = np.logspace( 3, 8, 5 )
    ism.mc.set_temperature( 30. )
    assert np.all ( ism.mc.extinction( ll ) == pytest.approx( [8.76146302e+02,
                                                               2.07767059e+01,
                                                               4.92693408e-01,
                                                               6.95947684e-03,
                                                               6.95947636e-05]
    ) )
#------------------------------------------------------------------------------#

def test_ism_mc_A_V (ism) :
    
    """
    Test the MC extiction value in the visible band at T = 30 K .

    """
    
    ism.mc.set_temperature( 30. )
    assert ism.mc.A_V() == pytest.approx( 95.52238805014926 )

#------------------------------------------------------------------------------#

def test_ism_mc_eta (ism) :

    """
    Check the eta behaviour at different time .

    """
    
    assert ism.mc.eta( 1.e+6 ) == 1.
    assert ism.mc.eta( 1.5e+7 ) == 0.5
    assert ism.mc.eta( 1.e+8 ) == 0.

################################################################################
############################## TOTAL ATTENUATION ###############################
################################################################################


"""
Total attenuation, as it depends on the time spent by the photons within the 
molecular cloud through the parameter tau_esc, is a function of time.
The function total_attenuation of the ISM class returns a tuple with 2 matrices, 
of shape (len(ll), len(tt)) where ll is the wavelenght grid and tt the time 
grid, containing the contribute to attenuation from MC-only and the global 
attenuation itself.

"""

def test_ism_total_attenuation (ism) :
    
    """
    Test the total attenuation at 2 different times for 2 certain wavelength.

    """

    ll=np.logspace( 3, 8, 2 )
    t_esc = ism.mc.params[ 'tau_esc' ]
    tunit = np.logspace( np.log10( 1 ), np.log10( 4 ), 2 )
    tt = tunit * t_esc
    attTot = ism.total_attenuation( ll, tt )[ 1 ]
    assert np.all ( attTot == pytest.approx( [ 0.        ,
                                               0.73803851,
                                               0.99993566,
                                               0.99999976 ] ) )


#------------------------------------------------------------------------------#

def test_ism_attenuation_f_MC_1 (ism) :

    """
    Test the total attenuation if the only contribution is given by the MC .
    The DD component attenuation factor must be equal to 1 .
    The total attenuation must be equal to the MC attenuation at any time .
    Furthermore, the total attenuation at time 0. must give the same result 
    of ism.mc.attenuation .

    """
    
   
    ll = np.logspace( 3, 8, 2 )
    ism.set_parameters( f_MC = 1. )
    tunit = np.logspace( np.log10( 1 ), np.log10( 10 ), 1 )
    t_esc = ism.mc.params[ 'tau_esc' ]
    tt = tunit * t_esc
    attTotMC, attTot = ism.total_attenuation( ll, tt )
    assert np.all ( ism.dd.attenuation( ll ) == pytest.approx( [ 1., 1.] ) )
    assert np.all ( ism.mc.attenuation( ll ) == pytest.approx( attTot ) )
    assert np.all( attTotMC == attTot )

#------------------------------------------------------------------------------#


def test_ism_attenuation_f_MC_0 (ism) :

    """
    Test the total attenuation if the only contribution is given by the DD .
    The MC component attenuation factor must be equal to 1 .
    Furthermore, the total attenuation at time 0. must give the same result 
    of ism.dd.attenuation .

    """
    ll = np.logspace( 3, 8, 2 )
    ism.set_parameters( f_MC = 0. )
    tunit = np.logspace( np.log10( 1 ), np.log10( 10 ), 1 )
    t_esc = ism.mc.params[ 'tau_esc' ]
    tt = tunit * t_esc
    attTotMC, attTot = ism.total_attenuation( ll, tt )
    assert np.all ( ism.mc.attenuation( ll ) == pytest.approx( [ 1., 1. ] ) )
    assert np.all ( ism.dd.attenuation( ll ) == pytest.approx( attTot ) )
