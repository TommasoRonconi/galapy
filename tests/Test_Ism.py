
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




def test_ism_init () :

    """
    Testing initialization of class "ISM" .
    
    """

    ism = gpism.ISM( )
    assert isinstance( ism, galapy.InterStellarMedium.ISM )
    
#------------------------------------------------------------------------------#

def test_ism_parameters () :

    """
    Available phases for this module are only 'mc' and 'dd' .
    
    """

    ism = gpism.ISM( )
    phase = 'new'
    with pytest.raises( ValueError, match = f'Phase "{phase}" chosen not '
                        'available. Available phases are: "mc" "dd"' ) :
        gpism.gen_params_dict( phase = phase )

        
################################################################################
################################## DIFFUSE DUST ################################
################################################################################

def test_ism_dd_parameters () :

    """
    Test the default parameters set for the DD component .
    
    """

    ism = gpism.ISM( )
    ism.set_parameters( )
    assert gpism.gen_params_dict( phase = 'dd' ) == { 'f_MC': 0.5,
                                                      'norm_DD': 1.0,
                                                      'Mdust': 10000000.0,
                                                      'Rdust': 1000.0,
                                                      'f_PAH': 0.2 }
    
#------------------------------------------------------------------------------#

def test_ism_dd_temperature () :

    """
    Test DD component temperature with a certain total emitted power .
    
    """
    ism = gpism.ISM( )
    ism.set_parameters( )
    E_tot = 4.e+43
    assert ism.dd.temperature( E_tot ) == pytest.approx( 37.1029777891508 )
    
#------------------------------------------------------------------------------#

def test_ism_dd_emission () :

    """
    Test DD emission at T = 0 and T = 38 K .
    
    """

    ll = np.logspace( 5, 8, 6 )
    ism = gpism.ISM( )
    ism.set_parameters( )
    ism.dd.set_temperature( 0. )
    assert all ( [ ism.dd.emission( l ) == pytest.approx( 0. ) for l in ll ] )
    ism.dd.set_temperature( 38. )
    assert np.all ( ism.dd.emission( ll ) == pytest.approx( [ 6.74123441e-06,
                                                              5.29600438e+03,
                                                              1.51184833e+03,
                                                              1.14933054e+00,
                                                              3.66540302e-04,
                                                              9.75030567e-08 ] )
                    )

#------------------------------------------------------------------------------#

def test_ism_dd_attenuation () :

    """
    Test the DD component attenuation factor at T = 38 K .
    
    """
    
    ll = np.logspace( 2, 6, 6 )
    ism = gpism.ISM( )
    ism.set_parameters( )
    ism.dd.set_temperature( 38. )
    assert np.all ( ism.dd.attenuation( ll ) == pytest.approx( [ 0.21818702,
                                                                 0.6575034,
                                                                 0.89093299,
                                                                 0.96869304,
                                                                 0.99127775,
                                                                 0.99759006 ] ))
    
#------------------------------------------------------------------------------#

def test_ism_dd_extinction () :

    """
    Test the extinction for the DD component at T = 38 K .
    
    """

    ll = np.logspace( 2, 6, 6 )
    ism = gpism.ISM( )
    ism.set_parameters( )
    ism.dd.set_temperature( 38. )
    assert np.all ( ism.dd.extinction( ll ) == pytest.approx( [1.65292772,
                                                               0.45525411,
                                                               0.1253874,
                                                               0.03453456,
                                                               0.00951161,
                                                               0.002619714 ] ) )

#------------------------------------------------------------------------------#

def test_ism_dd_A_V () :

    """
    Test the DD extiction value in the visible band at T = 38 K .

    """
    
    ism = gpism.ISM( )
    ism.dd.set_temperature( 38. )
    assert ism.dd.A_V() == 0.1


    
################################################################################
################################ MOLECULAR CLOUD ###############################
################################################################################


def test_ism_mc_parameters () :

    """
    Test the default parameters set for the MC component .
    
    """

    ism = gpism.ISM( )
    assert gpism.gen_params_dict( phase = 'mc' ) == { 'f_MC': 0.5,
                                                      'norm_MC': 100.0,
                                                      'N_MC': 1000.0,
                                                      'R_MC': 10.0,
                                                      'Zgas': 0.5,
                                                      'tau_esc': 10000000.0,
                                                      'Mgas': 1000000000.0 }
    
#------------------------------------------------------------------------------#

def test_ism_mc_temperature () :

    """
    Test MC component temperature with a certain total emitted power .
    
    """
    ism = gpism.ISM( )
    ism.set_parameters( )
    E_tot = 4.e+43
    assert ism.mc.temperature( E_tot ) == pytest.approx( 30.918524140490067 )
    
#------------------------------------------------------------------------------#

def test_ism_mc_emission () :

    """
    Test MC emission at T = 0 and T = 30 K .
    
    """

    ll = np.logspace( 5, 8, 6 )
    ism = gpism.ISM( )
    ism.set_parameters( )
    ism.mc.set_temperature( 0. )
    assert all ( [ ism.mc.emission( l ) == pytest.approx( 0. ) for l in ll ] )
    ism.mc.set_temperature( 30. )
    assert np.all ( ism.mc.emission( ll ) == pytest.approx( [ 1.71645096e-09,
                                                              1.84418104e+03,
                                                              2.54795668e+03,
                                                              4.88382067e+00,
                                                              2.90265588e-03,
                                                              1.36312555e-06 ] )
                    )

#------------------------------------------------------------------------------#

def test_ism_mc_attenuation () :

    """
    Test the MC component attenuation factor at T = 30 K .
    
    """

    ll = np.logspace( 5, 7, 6 )
    ism = gpism.ISM( )
    ism.mc.set_temperature( 30. )
    assert np.all ( ism.mc.attenuation( ll ) == pytest.approx( [ 0.25715255,
                                                                 0.66356041,
                                                                 0.88350449,
                                                                 0.96794675,
                                                                 0.99256455,
                                                                 0.99829173 ] )
                    )
#------------------------------------------------------------------------------#

def test_ism_mc_extinction () :

    """
    Test the extinction for the MC component at T = 30 K .
    
    """
    ll = np.logspace( 4, 6, 6 )
    ism = gpism.ISM( )
    ism.set_parameters( )
    ism.mc.set_temperature( 30. )
    assert np.all ( ism.mc.extinction( ll ) == pytest.approx( [ 29.4205971,
                                                                8.88487867,
                                                                2.68319058,
                                                                0.81031064,
                                                                0.24470991,
                                                                0.07390121 ] ) )
#------------------------------------------------------------------------------#

def test_ism_mc_A_V () :
    
    """
    Test the MC extiction value in the visible band at T = 30 K .

    """
    ism = gpism.ISM( )
    ism.set_parameters( )
    ism.mc.set_temperature( 30. )
    assert ism.mc.A_V() == pytest.approx( 64. )

#------------------------------------------------------------------------------#

def test_ism_mc_eta () :

    """
    Check the eta behaviour at different time .
    """
    ism = gpism.ISM( )
    ism.set_parameters( )
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

def test_ism_total_attenuation () :
    
    """
    Test the total attenuation at 5 different times for a certain wavelength.

    """
    
    ism = gpism.ISM( )
    ll = np.logspace( 2, 7, 10 )
    t_esc = ism.mc.params[ 'tau_esc' ]
    tunit = np.logspace( np.log10( .8 ), np.log10( 4 ), 5 )
    tt = tunit * t_esc
    attTot = ism.total_attenuation( ll, tt )[ 1 ]
    assert np.all (attTot[3] == pytest.approx( [ 0.00000000e+00, 1.76941143e-01,
                                                 7.11134549e-01, 9.01477593e-01,
                                                 9.01477593e-01 ] ) )


#------------------------------------------------------------------------------#

def test_ism_attenuation_f_MC_1 () :

    """
    Test the total attenuation if the only contribution is given by the MC .
    The DD component attenuation factor must be equal to 1 .
    The total attenuation must be equal to the MC attenuation at any time .
    Furthermore, the total attenuation at time 0. must give the same result 
    of ism.mc.attenuation .

    """
    
    ism = gpism.ISM()
    ll = np.logspace( 5, 7, 6 )
    ism.set_parameters( f_MC = 1. )
    tunit = np.logspace( np.log10( 1 ), np.log10( 10 ), 5 )
    t_esc = ism.mc.params[ 'tau_esc' ]
    tt = tunit * t_esc
    attTotMC, attTot = ism.total_attenuation( ll, tt )
    assert all ( [ ism.dd.attenuation( l ) == 1. for l in ll ] )
    assert all ( [ ism.mc.attenuation( ll[ i ] ) == pytest.approx( attTot[ i ]
                                                                   [ 0 ] )
                   for i in range( 0, 6 ) ] )
    assert np.all( attTotMC == attTot )

#------------------------------------------------------------------------------#


def test_ism_attenuation_f_MC_0 () :

    """
    Test the total attenuation if the only contribution is given by the DD .
    The MC component attenuation factor must be equal to 1 .
    Furthermore, the total attenuation at time 0. must give the same result 
    of ism.dd.attenuation .

    """
    
    ism = gpism.ISM()
    ll = np.logspace( 2, 6, 6 )
    ism.set_parameters( f_MC = 0. )
    tunit = np.logspace( np.log10( 1 ), np.log10( 10 ), 5 )
    t_esc = ism.mc.params[ 'tau_esc' ]
    tt = tunit * t_esc
    attTotMC, attTot = ism.total_attenuation( ll, tt )
    assert all ( [ ism.mc.attenuation( l ) == 1. for l in ll ] )
    assert all ( [ ism.dd.attenuation( ll[ i ] ) == pytest.approx( attTot[ i ]
                                                                   [ 0 ] )
                   for i in range( 0, 6 ) ] )