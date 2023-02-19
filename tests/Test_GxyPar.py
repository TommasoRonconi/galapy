
# External imports
import numpy as np
import pytest
import warnings

# Internal imports
import galapy
from galapy.Galaxy import GXY
from galapy import Handlers as gpgxp
from galapy.internal.utils import set_nested, unwrap_items



@pytest.fixture
def gxp () :

    gxy_model = GXY()
    sample_params = {}
    gxp = gpgxp.GXYParameters( gxy_model, sample_params  )
    return gxp

def test_gxp_init ( gxp ) :

    assert isinstance ( gxp, galapy.Handlers.GXYParameters ) 

def test_gxp_default_params ( gxp ) :

    assert gxp.parameters == { 'galaxy.age': 1000000.0,
                               'galaxy.redshift': 0.0,
                               'galaxy.sfh.tau_quench': 2e+20,
                               'galaxy.sfh.model': 'insitu',
                               'galaxy.sfh.psi_max': 100.0,
                               'galaxy.sfh.tau_star': 300000000.0,
                               'galaxy.ism.f_MC': 0.5,
                               'galaxy.ism.norm_MC': 100.0,
                               'galaxy.ism.N_MC': 1000.0,
                               'galaxy.ism.R_MC': 10.0,
                               'galaxy.ism.Zgas': 6.663073213343221e-05,
                               'galaxy.ism.tau_esc': 10000000.0,
                               'galaxy.ism.Mgas': 96758967.01542452,
                               'galaxy.ism.dMClow': 1.3,
                               'galaxy.ism.dMCupp': 1.6,
                               'galaxy.ism.norm_DD': 1.0,
                               'galaxy.ism.Mdust': 239.11118605751645,
                               'galaxy.ism.Rdust': 1000.0,
                               'galaxy.ism.f_PAH': 0.2,
                               'galaxy.ism.dDDlow': 0.7,
                               'galaxy.ism.dDDupp': 2.0 }

def test_gxp_invalid_model () :

    
    gxy_model = 'not_GXY_object'
    sample_params = {}

    with pytest.raises( RuntimeError,
                        match = f'The provided model ``{gxy_model}`` is not valid. Valid models are instances of class galapy.internal.abc.Model' ):
        gxp = gpgxp.GXYParameters( gxy_model, sample_params )
        

def test_gxp_invalid_value ( recwarn ) :

    value = 'string'
    gxp =  gpgxp.GXYParameters( GXY(), { 'age' : value } )
    assert gxp.parameters['galaxy.age'] == 1.e+6


def test_gxp_invalid_key () :

    key= 'wrong_key'
    with pytest.warns( UserWarning, match = f'Parameter "galaxy.{key}" is not present in the model requested and will be ignored.' ):
        gxp = gpgxp.GXYParameters( GXY(), { key : 1. } )


def test_gxp_params_fixed () :

    gxp =  gpgxp.GXYParameters( GXY(), { 'age' : 1.e+8 , 'ism.f_MC' : 0.6 } )
    assert gxp.parameters[ 'galaxy.age' ] == 1.e+8
    assert gxp.parameters[ 'galaxy.ism.f_MC' ] == 0.6 
    
def test_gxp_params_free () :

    value_age = ( [5., 9.] , True )
    value_fMC = ( [0., 1.] , False )
    gxp =  gpgxp.GXYParameters( GXY(), { 'age' : value_age, 'ism.f_MC' : value_fMC } )
    assert 1.e+5 <= gxp.parameters[ 'galaxy.age' ] <= 1.e+9
    assert 0. <= gxp.parameters[ 'galaxy.ism.f_MC' ] <= 1.
    assert np.all( gxp.par_free == [ 'galaxy.age', 'galaxy.ism.f_MC' ] )
    assert np.all( gxp.par_prior == [ [ 5., 9. ], [ 0., 1. ] ] )
    assert np.all( gxp.par_log == [ True, False ] )

def test_gxp_return_nested_err () :

     value_fMC = ( [0., 1.] , False )
     par = ['galaxy.age', 'galaxy.ism.f_MC']
     gxp =  gpgxp.GXYParameters( GXY(), { 'age' : 1.e+8 , 'ism.f_MC' : value_fMC } )
     with pytest.raises( RuntimeError, match = f'Provided {len(par)} but there are exactly '
                                f'{len(gxp.par_free)} free parameters to set.' ):
         
        gxp.return_nested( par )
        
def test_gxp_return_nested () :

     value_fMC = ( [0., 1.] , False )
     par = np.array([0.1])
     gxp =  gpgxp.GXYParameters( GXY(), { 'age' : 1.e+8 , 'ism.f_MC' : value_fMC } )
     print( gxp.return_nested( par ))
        
