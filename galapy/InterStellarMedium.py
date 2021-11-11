# External Imports
import numpy
# from abc import ABC, abstractmethod

# Internal imports
from .ISM_core import CMC, CDD, total_attenuation

ism_tunables = {
    'mc' : [ 'f_MC', 'norm_MC', 'N_MC', 'R_MC', 'Zgas', 'tau_esc', 'Mgas' ],
    'dd' : [ 'f_MC', 'norm_DD', 'Mdust', 'Rdust', 'f_PAH' ],
}

mc_dict = {
    'f_MC'    : 0.5,
    'norm_MC' : 1.e+02,
    'N_MC'    : 1.e+03,
    'R_MC'    : 10.,
    'Zgas'    : 0.5,
    'tau_esc' : 1.e+07,
    'Mgas'    : 1.e+09,
}

dd_dict = {
    'f_MC'    : 0.5,
    'norm_DD' : 1.e+00,
    'Mdust'   : 1.e+07,
    'Rdust'   : 1.e+03,
    'f_PAH'   : 0.2
}

def gen_params_dict ( phase, **kwargs ) :
    if phase == 'mc' :
        return mc_dict
    if phase == 'dd' :
        return dd_dict
    av_phases = ''
    for k in ism_tunables.keys() :
        av_phases += f'"{k}" '
    raise ValueError( f'Phase "{phase}" chosen not available. '
                      f'Available phases are: {av_phases}')
    
class ismPhase (): 

    def __init__ ( self, phase, builder, T = None, **kwargs ) :
        
        self.phase  = phase
        self.params = gen_params_dict( self.phase, **kwargs )
        self.core   = builder()
        self.set_parameters( **kwargs )
        if T :
            self.set_temperature( T )
                
        # steal docstrings from C-core:
        self.set_temperature.__func__.__doc__ = self.core.set_temperature.__doc__
        self.temperature.__func__.__doc__     = self.core.temperature.__doc__
        self.emission.__func__.__doc__        = self.core.emission.__doc__
        self.attenuation.__func__.__doc__     = self.core.attenuation.__doc__
        self.extinction.__func__.__doc__      = self.core.extinction.__doc__
        self.A_V.__func__.__doc__             = self.core.A_V.__doc__

    def set_parameters ( self, **kwargs ) :
        """ Function for ...
        """
        self.params.update( kwargs )
        self.core.set_params( numpy.asarray( [
            self.params[k]
            for k in ism_tunables[self.phase]
        ], dtype = float ) )
    
    def set_temperature ( self, T ) :
        self.T = T
        self.core.set_temperature( self.T )
        return;

    def temperature ( self, Etot ) :
        self.T = self.core.temperature( Etot )
        return self.T

    def emission ( self, ll, T = None ) :
        if T :
            self.set_temperature( T )
        return self.core.emission( ll )

    def attenuation ( self, ll ) :
        return self.core.attenuation( ll )

    def extinction ( self, ll ) :
        return self.core.extinction( ll )

    def A_V ( self ) :
        return self.core.A_V()

class MC ( ismPhase ) :

    def __init__ ( self, T = None, **kwargs ) :
        super().__init__( 'mc', CMC, T )

    def eta ( self, tt ) :
        return self.core.eta( tt )

class DD ( ismPhase ) :

    def __init__ ( self, T = None, **kwargs ) :
        super().__init__( 'dd', CDD, T )

class ISM () :

    def __init__ ( self ) :
        self.mc = MC()
        self.dd = DD()

    def set_parameters ( self ) :
        pass

    def total_attenuation ( self ) :
        pass


        
