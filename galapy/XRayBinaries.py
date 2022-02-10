""" The XRayBinaries module implements the contribution to the X-band emission from stellar binary systems.
"""

# External imports
from abc import ABC, abstractmethod
import numpy

# Internal imports
from galapy.internal.utils import trap_int, powerlaw_exp_cutoff, poly_N
from galapy.internal.constants import Ang_to_keV, Lsun, sunL
from galapy.internal.interp import lin_interp

xrb_tunables = {
    # High Mass Stars
    'hm' : [ 'psi', 'Zstar' ],
    
    # Low Mass Stars
    'lm' : [ 'Mstar', 'age' ],
}
""" Dictionary with the tunable parameters of the X-ray contribution from binary star systems

#. :code:`hm` : **High Mass Star**
#. :code:`lm` : **Low Mass Stars**

"""

def xrb_build_params ( startype, **kwargs ) :
    """ Builds the parameters dictionary for given star type.
    """
    if startype not in xrb_tunables.keys() :
        av_types = ''
        for k in xrb_tunables.keys() :
            av_types += f'"{k}" '
        raise ValueError( f'Star type "{startype}" chosen not available. '
                          f'Available phases are: {av_types}')
    
    if startype == 'hm' :
        out = {
            'psi'   : 1.,
            'Zstar' : 0.02,
        }
        
    if startype == 'lm' :
        out = {
            'Mstar' : 1.e+10,
            'age'   : 1.e+8,
        }

    for k in set( out.keys() ).intersection(kwargs.keys()) :
        out[k] = kwargs[k]
    return out

class _XRBase ( ABC ) :

    @property
    @abstractmethod
    def coeff ( self ) :
        pass

    def __init__ ( self, startype, gamma, Ecut,
                   lmin = 1., lmax = 1.e+10, **kwargs ) :
        self.startype = startype
        self.params   = xrb_build_params( self.startype, **kwargs )
        self.lmin     = lmin
        self.lmax     = lmax
        self.NormStarType = None

        ###########################################################
        # Compute X-template

        # generate wavelenght grid
        ll = numpy.logspace( numpy.log10( self.lmin ),
                             numpy.log10( self.lmax ),
                             256 )
        
        # convert wavelenght to energy
        El = Ang_to_keV( ll )
        
        # find interval for hard-X normalization
        wE = ( 2. <= El ) & ( El <= 10. )

        # compute emission law normalization
        Lnorm = -1. / trap_int( El[wE], powerlaw_exp_cutoff( El[wE], gamma=gamma, Ecut=Ecut ) )

        # Compute the normalized power-law:
        # - high energy cut-off at 300 keV
        # - low energy cut-off at 50 Angstrom (~0.25 keV)
        ret = numpy.zeros_like( ll )
        wL = ( ll <= 5.e+1 ) 
        ret[wL] = powerlaw_exp_cutoff( El[wL], gamma=gamma, Ecut=Ecut ) * Lnorm

        # store normalized interpolator object for X-ray emission
        self.f_norm_X = lin_interp( numpy.ascontiguousarray( ll ),
                                    numpy.ascontiguousarray( ret ) )

    @abstractmethod
    def set_parameters ( self ) :
        pass
    
    def _set_norm ( self, xx, coeff, fact ) :
        logNorm = poly_N( xx, coeff )
        self.NormStarType = fact * 10.**logNorm * sunL
        return;

    def emission ( self, ll ) :
        return self.NormStarType * self.f_norm_X( ll )

    
class HMXRB ( _XRBase ) :

    coeff = [ 1968.33, -1883.80, 569.44, -62.12, 40.28 ]

    def __init__ ( self, lmin, lmax, **kwargs ) :
        super().__init__( startype = 'hm',
                          gamma    = 2.0,
                          Ecut     = 100.,
                          lmin = lmin, lmax = lmax,
                          **kwargs )
        self.set_parameters()

    def set_parameters ( self, **kwargs ) :

        self.params.update( kwargs )
        super()._set_norm( self.params['Zstar'],
                           type( self ).coeff,
                           self.params['psi'] )
        return;

    
class LMXRB ( _XRBase ) :

    coeff = [ 0.136, 0.425, -0.423, -1.503, 40.276 ]

    def __init__ ( self, lmin, lmax, **kwargs ) :
        super().__init__( startype = 'lm',
                          gamma    = 1.6,
                          Ecut     = 100.,
                          lmin = lmin, lmax = lmax,
                          **kwargs )
        self.set_parameters()

    def set_parameters ( self, **kwargs ) :

        self.params.update( kwargs )
        super()._set_norm( numpy.log10( self.params['age'] * 1.e-9 ),
                           type( self ).coeff,
                           self.params['Mstar'] * 1.e-10 )
        return;

class XRB () :

    def __init__ ( self, lmin, lmax, **kwargs ) :

        self.hm = HMXRB(
            lmin, lmax,
            **{ k : v
                for k, v in kwargs.items()
                if k in xrb_tunables[ 'hm' ] }
        )
        
        self.lm = LMXRB(
            lmin, lmax,
            **{ k : v
                for k, v in kwargs.items()
                if k in xrb_tunables[ 'lm' ] }
        )

    def set_parameters ( self, **kwargs ) :
        self.hm.set_parameters( **{ k : v
                                    for k, v in kwargs.items()
                                    if k in xrb_tunables[ 'hm' ] } )
        self.lm.set_parameters( **{ k : v
                                    for k, v in kwargs.items()
                                    if k in xrb_tunables[ 'lm' ] } )
        return;

    def emission ( self, ll ) :
        return (
            self.hm.emission( ll ) +
            self.lm.emission( ll )
        )
