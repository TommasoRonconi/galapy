""" The ActiveGalacticNucleus module implements the AGN contribution by loading templates from Fritz et al., 2006
"""

# External imports
import numpy

# Internal imports
from galapy.internal.utils import find_nearest, trap_int
from galapy.internal.interp import lin_interp

_template_tunables = ( 'ct', 'al', 'be', 'ta', 'rm', 'ia' )

_template_pars = {
    'ct' : numpy.array( [ 20, 40, 60 ], dtype = int ),
    'al' : numpy.array( [ 0.0, 2.0, 4.0, 6.0 ], dtype = float ),
    'be' : numpy.array( [ -1.00, -0.75, -0.50, -0.25, 0.00 ], dtype = float ),
    'ta' : numpy.array( [ 0.1, 0.3, 0.6, 1.0, 2.0, 3.0, 6.0, 10.0 ], dtype = float ),
    'rm' : numpy.array( [ 10, 30, 60, 100, 150 ], dtype = int ),
    'ia' : numpy.array( [ 1.e-3, 10.1, 20.1, 30.1, 40.1, 50.1, 60.1, 70.1, 80.1, 89.99 ], dtype = float ),
}

def find_template_par ( key, value ) :
    try :
        return _template_pars[ key ][ find_nearest( _template_pars[ key ], value ) ]
    except KeyError :
        raise AttributeError( f"Parameter '{key}' is not a valid template parameter." )

def agn_build_params ( fAGN, **kwargs ) :
    out = {
        'fAGN' : fAGN,
        'template' : {
            'ct' : 40,
            'al' : 0.,
            'be' : -0.5,
            'ta' : 6.,
            'rm' : 60,
            'ia' : 0.001
        }
    }
    for k in set( out['template'].keys() ).intersection( kwargs.keys() ) :
        out['template'][k] = find_template_par( k, kwargs[k] )
        
    return out

#################################################################################
# AGN class
#################################################################################

class AGN () :
    
    def __init__ ( self, lmin, lmax, pad = 16, fAGN = 1.e-3, **kwargs ) :
        import galapy.internal.globs as GP_GBL
        import os

        self.lmin, self.lmax = lmin, lmax
        self._pad = pad
        
        # common name of all template files
        self._filebase = os.path.join( os.path.dirname( GP_GBL.__file__ ),
                                       GP_GBL.AGN_FILE )
        
        self.params = agn_build_params( fAGN, **kwargs)
        self.load_template()

    def load_template ( self ) :

        # load template from closest file to the parameters chosen
        self.ll, self.tot, self.ther, self.scat, self.disk = \
            numpy.array( [1.e+4,1.,1.,1.,1.] )[:,numpy.newaxis] * \
            numpy.genfromtxt( self._filebase.format( *( self.params[ 'template' ][k]
                                                        for k in _template_tunables  ) ),
                              unpack=True )

        # Extend the wavelenght domain by padding with zeros the emissions
        ltmp = numpy.pad( self.ll, self._pad, constant_values = 0. )
        ltmp[:self._pad+1]  = numpy.logspace( numpy.log10(self.lmin),
                                              numpy.log10(self.ll.min()),
                                              self._pad+1 )
        ltmp[-self._pad-1:] = numpy.logspace( numpy.log10(self.ll.max()),
                                              numpy.log10(self.lmax),
                                              self._pad+1 )
        self.ll  = ltmp
        self.ther = numpy.pad(self.ther, self._pad, constant_values = 0.)
        self.scat = numpy.pad(self.scat, self._pad, constant_values = 0.)
        self.disk = numpy.pad(self.disk, self._pad, constant_values = 0.)
        self.tot  = numpy.pad(self.tot,  self._pad, constant_values = 0.)

        # Compute constant normalization factor for the current template 
        self.Lnorm = 1. / trap_int( self.ll, self.tot )

        # Build a linear normalized interpolator
        self.f_norm_tot = lin_interp( numpy.ascontiguousarray( self.ll ),
                                      numpy.ascontiguousarray( self.tot * self.Lnorm ) )
        return;
        
    def set_parameters ( self, fAGN = None, **kwargs ) :

        if fAGN is not None :
            self.params.update( fAGN = fAGN )
        if len(kwargs) > 0 :
            for k in set( self.params['template'].keys() ).intersection( kwargs.keys() ) :
                self.params['template'][k] = find_template_par( k, kwargs[k] )
            self.load_template()
        return;

    def emission ( self, ll, Lscale ) :
        ll = numpy.ascontiguousarray( ll, dtype = numpy.float64 )
        fact = self.params['fAGN']/(1-self.params['fAGN'])
        return fact * Lscale * self.f_norm_tot( ll )
    
        

