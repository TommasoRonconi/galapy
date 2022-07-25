""" The ActiveGalacticNucleus module implements the AGN contribution by loading templates from Fritz et al., 2006
"""

# External imports
import numpy

# Internal imports
import galapy.internal.globs as GP_GBL
from galapy.internal.utils import find_nearest, trap_int, powerlaw_exp_cutoff
from galapy.internal.constants import Ang_to_keV
from galapy.internal.interp import lin_interp
from galapy.internal.data import DataFile

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
    """ AGN component class
    This class implements the templated emission from an Active Galactic Nucleus
    within the galaxy.
    We use templates from Fritz et al., 2006 to model the AGN emission normalized
    at build time. The total emission is then rescaled to a user-provided total emission.
    
    The final emission can be expressed as
    
    .. math::
    
       L_\text{AGN}(\lambda) = \dfrac{f_\text{AGN}}{1-f_\text{AGN}}L_\text{AGN}^\text{norm}(\lambda)\cdot L_\text{ref}
    
    where 
    :math:`f_\text{AGN}` is the fraction of the reference emission, :math:`L_\text{ref}`, 
    radiated by the AGN and
    :math:`L_\text{AGN}^\text{norm}` is the normalized total templated emission.
    
    Parameters
    ----------
    lmin, lmax : float
      minimum and maximum values of the wavelenght-domain. The original templates
      are computed within the 10-10^7 Angstrom interval. At build time this domain
      is extended with `pad` padding values before and after the limits of the 
      template's wavelenght domain to match the requested limits.
    pad : integer
      number of padding values to match the requested wavelenght domain
    fAGN : float
      fraction of the total reference emission

    Keyword Arguments
    -----------------
    ct : scalar
    al : scalar
    be : scalar
    ta : scalar
    rm : scalar
    ia : scalar   
    """
    
    def __init__ ( self, lmin, lmax, pad = 16, fAGN = 1.e-3, Xray = True, **kwargs ) :
        import galapy.internal.globs as GP_GBL
        import os

        # store the argument variables
        self.lmin, self.lmax = lmin, lmax
        self._pad = pad
        self._Xray = Xray
        
        # common name of all template files
        self._filebase = GP_GBL.AGN_FILE

        # build the parameter dictionary
        self.params = agn_build_params( fAGN, **kwargs)

        # load the template with physics nearest to the parameters value 
        self.load_template()

        # also compute the X-ray template if requested (default=True)
        if self._Xray :
            self.compute_X_template()

    def load_template ( self ) :

        # load template from closest file to the parameters chosen
        filename = self._filebase.format( *( self.params[ 'template' ][k]
                                             for k in _template_tunables  ) )
        self.ll, self.tot, self.ther, self.scat, self.disk = \
            numpy.array( [1.e+4,1.e-4,1.e-4,1.e-4,1.e-4] )[:,numpy.newaxis] * \
            numpy.genfromtxt( DataFile( filename, GP_GBL.AGN_DIR ).get_file(),
                              unpack=True )

        # Extend the wavelenght domain by padding with zeros the emissions
        ltmp = numpy.pad( self.ll, self._pad, constant_values = 0. )
        ltmp[:self._pad+1]  = numpy.logspace( numpy.log10(self.lmin),
                                              numpy.log10(self.ll.min()),
                                              self._pad+1 )
        ltmp[-self._pad-1:] = numpy.logspace( numpy.log10(self.ll.max()),
                                              numpy.log10(self.lmax),
                                              self._pad+1 )
        self.ll   = ltmp
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

    def compute_X_template ( self ) :

        # generate wavelenght grid
        ll = numpy.logspace( numpy.log10( self.lmin ),
                             numpy.log10( self.lmax ),
                             256 )
        # convert wavelenght to energy
        El = Ang_to_keV( ll )

        if El.min() > 2. :
            raise RuntimeError( "Cannot build the X-ray spectrum for "
                                "a wavelenght grid starting at lambda > "
                                "6 Angstrom ~ 2 keV! "
                                "Set a smaller `lmin` value." )
        
        # find interval for hard-X normalization
        wE = ( 2. <= El ) & ( El <= 10. )

        # compute emission law normalization
        Lnorm = -1. / trap_int( El[wE], powerlaw_exp_cutoff( El[wE], gamma=1.8, Ecut=3.e+2 ) )

        # Compute the normalized power-law:
        # - high energy cut-off at 300 keV
        # - spectral index fixed to 1.8
        # - low energy cut-off at 50 Angstrom (~0.25 keV)
        ret = numpy.zeros_like( ll )
        wL = ( ll <= 5.e+1 ) 
        ret[wL] = powerlaw_exp_cutoff( El[wL], gamma=1.8, Ecut=3.e+2 ) * Lnorm

        # store normalized interpolator object for X-ray emission
        self.f_norm_X = lin_interp( numpy.ascontiguousarray( ll ),
                                    numpy.ascontiguousarray( ret ) )
        return;
        

    def X_bolometric_correction ( self, Lref ) :
        """ Duras et al., 2020
        """
        return 10.96 * ( 1. + ( numpy.log10( Lref ) / 11.93 )**17.79 )
        
    def set_parameters ( self, fAGN = None, **kwargs ) :

        if fAGN is not None :
            self.params.update( fAGN = fAGN )
        if len(kwargs) > 0 :
            for k in set( self.params['template'].keys() ).intersection( kwargs.keys() ) :
                self.params['template'][k] = find_template_par( k, kwargs[k] )
            self.load_template()
        return;

    def emission ( self, ll, Lref ) :
        ll = numpy.ascontiguousarray( ll, dtype = numpy.float64 )
        fact = self.params['fAGN']/(1-self.params['fAGN'])
        if self._Xray :
            return fact * Lref * ( self.f_norm_tot( ll ) +
                                   self.f_norm_X( ll ) *
                                   self.X_bolometric_correction( Lref ) )
        return  fact * Lref * self.f_norm_tot( ll )
    
        

