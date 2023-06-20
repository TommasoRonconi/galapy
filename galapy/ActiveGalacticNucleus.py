r""" The ActiveGalacticNucleus module implements the AGN contribution by loading templates from Fritz et al., 2006.

Templated emission is divided into 3 components:

- Accretion disk around the central SMBH (accessible through the instance variable :code:`AGN().disk`)
- Scattered emission by the surrounding dusty torus (accessible through the instance variable :code:`AGN().scat`)
- Thermal dust emission associated to the dusty torus (accessible through the instance variable :code:`AGN().ther`)

The templates are computed accounting for the variation of 6 structural parameters:

- the covering angle of the torus :math:`\Gamma`, expressed in terms of half the aperture-angle with respect to the equatorial plane:

  .. math ::

    \Theta = 90^\circ - \dfrac{\Gamma}{2}

- the dust density distribution in spherical coordinates:

  .. math ::

    \rho \propto r^\beta e^{- \alpha|\cos\theta|}

  parameterized in terms of :math:`\alpha` and :math:`\beta`;
- the ratio :math:`R^\text{AGN}_\text{torus}` between the maximum to minimum radii of the dusty torus;
- the optical depth :math:`\tau^\text{AGN}_{9.7}` at :math:`9.7\ \mu m`;
- the viewing angle :math:`\Psi^\text{AGN}_\text{los}`, i.e. the angle between the rotation axis and the line-of-sight.
  (note that by assuming the unifed AGN model, with :math:`\Psi^\text{AGN}_\text{los} = 90^\circ` the geometry of 
  the object corresponds to a type 1 AGN, while :math:`\Psi^\text{AGN}_\text{los} = 0^\circ` to a type 2 AGN)
  
The above parameters in our implementation are controlled through the above parameters following the original naming schema used by the authors and summarized in the following table:

+------------------------------------+------------+-------------------------+----------------------------------------------------------------+
| Parameter                          | Key        | :math:`N_\text{values}` | Values                                                         |
+====================================+============+=========================+================================================================+
| :math:`\Theta`                     | :code:`ct` | 3                       | (:math:`20^\circ,\ 40^\circ,\ 60^\circ`)                       |
+------------------------------------+------------+-------------------------+----------------------------------------------------------------+ 
| :math:`\alpha`                     | :code:`al` | 4                       | (:math:`0,\ 2,\ 4,\ 6`)                                        |
+------------------------------------+------------+-------------------------+----------------------------------------------------------------+
| :math:`\beta`                      | :code:`be` | 5                       | (:math:`-1,\ -0.75,\ -0.5,\ -0.25,\ 0`)                        |
+------------------------------------+------------+-------------------------+----------------------------------------------------------------+
| :math:`R^\text{AGN}_\text{torus}`  | :code:`rm` | 5                       | (:math:`10,\ 30,\ 60,\ 100,\ 150`)                             |
+------------------------------------+------------+-------------------------+----------------------------------------------------------------+
| :math:`\tau^\text{AGN}_{9.7}`      | :code:`ta` | 8                       | (:math:`0.1,\ 0.3,\ 0.6,\ 1,\ 2,\ 3,\ 6,\ 10`)                 |
+------------------------------------+------------+-------------------------+----------------------------------------------------------------+
| :math:`\Psi^\text{AGN}_\text{los}` | :code:`ia` | 10                      | :math:`0^\circ` to :math:`90^\circ` with step :math:`10^\circ` |
+------------------------------------+------------+-------------------------+----------------------------------------------------------------+

Which results in 24000 available templates. The user can nevertheless set the parameters to whatever value and the template whose parameters are closer to the selection is chosen.
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
    """ Given a key-value pair returns the parameter corresponding to `key`
    included in the template-library and closest to `value`
    """
    try :
        return _template_pars[ key ][ find_nearest( _template_pars[ key ], value ) ]
    except KeyError :
        raise AttributeError( f"Parameter '{key}' is not a valid template parameter." )

def agn_build_params ( fAGN, **kwargs ) :
    """ Standard function for building the parameters dictionary of class AGN()
    """
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
    r""" AGN component class
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
      minimum and maximum values of the wavelength-domain. The original templates
      are computed within the 10-10^7 Angstrom interval. At build time this domain
      is extended with `pad` padding values before and after the limits of the 
      template's wavelength domain to match the requested limits.
    pad : integer
      number of padding values to match the requested wavelength domain
    fAGN : float
      fraction of the total reference emission
    do_Xray : bool
      If :code:`True` compute the X-ray template along with the closest template for
      the lower energy part of the spectrum (default is :code:`True`)

    Keyword Arguments
    -----------------
    ct : scalar
      the covering angle of the torus :math:`\Gamma`, expressed in terms of half the aperture-angle with respect to the equatorial plane
    al : scalar
      parameter regulating the dust density distribution in spherical coordinates, corresponds to :math:`\alpha`
    be : scalar
      parameter regulating the dust density distribution in spherical coordinates, corresponds to :math:`\beta`
    ta : scalar
      the optical depth at :math:`9.7\ \mu m`
    rm : scalar
      the ratio between the maximum to minimum radii of the dusty torus      
    ia : scalar
      the inclination angle between the rotation axis of the AGN torus and the line of sight.
      Note that by assuming the unifed AGN model, with :code:`ia = 90` the geometry of the object corresponds to a type 1 AGN, 
      while :code:`ia = 0` to a type 2 AGN.
    """
    
    def __init__ ( self, lmin, lmax, pad = 16, fAGN = 1.e-3, do_Xray = True, **kwargs ) :
        import galapy.internal.globs as GP_GBL
        import os

        # check inputs
        if not fAGN < 1.0 :
            raise RuntimeError( 'trying to set fAGN>=1.0 not allowed.' )
        
        # store the argument variables
        self.lmin, self.lmax = lmin, lmax
        self._pad = pad
        self.do_Xray = do_Xray
        
        # common name of all template files
        self._filebase = GP_GBL.AGN_FILE

        # build the parameter dictionary
        self.params = agn_build_params( fAGN, **kwargs)

        # load the template with physics nearest to the parameters value 
        self.load_template()

        # also compute the X-ray template if requested (default=True)
        if self.do_Xray :
            self.compute_X_template()

    def load_template ( self ) :
        """ Loads the template corresponding to the current parameters (mostly intended for internal usage).
        """

        # load template from closest file to the parameters chosen
        filename = self._filebase.format( *( self.params[ 'template' ][k]
                                             for k in _template_tunables  ) )
        self.ll, self.tot, self.ther, self.scat, self.disk = \
            numpy.array( [1.e+4,1.e-4,1.e-4,1.e-4,1.e-4] )[:,numpy.newaxis] * \
            numpy.genfromtxt( DataFile( filename, GP_GBL.AGN_DIR ).get_file(),
                              unpack=True )

        # Extend the wavelength domain by padding with zeros the emissions
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
        """ Pre-computes the not-normalized and not-bolometric-corrected X-ray spectrum. 
        """

        # generate wavelength grid
        ll = numpy.logspace( numpy.log10( self.lmin ),
                             numpy.log10( self.lmax ),
                             256 )
        # convert wavelength to energy
        El = Ang_to_keV( ll )

        if El.max() < 2. :
            raise RuntimeError( "Cannot build the X-ray spectrum for "
                                "a wavelength grid starting at lambda > "
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
        r""" Computes the bolometric correction from Duras et al., 2020 (Eq. 2):
        
        .. math ::
        
          \dfrac{L_\text{bol}^\text{AGN}}{L_\text{X}^\text{AGN}} \approx 10.96\;\biggl[1 + \biggl(\dfrac{\log L_\text{bol}^\text{AGN}/L_\odot}{11.93}\biggr)^17.79\biggr]

        Parameters
        ----------
        Lref : scalar
          Bolometric luminosity of the AGN in units of solar luminosities :math:`L_\odot`
        
        Returns
        -------
        : scalar
          Value of the bolometric correction.
        """
        return 10.96 * ( 1. + ( numpy.log10( Lref ) / 11.93 )**17.79 )
        
    def set_parameters ( self, fAGN = None, **kwargs ) :
        r""" Sets the internal parameters defining the emission coming from the AGN.
        All the parameters that are not passed to this function will be left to the
        value they already have.
        
        Parameters
        ----------
        fAGN : scalar
          sets the normalization of the templated emission. The parameters fAGN 
          provides the amount of energy radiated by the AGN in units of some other
          emission contributor, typically the energy radiated by dust in the same band.
        **kwargs 
          These are the parameters passed to function :code:`find_template_par`. Used
          to find the nearest AGN emission template in the Fritz et al. (2006) library.
        """

        if fAGN is not None :
            if fAGN < 1.0 :
                self.params.update( fAGN = fAGN )
            else :
                raise RuntimeError( 'trying to set fAGN>=1.0 not allowed.' )
        if len(kwargs) > 0 :
            for k in set( self.params['template'].keys() ).intersection( kwargs.keys() ) :
                self.params['template'][k] = find_template_par( k, kwargs[k] )
            self.load_template()
        return;

    def emission ( self, ll, Lref ) :
        """ Computes the emission coming from the AGN adding the X-ray part 
        if the internal variable do_Xray has been set to :code:`True` at build time
        
        Parameters
        ----------
        ll : sequence
          wavelength list or array or iterable
        Lref : float
          bolometric reference luminosity used to compute the bolometric correction
          of the high-energy part of the spectrum
        
        Returns
        -------
        : numpy-array
          Luminosity per unit wavelength in units of solar luminosities emitted by
          the AGN.
        """
        ll = numpy.ascontiguousarray( ll, dtype = numpy.float64 )
        fact = self.params['fAGN']/(1-self.params['fAGN'])
        if self.do_Xray :
            return fact * Lref * ( self.f_norm_tot( ll ) +
                                   self.f_norm_X( ll ) *
                                   self.X_bolometric_correction( Lref ) )
        return  fact * Lref * self.f_norm_tot( ll )
    
        

