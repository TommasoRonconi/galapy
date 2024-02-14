""" The Synchrotron module implements a generic parameterized synchrotron emission.
"""

########################################################################################

# External imports
import numpy

# Internal imports
from galapy.internal.constants import clight
from .SYN_core import CSYN

########################################################################################

__all__ = [ 'SYN', 'SNSYN' ]

syn_tunables = ( 'alpha_syn', 'nu_self_syn' )
"""Tunable parameters for synchrotron emission."""

def syn_build_params ( **kwargs ) :
    """ Builds the parameters dictionary for the generic synchrotron.
    """
    
    out = {
        'alpha_syn'   : 0.75,
        'nu_self_syn' : 0.2,  # [GHz]
    }
    for k in set( out.keys() ).intersection(kwargs.keys()) :
        out[k] = kwargs[k]
        
    return out

class SYN () :
    r"""Class wrapping the C-core implementation of a generic Synchrotron emission type.

    The shape of the emission is given by
    
    .. math::
       
       L_\text{syn}(\lambda,\tau) = L_\text{norm} \cdot \left(\dfrac{\nu_\lambda}{G\text{Hz}}\right)^{-\alpha_\text{syn}}\cdot
       \left[1+\left(\dfrac{\nu}{20\, G\text{\text{Hz}}}\right)^{0.5}\,\right]^{-1}\, F[\tau_\text{syn}(\nu_\lambda)]

    where :math:`\alpha_\text{syn}` is the spectral index, the term in square brackets takes into account spectral-ageing effects, 
    and the function :math:`F(x)=(1-e^{-x})/x` incorporates synchrotron self-absorption in terms of the optical depth 
    :math:`\tau_{\rm sync}\approx (\nu_\lambda/\nu_\text{self})^{-\alpha_\text{syn}-5/2}` that is thought to become relevant 
    at frequencies :math:`\nu\lesssim \nu_\text{self}\approx 0.2\ [GHz]`.
    
    Parameters
    ----------
    ll : float array
      the wavelength grid where the synchrotron emission is computed
    **kwargs : keyword arguments
      Parameters for the synchrotron model.
      The available free-parameters are listed in the galapy.Synchrotron.syn_tunables tuple

    Keyword Arguments
    -----------------
    alpha_syn : float
       synchrotron spectral index (i.e. :math:`\alpha_\text{syn}`)
    nu_self_syn : float
       synchrotron self-absorption threshold (i.e. :math:`\nu_\text{self}`)
    """

    def __init__ ( self, ll, **kwargs ) :
        """Initialize SYN object."""

        self.l = numpy.ascontiguousarray( ll )
        self.AngUnit = clight['A/s'] / self.l**2
        self.core = CSYN( self.l )
        self.params = syn_build_params( **kwargs )
        self.set_parameters()
        
    def set_parameters ( self, **kwargs ) :
        r""" Function for setting the parameters of the model.
        
        Parameters
        ----------
        **kwargs : keyword arguments
            Parameters to set for the synchrotron model.
        
        Returns
        -------
        : None
        """

        self.params.update( kwargs )
        self.core.set_params( numpy.asarray( [
            self.params[k]
            for k in syn_tunables
        ], dtype=float) )
        return;

    def opt_depth ( self, il = None ) :
        """Compute optical depth on the default wavelength grid.
        
        Parameters
        ----------
        il : array-like or int or None, optional
            Indices of wavelengths in the default grid to compute optical depth for.
        
        Returns
        -------
        numpy.ndarray or float
            Array of optical depths or a single optical depth if il is a singl integer.
            If il is None return optical depth on the whole wavelength grid.
        """

        if il is None :
            il = numpy.arange( self.l.size, dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True

        ret = numpy.array( [ self.core.opt_depth( _i ) for _i in il ] )
        if scalar_input :
            return ret.item()
        return ret

    def energy ( self, SynNorm, il = None ) :
        """Returns the normalized total energy radiated at given wavelength.
        
        Parameters
        ----------
        SynNorm : float
            Normalization factor for synchrotron energy.
        il : array-like or int or None, optional
            Indices of wavelengths in the default grid  to compute energy for.
        
        Returns
        -------
        numpy.ndarray or float
            Array of energy values or a single energy value if il is a single integer.
            If il is None return energy on the whole wavelength grid.
        """
        
        if il is None :
            il = numpy.arange( len(self.l), dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True
        
        ret = self.core.energy( il, SynNorm )
        if scalar_input :
            return ret.item()
        return ret

    def emission ( self, SynNorm, il = None ) :
        """Compute emission.
        
        Parameters
        ----------
        SynNorm : float
            Normalization factor for synchrotron energy.
        il : array-like or int or None, optional
            Indices of wavelengths to compute emission for.
        
        Returns
        -------
        numpy.ndarray or float
            Array of emission values or a single emission value if il is a single integer
            If il is None return the emission on the whole wavelength grid.
        """
        
        return self.energy( SynNorm, il ) * self.AngUnit

########################################################################################
    
class SNSYN ( SYN ) :
    r"""Specialised class for handling Synchrotron emission due to Core-Collapse Supernovae.

    The implementation inherits from the generic galapy.Synchrotron.SYN class and computes the
    emission as

    .. math::
       
       L_\text{syn}^\text{corr}(\lambda,\tau)=\dfrac{L_\text{syn}(\lambda,\tau)}{1+[L_\text{syn}^0/L_\text{syn}(\lambda,\tau)]^\zeta}
    
    with :math:`L_\text{syn}^0 = 3\cdot10^{28}\ [erg]` (fixed by the class-attribute SNSYN.Lsyn0), 
    :math:`\zeta=2` (fixed by the class-attribute SNSYN.zeta) and with emission shape normalised
    by :math:`L_\text{norm} = \text{const} \cdot \dfrac{\mathcal{R}_\text{CCSN}(\tau)}{\text{yr}^{-1}}`.
    The :math:`\text{const} = 10^{30}\ [erg]` is fixed (by the class-attribute SNSYN.NormFact) 
    while the value :math:`\mathcal{R}_\text{CCSN}` is the production rate of CC Supernovae and is 
    a free parameter of the model and can be computed
    with :py:func:`galapy.CompositeStellarPopulation.CSP.RCCSN`
    
    Parameters
    ----------
    ll : float array
      the wavelength grid where the synchrotron emission is computed
    **kwargs : keyword arguments
      Parameters for the synchrotron model.

    Keyword Arguments
    -----------------
    alpha_syn : float
       synchrotron spectral index (i.e. :math:`\alpha_\text{syn}`)
    nu_self_syn : float
       synchrotron self-absorption threshold (i.e. :math:`\nu_\text{self}`)
    RCCSN : float 
       Core-collapse Super-Nova production rate (unit :math:`\text{yr}^{-1}`
    """

    NormFact = 1.e+30 # [ erg s^-1 Hz^-1 ]
    Lsyn0    = 3.e+28 # [ erg s^-1 Hz^-1 ]
    zeta     = 2.     # [ adimensional ] 

    def __init__ ( self, *args, **kwargs ) :
        """Initialize SNSYN object."""

        super().__init__( *args, **kwargs )
        self.params['RCCSN'] = 0.
        self.params.update(kwargs)

    def emission ( self, il = None ) :
        """Compute emission.
        
        Parameters
        ----------
        il : array-like or int or None, optional
            Indices of wavelengths to compute emission for.
        
        Returns
        -------
        numpy.ndarray or float
            Array of emission values or a single emission value if il is a single integer
            If il is None return the emission on the whole wavelength grid.
        """
        
        Lsyn = super().energy( self.params['RCCSN'] * SNSYN.NormFact, il )
        wn0 = Lsyn > 0.
        Lsyn[wn0] /= 1 + ( SNSYN.Lsyn0 / Lsyn[wn0] )**SNSYN.zeta

        if il is None :
            return Lsyn * self.AngUnit
        return Lsyn * self.AngUnit[ il ]

########################################################################################
