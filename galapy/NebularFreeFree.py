"""Implements the Nebular Continuum Emission (Bremsstrahlung) component.
"""

# External imports
import numpy

# Internal imports
from .NFF_core import CNFF

nff_tunables = ( 'Zgas', 'Zi' )

def nff_build_params ( **kwargs ) :
    
    out = {
        'Zgas' : 0.01,
        'Zi'   : 1.,
    }
    for k in set( out.keys() ).intersection( kwargs.keys() ) :
        out[ k ] = kwargs[ k ]

    return out

class NFF () :
    """ Class wrapping the C-core implementation of the Nebular Free-Free emission type.   
    
    Parameters
    ----------
    ll : float array
      the wavelength grid where the free-free emission is computed
    **kwargs : dictionary, optional
      All possible free-parameters of the model, refer to the
      Keyword Arguments section.
    
    Keyword Arguments
    -----------------
    Zgas : float
      (Optional, default = 0.01) Absolute metallicity of gas in the hosting galaxy.
    Zi : float
      (Optional, default = 1) Average atomic number of gas in the hosting galaxy.
    """
    
    def __init__ ( self, ll, **kwargs ) :

        self.l = numpy.ascontiguousarray( ll )
        self.core = CNFF( self.l )
        self.params = nff_build_params( **kwargs )
        self.set_parameters()

        
    def set_parameters ( self, **kwargs ) :
        r""" Function for setting the parameters of the model.

        Parameters
        ----------
        **kwargs : dictionary, optional
          All possible free-parameters of the model, refer to the
          Keyword Arguments section.
    
        Keyword Arguments
        -----------------
        Zgas : float
          Absolute metallicity of gas in the hosting galaxy.
        Zi : float
          Average atomic number of gas in the hosting galaxy.
        
        Returns
        -------
        : None
        """

        self.params.update( kwargs )
        self.core.set_params( numpy.asarray( [
            self.params[k]
            for k in nff_tunables
        ], dtype=float) )
        return;

    def gff ( self, il = None, Te = None ) :
        r"""Computes the Gaunt Factor as given in Draine, 2011 (Eq.10.9):

        .. math ::
            
            g_\text{NFF}=\ln\left\{\exp\left[5.96-\frac{\sqrt{3}}{\pi}\, 
            \ln\left( Z_i\, \frac{\nu_\lambda}{G\text{Hz}}\, \left(\frac{T_e}{10^4\, K}\right)^{-1.5}\right)\right]+
            \exp(1)\right\}
        
        where :math:`T_e` is the electron temperature (function argument), 
        :math:`Z_i` the average atomic number (class parameter) and 
        :math:`\nu_\lambda` is the frequency at given wavelength 
        (:math:`\lambda` as defined in the grid used to instantiate the class)
        
        Parameters
        ----------
        il : 1d array or scalar integer
          indices in the instance wavelength grid. If None, the gaunt factor
          is computed over all the wavelength grid.
        Te : float
          electron temperature of the medium, if None, the internal function :code:`NFF.Te()` is used.
        
        Returns
        -------
        : 1d array or scalar float
          Gaunt Factor at given wavelength
        """
        
        if il is None :
            il = numpy.arange( self.l.size, dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True
        if Te is None :
            Te = self.Te()

        ret = numpy.array( [ self.core.gff( _i, Te ) for _i in il ] )
        if scalar_input :
            return ret.item()
        return ret

    def Te ( self, Zgas = None ) :
        r"""Electron temperature as empirically derived in Vega et al., 2008:
        
        .. math ::
            
            T_e = \exp\left\{3.89-0.4802\, \log(Z_\text{gas}/0.02)-0.0205\,[\log(Z_\text{gas}/0.02)]^2\right\}
        
        where :math:`Z_\text{gas}` is the absolute gas metallicity.
        
        Parameters
        ----------
        Zgas : float
          (optional, default = same as parameters set) the average absolute metallicity of the gas.
          If None, it will use the internally stored parameter (can be set with :code:`NFF.set_parameters( Zgas = ... )`)
        
        Returns
        -------
        : float
          Electron temperature at given metallicity
        """
        
        if Zgas is None :
            Zgas = self.params[ 'Zgas' ]
        Zgas = numpy.asarray( Zgas )
        scalar_input = False
        if Zgas.ndim == 0 :
            Zgas = Zgas[None] # makes il 1D
            scalar_input = True
            
        ret = numpy.array( [ self.core.Te( _Z ) for _Z in Zgas ] )
        if scalar_input :
            return ret.item()
        return ret

    def emission ( self, Q_H, il = None ) :
        r"""Intrinsic Free-free emission as given in Draine, 2011 (Chapter 10):
        
        .. math::
            
            L_\text{NFF}(\lambda,\tau) \approx 1.8\times 10^{-27}\, \text{erg}\,s^{-1} \text{Hz}^{-1}\cdot
            \dfrac{\mathcal{Q}_\text{H}(\tau)}{s^{-1}}\, \left(\dfrac{T_e}{10^4\, K}\right)^{0.3}\, g_\text{NFF}\, 
            \exp\left(\dfrac{-h_P\, \nu_\lambda}{k_B\, T_e}\right)
        
        Parameters
        ----------
        Q_H : float
          intrinsic photo-ionisation rate, can be computed from the
          intrinsic stellar luminosity:
          
          .. math ::
              
              \mathcal{Q}_\text{H}(\tau) = \int_0^{912\mathring{A}}\text{d}\lambda\, 
              \frac{L_\text{CSP}^\text{i}(\lambda,\tau)}{h_P\, \nu_\lambda}
          
          where :math:`h_P` is the Planck constant and :math:`\nu_\lambda` is the frequency
          at given wavelength :math:`\lambda`.
          The intrinsic stellar luminosity, :math:`L_\text{CSP}^\text{i}(\lambda,\tau)`,
          can be computed using function :code:`galapy.CompositeStellarPopulation.CSP.emission()`.
        il : 1d array or scalar integer
          indices in the instance wavelength grid. If None, the emission
          is computed over all the wavelength grid.
        
        Returns
        -------
        : 1d array or scalar float
          intrinsic emission due to Free-Free transitions.
        """
        
        if il is None :
            il = numpy.arange( len(self.l), dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True
        
        ret = self.core.emission( il, Q_H )
        if scalar_input :
            return ret.item()
        return ret
        
    
