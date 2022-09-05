""" GalaPy module for combining Simple Stellar Populations(SSP) into a 
    Composite Stellar Population (CSP).



"""

# External imports
import numpy
import os
import warnings

# Internal imports
from .CSP_core import loadSSP, CCSP
from .SYN_core import CSYN
import galapy.internal.globs as GP_GBL
from galapy.internal.data import DataFile

_SSP_LIB = {
    'bc03.basel.chab.extend'  : GP_GBL.bc03_basel_chab_zeros,
    'bc03.basel.chab.refined' : GP_GBL.bc03_basel_chab_zeros_refined,
    'bc03.stelib.chab.extend' : GP_GBL.bc03_stelib_chab_zeros,
    'bc03.stelib.chab.extrap' : GP_GBL.bc03_stelib_chab_extrap,
    'br22.NT'  : GP_GBL.br22_NT,
    'br22.NTL' : GP_GBL.br22_NTL,
    'br22.NT.refined'  : GP_GBL.br22_NT_refined,
    'br22.NTL.refined' : GP_GBL.br22_NTL_refined,
    }

def print_ssp_libs () :
    """ Method for printing on screen the possible choices for the SSP library.
    
    This has the only purpose of listing the possible choices.

    Naming convention: **author.method.imf.filling_schema**

    #. **author**: An achronym for the authors of the SSP library
    #. **method**: if present, shows which method was used to compute the SSP
    #. **imf**: initial mass function
    #. **filling_schema**: all the SSPs' wavelenght domain has been extended to
       span from :math:`1\ \mathring{A}` to :math:`10^{10}\ \mathring{A}`. 
       This code provides the strategy used (:code:`extend` = filled with zeros, 
       :code:`extrap` = extrapolated linearly in the logarithm, 
       :code:`refined` = thinner lambda grid obtained by 
       linearly-interpolating the :code:`extend` equivalent)
    """
    list_libs = ''
    for k in _SSP_LIB.keys() :
        list_libs += f'* {k};\n'
    print( 'Available SSP formats\n'
           '---------------------\n'
           + list_libs )
    return;

class CSP () :
    """ The Composite Stellar Population object handles and manipulates SSP libraries.

    It provides convenient functions for reading the SSPs from the provided formatted 
    binary files. 
    SSPs can be extracted for given age/wavelenght/metallicity.
    Providing an object of type galapy.StarFormationHistory the self.emission function
    computes the luminosity of the resulting Composite Stellar Population.
    
    Keyword Arguments
    -----------------
    ssp_lib : string
      which SSP library to load. The default is :code:`bc03.basel.chab.extend`-
      To see the list of available libraries run 
      :code:`galapy.CompositeStellarPopulation.print_ssp_libs()`
    CCSN : bool
      whether to allow for Core-Collapse Supernova support (default=False).

    See Also
    --------
    :func:`galapy.CompositeStellarPopulation.print_ssp_libs()`
    """

    def __init__ ( self, ssp_lib = 'bc03.basel.chab.extend', CCSN = False ) :
        
        self.ssp_lib = None
        if ssp_lib in set(_SSP_LIB.keys()) :
            self.ssp_lib = ssp_lib
        else :
            raise ValueError( f'SSP library "{ssp_lib}" not available; '
                              'you can see the list of available libraries by running '
                              'galapy.CompositeStellarPopulation.print_ssp_libs()' )
        
        self.l, self.t, self.Z, self.L = loadSSP( DataFile( *_SSP_LIB[self.ssp_lib] ).get_file() )
        self.shape = (self.l.size, self.t.size, self.Z.size)

        self.CCSN = CCSN
        self.core = CCSP( self.l, self.t, self.Z, self.L, self.CCSN )
        self._timetuple = None

        # steal docstrings from C-core:
        self.SSP.__func__.__doc__      = self.core.SSP.__doc__
        self.emission.__func__.__doc__ = self.core.emission.__doc__

    def set_parameters ( self, age, sfh ) :
        """ Function for setting the internal quantities required for combining SSPs into CSPs.

        Parameters
        ----------
        age : float
          The age in years of the CSP.
        sfh : object of type SFH()
          The chosen star formation history. This has to be (or inherit from) 
          the instance of an object of type galapy.StarFormationHistory.SFH().
        """
        from .StarFormationHistory import SFH
        from .SFH_core import CSFH
        
        if isinstance( sfh, CSFH ) :
            self._timetuple = sfh.time_grid( age, self.t, self.Z )
        elif isinstance( sfh, SFH ) :
            self._timetuple = sfh.core.time_grid( age, self.t, self.Z )
        else :
            raise ValueError( 'Argument sfh must be an instance of either '
                              'galapy.StarFormationHistory.SFH or '
                              'galapy.SFH_core.' )
        
        self.core.set_params( *self._timetuple )
        return;
            
    def SSP ( self, il, it, iz ) :
        """ Indexes the SSP library.
        
        Parameters
        ----------
        il : uint or array-like of uints
        it : uint or array-like of uints
        iz : uint or array-like of uints
        
        Returns
        -------
        SSP[il,it,iz] : float or array of floats
          Element (il,it,iz) of the SSP table.
        """
        if il < 0 or it < 0 or iz < 0 :
            raise IndexError( 'Negative index provided' )
        if il >= self.shape[ 0 ] or it >= self.shape[ 1 ] or iz >= self.shape[ 2 ] :
            raise IndexError( 'Luminosity index out of range' )
        return self.core.SSP( il, it, iz )

    def emission ( self, age, sfh, il = None, ftau = None ) :
        """ Computes the CSP emission at given index in the wavelenght-grid.

        It approximates the integral:

        .. math::
           
           L_\\lambda^\\text{CSP}(\\tau') = 
           \\int_0^{\\tau'}\\text{d}\\tau F(\\lambda,\\tau)\\cdot 
           L_\\lambda^\\text{SSP}\\bigl[\\tau, Z_\\ast(\\tau'-\\tau)\\bigr]\\psi(\\tau'-\\tau)


        where :math:`\\tau'` is the age of the CSP, :math:`L_\\lambda^\\text{SSP}[\\tau, Z\\ast]`
        is the luminosity of 
        the Simple Stellar Population at given time :math:`\\tau` and at given stellar metallicity
        :math:`Z_\\ast`, :math:`\\psi(\\tau)` is the Star Formation History and
        :math:`F(\\lambda, \\tau)` is aa attenuating function of choice.
        The latter is passed as the keyword argument :code:`ftau` and should be a matrix
        with dimension :math:`(N_\\lambda, N_\\tau)`, where :math:`N_\\lambda` is the
        size of the wavelenght grid (chosen through the keyword argument :code:`il`) and 
        :math:`N_\\tau` is the size of the time grid of the SSP (:code:`self.shape[1]`).
        :code:`ftau` defaults to :code:`None`, in which case no attenuation is applied.

        Parameters
        ----------
        age : float
          The age in years of the CSP.
        sfh : object of type SFH()
          The chosen star formation history. This has to be (or inherit from) 
          the instance of an object of type galapy.StarFormationHistory.SFH().        
        
        Keyword Arguments
        -----------------
        il : array of int
           array of indexes of the positions in the wavelenght-grid
           for which to compute the emission. The default value is :code:`None`, for which 
           the function will return the emission all over the wavelenght grid.
        Ftau : array
           array containing a function of time to convolve the integrand with 
           (must have same dimension of the `il` array size times the SSP's time-grid size)

        Returns
        -------
        L_CSP : array or scalar float
           the emission of the galaxy's CSP filtered by the function of time :math:`F(\\tau)`
        """

        self.set_parameters( age, sfh )
        if il is None :
            il = numpy.arange( len(self.l), dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True
        if ftau is None :
            ftau = numpy.ones( ( self.t.size * len( il ), ),
                               dtype = numpy.float64 )
        ftau = numpy.ascontiguousarray( ftau.ravel() )

        ret = self.core.emission( il, ftau )
        if scalar_input :
            return ret.item()
        return ret

    def RCCSN ( self, age, sfh ) :
        """
        Parameters
        ----------
        age : float
          The age in years of the CSP.
        sfh : object of type SFH()
          The chosen star formation history. This has to be (or inherit from) 
          the instance of an object of type galapy.StarFormationHistory.SFH().
        
        Returns
        -------
        RCCSN : scalar float
           The Core-Collapse Supernova rate at given age for a given SFH.

        Raises
        ------
        Warning
          if the object has been built without CCSN support.
        """

        if self.CCSN :
            self.set_parameters( age, sfh )
            return self.core.RCCSN()
        else :
            warnings.warn( f"This instance of the {type(self).__name__} "
                           "has been built without CCSN support. "
                           "Build with `CCSN = True` " )
            return 0.
        
    
