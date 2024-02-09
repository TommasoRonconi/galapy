""" GalaPy module for combining Simple Stellar Populations (SSP) into a 
Composite Stellar Population (CSP).
"""

# External imports
import numpy
import os
import warnings

# Internal imports
from .CSP_core import loadSSP as _loadSSP, CCSP
from .SYN_core import CSYN
import galapy.internal.globs as GP_GBL
from galapy.internal.data import DataFile

def _recursive_list_ssp_libs ( root, pathline = [],
                               outlist = [], outpath = None ) :
    if os.path.isdir(root) :
        rootpath, subd, files = next(os.walk(root))
        for f in files :
            outlist += ['.'.join(pathline+f.split('.')[:-1])]
            if outpath is not None :
                outpath += [ os.path.join(rootpath, f) ]
        for s in subd : _recursive_list_ssp_libs(os.path.join(root, s),
                                                 pathline+[s], outlist,
                                                 outpath)
    return;

def print_ssp_libs () :
    """Print on screen a list of the available SSP libraries.
    
    This has the only purpose of listing the possible choices.

    Naming convention: **author.method[.imf.filling_schema]**

    #. **author**: An achronym for the authors of the SSP library
    #. **method**: if present, shows which method was used to compute the SSP
    #. **imf**: if present shows the initial mass function used to compute the SSP
    #. **filling_schema**: all the SSPs' wavelength domain has been extended to
       span from :math:`1\ \mathring{A}` to :math:`10^{10}\ \mathring{A}`. 
       This code provides the strategy used (plain = filled with zeros, 
       :code:`extrap` = extrapolated linearly in the logarithm, 
       :code:`refined` = thinner lambda grid obtained by 
       linearly-interpolating the plain equivalent)
    """
    from galapy.configuration import rcParams
    outlist = []
    for path in rcParams['datapath'] :
        _recursive_list_ssp_libs( root = os.path.join( path, *GP_GBL.SSP_DIR ),
                                  pathline = [], outlist = outlist )
    list_libs = ''
    for k in outlist :
        list_libs += f'* {k};\n'
    print( 'Available SSP formats\n'
           '---------------------\n'
           + list_libs )

def list_ssp_libs ( return_path = False ) :
    """Return a list of the available SSP libraries.
    
    Parameters
    ----------
    return_path : bool
        (default = False) if True, also return a list
        with absolute paths to each of the available SSP libs.
    
    Returns
    -------
    : list
        available SSP libraries
    : list
        available SPP libraries as a list of absolute paths.
    """
    
    from galapy.configuration import rcParams
    outlist = []
    outpath = None
    if return_path :
        outpath = []
    for path in rcParams['datapath'] :
        _recursive_list_ssp_libs( root = os.path.join( path, *GP_GBL.SSP_DIR ),
                                  pathline = [], outlist = outlist,
                                  outpath = outpath )
    if outpath is not None :
        return outlist, outpath
    return outlist

_SSP_LIB = dict( zip( *list_ssp_libs(True) ) )

def load_SSP_table ( which ) :
    """Load a SSP library corresponding to a given name.
    """

    filepath, filename = os.path.split( _SSP_LIB[ which ] )
    path_line = []
    level = None
    while level != GP_GBL.DATA_DIR :
        filepath, level = os.path.split( filepath )
        path_line = [ level ] + path_line
    
    return _loadSSP( DataFile( filename, tuple(path_line) ).get_file() )

def format_SSP_table ( L, lidx, tidx, Zidx, lsize, tsize, Zsize, flat = True ) :
    """Given an input SSP table returns a version with the formatting internally
    required by the library.
        
    Parameters
    ----------
    L : array-like
        input SSP table. It can be either a continuous array indexed on 3 coordinates
        (i.e. wavelength, age and metallicity) in whatever ordering, or a 3-dimensional
        matrix for which the 3 dimensions run over the wavelength, age and metallicity 
        (also in this case the order is not important).
    lidx : integer
        axis number in the array `L` corresponding to the wavelength dimension 
    tidx : integer
        axis number in the array `L` corresponding to the age dimension 
    Zidx : integer
        axis number in the array `L` corresponding to the metallicity dimension 
    lsize : integer
        dimension of the wavelength grid 
    tsize : integer
        dimension of the age grid 
    Zsize : integer
        dimension of the metallicity grid
    flat : bool
        if True, the input array is assumed to be flattened, if False it should instead
        be a 3D matrix with axes corresponding to the 3 grids in wavelength, age and metallicity
        
    Returns
    -------
    newL : 1d contiguous array
        output SSP table in the correct format
    """
    input_shape = numpy.empty(3, dtype=int)
    input_shape[lidx] = lsize
    input_shape[tidx] = tsize
    input_shape[Zidx] = Zsize
    if flat :
        L = numpy.array(L).reshape(input_shape)
    if tuple(input_shape) != L.shape :
        raise RuntimeError( f'Wrong input shape, expected {tuple(input_shape)}, got {L.shape}' )
        
    return numpy.ascontiguousarray( 
        numpy.transpose(L, axes=(Zidx, lidx, tidx)).ravel()
    )

def store_SSP_table ( outfile, l, t, Z, L, endianism = 'little', force = False ) :
    """Save a new SSP table
    
    Parameters
    ----------
    outfile : string
        path to output file name.
    l : array-like
        wavelength grid
    t : array-like
        age grid
    Z : array-like
        metallicity grid
    L : array-like
        flattened SSP table in the correct format accepted by the library
        (check function ``format_SSP_table()``)
    endianism : string
        endianism of your machine, you can check this by calling ``import sys; print(sys.byteorder)``
    force : bool
        whether to overwrite an already existing ``outfile``
    """
    from sys import byteorder
    if endianism != byteorder :
        raise RuntimeError('wrong endianism selected')
    if os.path.isfile(outfile) and not force :
        warnings.warn( 
            f'Selected filename {outfile} is an existing file, to overwrite it set parameter `force=True`')
        return;
    if not os.path.isdir(os.path.dirname(outfile)) :
        raise RuntimeError(
            'The selected output file name is pointing to path that does not exist in the filesystem'
        )
    l = numpy.ascontiguousarray(l)
    t = numpy.ascontiguousarray(t)
    Z = numpy.ascontiguousarray(Z)
    L = numpy.ascontiguousarray(L)
    with open( outfile, 'wb' ) as f :
        # write wavelength block
        f.write( l.size.to_bytes(8, endianism) )
        f.write( l.tobytes() )
        # write age block
        f.write( t.size.to_bytes(8, endianism) )
        f.write( t.tobytes() )
        # write metallicity block
        f.write( Z.size.to_bytes(8, endianism) )
        f.write( Z.tobytes() )
        # write luminosity block
        f.write( L.size.to_bytes(8, endianism) )
        f.write( L.tobytes() )
    return;

def reshape_SSP_table ( L, shape ) :
    """Reshapes a flattened SSP table to the original 3-dimensional shape.
    
    The 3D matrix has the dimensions of the 3 grids over which it varies:
    wavelength (l_grid), age (t_grid) and metallicity (Z_grid).
    
    Parameters
    ----------
    L : 1d numpy array
        Input flattened SSP table with length l_grid*t_grid*Z_grid
    shape : tuple
        Tuple of 3 integers with the sizes of
        
        1. wavelength grid 
        2. age grid
        3. metallicity grid
    
        The above ordering is mandatory
    
    Returns
    -------
    : 3-d numpy array
        3D copy of the original array with shape (l_grid, t_grid, Z_grid).
    """
    return numpy.transpose(
        L.reshape(
            (shape[2], shape[0], shape[1])
        ),
        axes = ( 1, 2, 0 )
    )

class CSP () :
    """ The Composite Stellar Population object handles and manipulates SSP libraries.

    It provides convenient functions for reading the SSPs from the provided formatted 
    binary files. 
    SSPs can be extracted for given age/wavelength/metallicity.
    Providing an object of type galapy.StarFormationHistory the self.emission function
    computes the luminosity of the resulting Composite Stellar Population.
    
    Parameters
    -----------
    ssp_lib : string
      which SSP library to load. The default is :code:`parsec22.NT`-
      To see the list of available libraries run 
      :code:`galapy.CompositeStellarPopulation.print_ssp_libs()`
    CCSN : bool
      whether to allow for Core-Collapse Supernova support (default=False).

    See Also
    --------
    :func:`galapy.CompositeStellarPopulation.print_ssp_libs()`
    """

    def __init__ ( self, ssp_lib = 'parsec22.NT', CCSN = False ) :
        
        self.ssp_lib = None
        if ssp_lib in set(_SSP_LIB.keys()) :
            self.ssp_lib = ssp_lib
        else :
            raise ValueError( f'SSP library "{ssp_lib}" not available; '
                              'you can see the list of available libraries by running '
                              'galapy.CompositeStellarPopulation.print_ssp_libs()' )
        
        self.l, self.t, self.Z, self.L = load_SSP_table( self.ssp_lib )
        self.shape = (self.l.size, self.t.size, self.Z.size)

        self.CCSN = CCSN
        self.core = CCSP( self.l, self.t, self.Z, self.L, self.CCSN )
        self.L = reshape_SSP_table( self.L, self.shape )
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
          index in the wavelength-grid of pre-computed SSPs
        it : uint or array-like of uints
          index in the time-grid of pre-computed SSPs
        iz : uint or array-like of uints
          index in the metallicity-grid of pre-computed SSPs
        
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
        """ Computes the CSP emission at given index in the wavelength-grid.

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
        size of the wavelength grid (chosen through the keyword argument :code:`il`) and 
        :math:`N_\\tau` is the size of the time grid of the SSP (:code:`self.shape[1]`).
        :code:`ftau` defaults to :code:`None`, in which case no attenuation is applied.

        Parameters
        ----------
        age : float
          The age in years of the CSP.
        sfh : object of type SFH()
          The chosen star formation history. This has to be (or inherit from) 
          the instance of an object of type galapy.StarFormationHistory.SFH().
        il : array of int
           array of indexes of the positions in the wavelength-grid
           for which to compute the emission. The default value is :code:`None`, for which 
           the function will return the emission all over the wavelength grid.
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
        
    
