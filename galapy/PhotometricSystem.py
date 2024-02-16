"""Class and tools for building the photometric system and inspecting the filters database

.. Tip::
   The ``PMS`` class is implemented to receive a variable number of arguments 
   (i.e. all the filters necessary, either from the database or custom).
   Given that
   
   - the filters from the database can be given as input to the ``PMS`` class as a sequence 
     of strings, each with the unique name of a filter;
   - the function :py:func:`galapy.PhotometricSystem.list_filters` returns
     a list of filters that can be selected by experiment (i.e. instrument) 
   
   a convenient way for passing all the filters from a selected instrument as arguments
   of the ``PMS`` class builder is to unpack the list returned by the list_filters function,
   e.g.:
   
   .. code:: python
      
      pms = PMS( *list_filters( 'Herschel.PACS' ), 'Herschel.SPIRE.PMW' )

   will build a PMS object with all the PACS filters and the PMW filter from the SPIRE
   experiment.
"""

########################################################################################

# External imports
import numpy
import os
from collections.abc import MutableMapping as MM

# Internal imports
import galapy.internal.globs as GP_GBL
from galapy.internal.data import DataFile
from galapy.internal.utils import FlagVal
from galapy.BandpassTransmission import BPT

########################################################################################

def _recursive_print_filters ( root, pathline = [] ) :
    _, subd, files = next(os.walk(root))
    for f in files : print('.'.join(pathline+f.split('.')[:-1]))
    for s in subd : _recursive_print_filters(os.path.join(root, s), pathline+[s])
    return;    

def print_filters ( experiment = None ) :
    """Prints on screen all filters available in database.
    
    Parameters
    ----------
    experiment : string or None
        (Optional, default = None) filter the results by experiment.
    
    Examples
    --------
    >>> print_filters(experiment='Herschel.SPIRE')
    # Path: /path/to/.galapy
    Herschel.SPIRE.PSW
    Herschel.SPIRE.PLW
    Herschel.SPIRE.PMW
    """
    from galapy.configuration import rcParams
    for path in rcParams['datapath'] :
        print(f"\n# Path: {path}\n")
        if experiment is None or len(experiment) == 0 :
            _recursive_print_filters( root = os.path.join(path, *GP_GBL.FILT_DIR ),
                                      pathline = [] )
        else :
            try :
                exppath = experiment.split('.')
                _recursive_print_filters( root = os.path.join(
                    path, *GP_GBL.FILT_DIR, *experiment.split('.') ),
                                          pathline = exppath )
            except :
                raise AttributeError( f'experiment {experiment} not in database.' )
    return;

def _recursive_list_filters ( root, pathline = [], outlist = [] ) :
    _, subd, files = next(os.walk(root))
    for f in files : outlist += ['.'.join(pathline+f.split('.')[:-1])]
    for s in subd : _recursive_list_filters(os.path.join(root, s), pathline+[s], outlist)
    return;    

def list_filters ( experiment = None ) :
    """Return a list of filters available in database.
    
    Parameters
    ----------
    experiment : string or None
        (Optional, default = None) filter the results by experiment.

    Returns
    -------
    : list
    a list of strings, containing the name of each filter in the database
    
    Examples
    --------
    >>> list_filters(experiment='Herschel.SPIRE')
    ['Herschel.SPIRE.PSW', 'Herschel.SPIRE.PLW', 'Herschel.SPIRE.PMW' ]
    """
    from galapy.configuration import rcParams
    outlist = []
    for path in rcParams['datapath'] :
        if experiment is None or len(experiment) == 0 :
            _recursive_list_filters( root = os.path.join(path, *GP_GBL.FILT_DIR ),
                                     pathline = [], outlist = outlist )
        else :
            try :
                exppath = experiment.split('.')
                _recursive_list_filters(
                    root = os.path.join(
                        path, *GP_GBL.FILT_DIR, *experiment.split('.') ),
                    pathline = exppath, outlist = outlist )
            except :
                raise AttributeError( f'experiment {experiment} not in database.' )
    return outlist

########################################################################################

class PMS () :
    """ Builds the Photometric System 
    
    Parameters
    ----------
    *args : tuple of strings
       As many strings as necessary. Each string has to name
       one of the available filters already present in the database.
    **kwargs : dictionary
       Keyword arguments can be used for providing user-defined 
       transmissions. 
       Custom transmissions are passed as properly formatted dictionaries:
       ``keyword = { 'wavelengths' : array-like,'photons' : array-like }``
       The two arrays provided, as the keys suggest, must define the
       wavelength grid and the corresponding transmission in photon units
       (and thus they must have the same size).
       Note that the chosen keyword will be used as unique identifier of
       the custom transmission.

    Examples
    --------
    >>> from galapy.PhotometricSystem import PMS
    >>> pms = PMS( 'Herschel.PACS.green', 'ALMA.B3' )
    
    builds a photometric system from two filters in the database, namely
    ALMA Band 3 and the green filter from the PACS experiment of the Herschel satellite
    
    >>> pms = PMS(
    ...     custom1 = {
    ...             'wavelengths' : [999., 1000., 1001., 1999., 2000., 2001.],
    ...             'photons' : [0., 1., 1., 1., 1., 0.]
    ...     }
    ... )
    
    builds a photometric system from a custom top-hat filter between 1000 to 2000 Angstrom,
    It is also possible to build heterogeneous objects, with the only caveat of providing
    the names of filters from the database first and then all the custom ones:
    
    >>> pms = PMS(
    ...     'Herschel.PACS.green', 'ALMA.B3',
    ...     custom1 = {
    ...             'wavelengths' : [999., 1000., 1001., 1999., 2000., 2001.],
    ...             'photons' : [0., 1., 1., 1., 1., 0.]
    ...     }
    ... )
    """

    def __init__ ( self, *args, **kwargs ) :

        # Build the Bandpass-Transmissions dictionary:
        self.bpt = {}

        # allocate an array to append filter limits
        extremes = numpy.array([], dtype=object)
        
        # From positional arguments pass band-names
        if len( args ) > 0 :
            for arg in args :
                try :
                    *path, name = arg.split('.')
                    ll, fl = numpy.loadtxt( DataFile(
                        f'{name:s}.dat',
                        numpy.append( GP_GBL.FILT_DIR, path )
                    ).get_file(), unpack=True )
                except OSError :
                    raise ValueError(
                        f'Filter "{arg}" provided is not present in galapy database.'
                    )
                except :
                    raise ValueError(
                        f'Argument "{arg}" provided does not name a known filter.'
                    )

                idx_sort = numpy.argsort( ll )
                ll = numpy.ascontiguousarray(ll[idx_sort])
                fl = numpy.ascontiguousarray(fl[idx_sort])

                try :
                    self.bpt[ arg ] = BPT( ll, fl )
                except :
                    ValueError()

                # append interval extremes to flattened array of flagged values
                # ('s' = start, 'e' = end)
                extremes = numpy.append( FlagVal( self.bpt[arg].get_lmin(), 's' ), extremes )
                extremes = numpy.append( FlagVal( self.bpt[arg].get_lmax(), 'e' ), extremes )

        # From keyword arguments pass custom bands
        for key, value in kwargs.items() :
            if isinstance( value, MM ) :
                try :
                    ll = numpy.asarray( value[ 'wavelengths' ] )
                    fl = numpy.asarray( value[ 'photons' ] )
                    idx_sort = numpy.argsort( ll )
                    ll = numpy.ascontiguousarray(ll[idx_sort])
                    fl = numpy.ascontiguousarray(fl[idx_sort])
                    self.bpt[ key ] = BPT( ll, fl )
                except KeyError :
                    raise ValueError(
                        f'The provided dictionary "{key}" does not contain the necessary items, '
                        'please provide a "wavelengths" and a "photons" array.'
                    )
                except :
                    raise ValueError()
                
                # append interval limits to flattened array of flagged values
                # ('s' = start, 'e' = end)
                extremes = numpy.append( FlagVal( self.bpt[key].get_lmin(), 's' ), extremes )
                extremes = numpy.append( FlagVal( self.bpt[key].get_lmax(), 'e' ), extremes )
                
            else :
                raise ValueError( 'Keyword arguments must be correctly formatted dictionaries: '
                                  'k = {"wavelengths":array-like,"photons":array-like}' )

        # Sort the flattened extremes array
        extremes = extremes[numpy.argsort([ext.value for ext in extremes])]

        # remove overlappings in O(NlogN) time (with N = number of intervals):
        self.olims = []
        count = 0
        for ext in extremes :
            if ext.flag == 's' :
                if count == 0 :
                    self.olims += [[ext.value]]
                count += 1 
            if ext.flag == 'e' :
                count -= 1
                if count == 0 :
                    self.olims[-1] += [ext.value]
        self.olims = numpy.asarray( self.olims )

        # Store array with pivot wavelengths
        self.lpiv = numpy.asarray( [ bpt.get_lpiv()
                                     for bpt
                                     in self.bpt.values() ] )
        self.keys = numpy.asarray( [ k for k in self.bpt.keys() ] )
        
        idxsort = numpy.argsort( self.lpiv )
        self.lpiv = self.lpiv[idxsort]
        self.keys = tuple( self.keys[idxsort] )

    def dump ( self ) :
        """Dump the Photometric System to formatted dictionary."""
        return {
            k : { 'wavelengths' : numpy.asarray( v.get_xaxis() ),
                  'photons'     : (
                      numpy.asarray( v.get_yaxis() ) *
                      numpy.asarray( v.get_xaxis() ) /
                      v.get_norm() ) }
            for k, v in self.bpt.items()
        }

    @classmethod
    def load ( cls, dictionary ) :
        """Load the Photometric System from a formatted dictionary"""
        return cls( **dictionary )

    def __len__ ( self ) :
        """Return number of filters"""
        return len( self.bpt )

    def get_intervals ( self ) :
        # not implemented yet:
        # should return the intervals within which there are filters
        # avoiding overlaps
        pass
    
    def get_fluxes ( self, ll, fl ) :
        """ Get a list of fluxes, one per each band in the photometric system.
        
        Parameters
        ----------
        ll : array-like
           wavelength grid 
        fl : array-like
           fluxes defined on the wavelength grid
        
        Returns
        -------
        : list
        the fluxes from the input array trasmitted by each
        filter in the photometric system
        """
        fluxes = []
        for key in self.keys :
            wl, wu = [ numpy.argmax( self.bpt[key].get_lmin() < ll ),
                       numpy.argmin( ll < self.bpt[key].get_lmax() ) ]
            fluxes += [self.bpt[key].get_bandpass_flux( ll[wl:wu], fl[wl:wu] )]
    
        return fluxes

