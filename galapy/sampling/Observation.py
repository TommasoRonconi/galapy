import numpy
from collections.abc import MutableMapping, Sequence
from galapy.PhotometricSystem import PMS

class Observation () :
    
    def __init__ ( self, bands, fluxes, errors, uplims, filters ) :
        """ Class for managing photometric observations of galaxies.
        It builds the photometric system from arguments and stores the 
        observational data in ordered numpy arrays. The order is set
        by the photometric system chosen. By construction, it follows
        the pivot wavelenghts of the transmission filters provided.
        
        Parameters
        ----------
        bands : array-like
          sequence of strings containing the name of the bandpass 
          filters associated to the obseved dataset

        fluxes : array-like
          sequence of fluxes (or upper limits) measured
          in each of the bandpass filters listed in argument `bands`

        errors : array-like
          sequence of errors associated to the fluxes (or upper limits)
          listed in argument `fluxes`

        uplims : array-like
          sequence of booleans identifying whether the values listed 
          in argument `fluxes` should be considered non-detection (`True`)
          or a detection (`False`)

        filters : either sequence or dictionary
          the argument to be passed to an object of type 
          `galapy.PhotometricSystem.PMS`, this can either be a 
          sequence of strings naming bandpass filters already present in
          the database, or a nested dictionary with formatted key-value pairs
          (for further details see the documentation of `galapy.PhotometricSystem.PMS`)
        """
    
        if isinstance( filters, list ) or isinstance( filters, tuple ) :
            self.pms = PMS( *filters )
        elif isinstance( filters, MutableMapping ) :
            self.pms = PMS( **filters )
        else :
            raise RuntimeError( 'filters argument must be either a list (or a tuple) of strings '
                                'or a properly formatted dictionary.' )
            
        if any( [ len( lst ) != len( self.pms ) 
                  for lst in [ bands, fluxes, errors, uplims ] ] ) :
            raise ValueError( 'All arguments should have the same lenght.' )
        
        # create temporary data dictionary to iterate by-key in the order provided by the pms object
        data = { k : v for k, v in zip( bands, zip( fluxes, errors, uplims ) ) }
        
        # store internally the dataset with the same order
        # of the bands defined in the photometric system
        self.fluxes = numpy.asarray( [ data[k][0] for k in self.pms.keys ] )
        self.errors = numpy.asarray( [ data[k][1] for k in self.pms.keys ] )
        self.uplims = numpy.asarray( [ data[k][2] for k in self.pms.keys ] )
