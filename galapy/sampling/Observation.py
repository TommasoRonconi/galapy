import numpy
from collections.abc import MutableMapping, Sequence
from galapy.PhotometricSystem import PMS

class Observation () :
    
    def __init__ ( self, bands, fluxes, errors, uplims, pms ) :
        """ Class for managing photometric observations of galaxies.
        It builds the photometric system from arguments and stores the 
        observational data in ordered numpy arrays. The order is set
        by the photometric system chosen. By construction, it follows
        the pivot wavelengths of the transmission filters provided.
        
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

        pms : galapy.PhotometricSystem.PMS
          instance of type ``galapy.PhotometricSystem.PMS`` with 
          the collection of bandpass transmissions producing the 
          integrated fluxes for the observation (same as listed in
          argument ``bands``)
        """

        if not isinstance( pms, PMS ) :
            try :
                self.pms = PMS( *bands )
            except Exception as err :
                print(
                    'Not possible to deduce the photometric system '
                    'from the list of bands provided. Exited with exception: '
                    f'{err:s}' )
                raise
            else :
                raise RuntimeError(
                    'pms argument must be an instance of type ``galapy.PhotometricSystem.PMS``'
                )
        self.pms = pms
            
        if any( [ len( lst ) != len( self.pms ) 
                  for lst in [ bands, fluxes, errors, uplims ] ] ) :
            raise ValueError( 'All arguments should have the same length.' )
        
        # create temporary data dictionary to iterate by-key in the order provided by the pms object
        data = { k : v for k, v in zip( bands, zip( fluxes, errors, uplims ) ) }
        
        # store internally the dataset with the same order
        # of the bands defined in the photometric system
        self.fluxes = numpy.asarray( [ data[k][0] for k in self.pms.keys ] )
        self.errors = numpy.asarray( [ data[k][1] for k in self.pms.keys ] )
        self.uplims = numpy.asarray( [ data[k][2] for k in self.pms.keys ] )

    def dump ( self ) :
        return dict(
            bands  = '|'.join(self.pms.keys),
            fluxes = self.fluxes,
            errors = self.errors,
            uplims = self.uplims,
            pms_kwargs = self.pms.dump(),
        )

    @classmethod
    def load ( cls, dictionary ) :
        # ensure deep-copy
        dictionary = dict( dictionary )
        return cls(  
            dictionary.pop('bands').split('|'),
            pms = PMS.load( dictionary.pop( 'pms_kwargs' ) ),
            **dictionary
        )
