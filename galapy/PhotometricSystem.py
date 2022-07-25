# External imports
import numpy
import os
from collections.abc import MutableMapping as MM

# Internal imports
import galapy.internal.globs as GP_GBL
from galapy.internal.data import DataFile
from galapy.internal.utils import FlagVal
from galapy.BandpassTransmission import BPT

class PMS () :
    """ Builds the Photometric System 
    
    Parameters
    ----------
    *args : strings
       As many strings as necessary. Each string has to name
       one of the available filters already present in the database.
    **kwargs : 'named' dictionaries
       Keyword arguments can be used for providing user-defined 
       transmissions. 
       Custom transmissions are passed as properly formatted dictionaries:
       keyword = { 'wavelenghts' : array-like,
                   'photons' : array-like }
       The two arrays provided, as the keys suggest, must define the
       wavelenght grid and the corresponding transmission in photon units
       (and thus they must have the same size).
       Note that the chosen keyword will be used as unique identifier of
       the custom transmission.
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
                    ll, fl = numpy.loadtxt( DataFile( f'{arg:s}.dat',
                                                      GP_GBL.FILT_DIR ).get_file(),
                                            unpack=True )
                except OSError :
                    raise ValueError( f'Filter "{arg}" provided is not present in galapy database.' )
                except :
                    raise ValueError( f'Argument "{arg}" provided does not name a known filter.' )
                    
                ll = numpy.ascontiguousarray(ll)
                fl = numpy.ascontiguousarray(fl)

                try :
                    self.bpt[ arg ] = BPT( ll, fl )
                except :
                    ValueError()

                # append interval extremes to flattened array of flagged values ('s' = start, 'e' = end)
                extremes = numpy.append( FlagVal( self.bpt[arg].get_lmin(), 's' ), extremes )
                extremes = numpy.append( FlagVal( self.bpt[arg].get_lmax(), 'e' ), extremes )

        # From keyword arguments pass custom bands
        for key, value in kwargs.items() :
            if isinstance( value, MM ) :
                try :
                    ll = numpy.ascontiguousarray(value[ 'wavelenghts' ])
                    fl = numpy.ascontiguousarray(value[ 'photons' ])
                    self.bpt[ key ] = BPT( ll, fl )
                except KeyError :
                    raise ValueError( f'The provided dictionary {k} does not contain the necessary items, '
                                      'please provide a "wavelenghts" and a "photons" array.' )
                except :
                    raise ValueError()
                
                # append interval limits to flattened array of flagged values ('s' = start, 'e' = end)
                extremes = numpy.append( FlagVal( self.bpt[key].get_lmin(), 's' ), extremes )
                extremes = numpy.append( FlagVal( self.bpt[key].get_lmax(), 'e' ), extremes )
                
            else :
                raise ValueError( 'Keyword arguments must be correctly formatted dictionaries: '
                                  'k = {"wavelenghts":array-like,"photons":array-like}' )

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

        # Store array with pivot wavelenghts
        self.lpiv = numpy.asarray( [ bpt.get_lpiv()
                                     for bpt
                                     in self.bpt.values() ] )

    def get_intervals ( self ) :
        """
        """
        pass
    
    def get_fluxes ( self, ll, fl ) :
        """ Get a list of fluxes, one per each band in the photometric system.
        
        Parameters
        ----------
        ll : array-like
        
        fl : array-like
        
        Returns
        -------
        : list
        """
        fluxes = []
        for filt in self.bpt.values() :
            wl, wu = [ numpy.argmax( filt.get_lmin() < ll ),
                       numpy.argmin( ll < filt.get_lmax() ) ]
            fluxes += [filt.get_bandpass_flux( ll[wl:wu], fl[wl:wu] )]
    
        return fluxes
