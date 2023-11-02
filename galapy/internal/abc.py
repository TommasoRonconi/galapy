from abc import ABC, abstractmethod, abstractclassmethod

class GalaPyObject ( ABC ) :
    
    @abstractmethod
    def _dump_hdf5 ( self, outfile ) :
        pass

    @abstractclassmethod
    def _load_hdf5 ( cls, indict ) :
        pass
    
    def dump ( self, outfile, method = 'hdf5' ) :
        if method == 'hdf5' :
            self._dump_hdf5( outfile )
            return;
        if method == 'pickle' :
            import pickle
            with open( outfile, 'wb' ) as f :
                pickle.dump( self, f )
            return;
        return;

    @classmethod
    def load ( cls, infile, method = 'hdf5' ) :
        if method == 'hdf5' :
            from galapy.io import load_from_hdf5
            return _load_hdf5( load_from_hdf5( infile ) )
        if method == 'pickle' :
            import pickle
            with open( infile, 'rb' ) as f :
                return pickle.load( f )

#class Model ( GalaPyObject ) :
class Model ( ABC ) :
    
    def __init__ ( self ) :
        """ Abstract base interface for all the objects that 
        should be considered models.
        """
        self.params = {}
        
    @abstractmethod
    def set_parameters ( self ) :
        pass
