from abc import ABC, abstractmethod

class Model (ABC) :
    
    def __init__ ( self ) :
        """ Abstract base interface for all the objects that 
        should be considered models.
        """
        self.params = {}
        
    @abstractmethod
    def set_parameters ( self ) :
        pass


