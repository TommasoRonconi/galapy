""" Noise module
"""

# External imports
import numpy
from abc import abstractmethod
from collections.abc import MutableMapping as MM

# Internal imports
from galapy.internal.abc import Model

########################################################################################

noise_params_defaults = {

    # Calibration error
    'f_cal' : ( 'Fractional error done in calibrating data',
                [-10., 1], True, 'f_\\mathrm{cal}' )

}

########################################################################################
# Base abstract class

class Noise ( Model ) :
    
    def __init__ ( self ) :
        super().__init__()

    @abstractmethod
    def dump ( self ) :
        pass
        
    @abstractmethod
    def apply ( self ) :
        pass

########################################################################################
# Calibration error

class CalibrationError ( Noise ) :
    
    def __init__ ( self, f_cal = 0.0 ) :
        """
        """
        super().__init__()
        self.f_cal = f_cal
        self.params['f_cal'] = self.f_cal

    def dump ( self ) :
        return self.params

    @classmethod
    def load ( cls, dictionary ) :
        return cls( **dictionary )
        
    def set_parameters ( self, f_cal = None ) :
        """
        """
        if f_cal is not None :
            self.f_cal = f_cal
            self.params['f_cal'] = self.f_cal
        return;
        
    def apply ( self, error, model ) :
        """
        """
        error = numpy.asarray( error )
        model = numpy.asarray( model )
        if len(error) != len(model) :
            raise AttributeError( 'Attributes error and model should have same length' )
        
        return numpy.sqrt( error**2 + ( model * self.f_cal )**2 )

########################################################################################

