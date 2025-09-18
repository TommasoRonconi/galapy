########################################################################################

import numpy

########################################################################################

from galapy.internal.constants import sunL
from galapy import AGN_core
# from galapy.internal.abc import Model

########################################################################################

class AGN () :

    def __init__ (
            self, lgrid, model = 'Fritz2006',
            **kwargs
    ) :

        if model in { 'Fritz2006', 'Panchromatic' } :
            self.model = model
        else :
            raise RuntimeError(
                f'Requested model "{model}" not available. ' +
                'Available models are "Fritz2006" and "Panchromatic"'
            )
        try :
            self.core = getattr( AGN_core, self.model )( lgrid, **kwargs )
        except :
            raise
        self.params = dict( self.core.params ) # deep copy
        self.params['model'] = self.model

    def set_parameters ( self, **kwargs ) :
        self.params.update( **kwargs )
        return self.core.set_parameters( **kwargs )

    def emission ( self, *args, **kwargs ) :
        return self.core.emission( *args, **kwargs ) * sunL

    def __call__ ( self, *args, **kwargs ) :
        return self.core.__call__( *args, **kwargs )
        

########################################################################################
