# External imports
import numpy

# Internal imports
from galapy.StarFormationHistory import SFH
from galapy.CompositeStellarPopulation import CSP
from galapy.InterStellarMedium import ISM
from galapy.internal.utils import trap_int
from galapy.internal.constants import Lsun

class GXY () :
    """ Wrapping everything up
    * It has some own parameters (e.g. the age and redshift)
    * It shall give the possibility of selecting which components to consider
      (even though sfh, csp and maybe ism should be mandatory)
    * Should figure out a smart way of passing arguments to constructor
    """
    
    def __init__ ( self, age, redshift, lstep = None, 
                   sfh_kw = {}, csp_kw = {}, ism_kw = {} ) :
        
        self.age      = age
        self.redshift = redshift
        self.sfh = SFH( **sfh_kw )
        self.csp = CSP( **csp_kw )
        self.ism = ISM( **ism_kw )
        if lstep is not None :
            self.lgrid = self.get_wavelenght_grid(lstep)
        else : 
            self.lgrid = numpy.arange(self.csp.l.size)
        
    def wl ( self ) :
        """ Wavelenght grid with mask applied
        """
        return self.csp.l[self.lgrid]
        
    def get_wavelenght_grid ( self, lstep ) :
        """ Reduces the granularity of the wavelenght grid [optimization]
        
        Parameters
        ----------
        lstep : scalar int or boolean mask
            If the input is scalar, select one any `lstep` wavelenghts.
            If the input is a boolean mask, only select the indexes of the
            array masked with lstep.
            
        Returns
        -------
        : integer array
            A list of indices of the wavelenght grid
        """
        if isinstance(lstep, int) :
            return numpy.arange(self.csp.l.size, step=lstep)
        try :
            return numpy.arange(self.csp.l.size)[lstep]
        except IndexError :
            raise TypeError('Argument lstep should be either an integer or a boolean mask!')

    # def set_wavelenght_grid ( self, lstep = None ) :
    #     """ Reduces the granularity of the wavelenght grid [optimization]
        
    #     Parameters
    #     ----------
    #     lstep : scalar int or boolean mask
    #         If the input is None (default), return a list with all the indices 
    #         in the wavelenght grid.
    #         If the input is scalar, select one any `lstep` wavelenghts.
    #         If the input is a boolean mask, only select the indexes of the
    #         array masked with lstep.
            
    #     Returns
    #     -------
    #     : integer array
    #         A list of indices of the wavelenght grid
    #     """
    #     if lstep is None :
    #         self.lgrid = numpy.arange( self.csp.l.size )
    #         return self.lgrid
    #     if isinstance(lstep, int) :
    #         self.lgrid = numpy.arange(self.csp.l.size, step=lstep)
    #         return self.lgrid
    #     try :
    #         self.lgrid =  numpy.arange(self.csp.l.size)[lstep]
    #         return self.lgrid
    #     except IndexError :
    #         raise TypeError('Argument lstep should be either an integer or a boolean mask!')
    
    def set_parameters ( self, age = None, redshift = None, sfh_kw = None, ism_kw = None ) :
        """divided in nested dictionaries or not?"""
        if age is not None :
            self.age = age
        if redshift is not None :
            self.redshift = redshift
            self._z = 1. / ( 1 + self.redshift )
        ism_d = {}
        if sfh_kw is not None :
            self.sfh.set_parameters(**sfh_kw)
            ism_d = { 'Zgas'  : self.sfh.core.Zgas(self.age), 
                      'Mgas'  : self.sfh.core.Mgas(self.age),
                      'Mdust' : self.sfh.core.Mdust(self.age) }
        if ism_kw is not None :
            ism_kw.update(ism_d)
            self.ism.set_parameters(**ism_kw)
        return;
    
    def Lstellar ( self ) :
        
        # set derived stellar properties
        self.csp.set_parameters( self.age, self.sfh )
        
        # emission from stars (using directly the core function for performance)  
        return self.csp.core.emission( self.lgrid )
    
    def get_emission ( self, **kwargs ) :
        """ Temporary implementation, final should
        * the keyword arguments are the free parameters, 
          (so this should call self.set_parameters(**kwargs) ftf)
        * what about the wavelenghts?
        
        Returns
        -------
        : array-like
            the emission on the selected wavelenght grid 
            in units of solar luminosities (:math:`[L_\odot]`)
        """
        # set derived stellar properties
        self.csp.set_parameters( self.age, self.sfh )
        
        # attenuation from ISM
        attTotMC, attTot = self.ism.total_attenuation( self.wl(), self.csp.t )
        
        # emission from stars (using directly the core function for performance)  
        Lunatt = self.csp.core.emission( self.lgrid )
        LattMC = self.csp.core.emission( self.lgrid, attTotMC )
        Ltot   = self.csp.core.emission( self.lgrid, attTot )
        
        # set the resulting temperature of ISM
        EDD = 0.5 * Lsun * trap_int( 
            self.wl(), ( LattMC - Ltot ) 
        ) 
        EMC = 0.5 * Lsun * trap_int( 
            self.wl(), ( Lunatt - LattMC ) 
        )
        _ = self.ism.dd.temperature( EDD )
        _ = self.ism.mc.temperature( EMC )
         
        # emission from ISM
        Ltot += self.ism.mc.emission( self.wl() )
        Ltot += self.ism.dd.emission( self.wl() )

        return Ltot
