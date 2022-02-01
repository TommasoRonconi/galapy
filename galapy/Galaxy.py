# External imports
import numpy
import os
from collections.abc import MutableMapping as MM

# Internal imports
from galapy.StarFormationHistory import SFH
from galapy.CompositeStellarPopulation import CSP
from galapy.InterStellarMedium import ISM
from galapy.ActiveGalacticNucleus import AGN
from galapy.PhotometricSystem import PMS
from galapy.internal.utils import trap_int
from galapy.internal.interp import lin_interp
from galapy.internal.constants import Lsun, clight, Mpc_to_cm
import galapy.internal.globs as GP_GBL

DL_DIR = os.path.join( os.path.dirname( GP_GBL.__file__ ),
                       GP_GBL.DL_DIR )

class GXY () :
    """ Wrapping everything up
    * It has some own parameters (e.g. the age and redshift)
    * It shall give the possibility of selecting which components to consider
      (even though sfh, csp and maybe ism should be mandatory)
    * Should figure out a smart way of passing arguments to constructor
    """
    
    def __init__ ( self, age, redshift, lstep = None, cosmo = 'Planck18',
                   sfh_kw = {}, csp_kw = {}, ism_kw = {}, agn_kw = None ) :
        
        self.age      = age
        self.redshift = redshift
        self.sfh = SFH( **sfh_kw )
        ism_kw.update( { 'Zgas'  : self.sfh.core.Zgas(self.age), 
                         'Mgas'  : self.sfh.core.Mgas(self.age),
                         'Mdust' : self.sfh.core.Mdust(self.age) } )
        self.csp = CSP( **csp_kw )
        self.ism = ISM( **ism_kw )
        self.agn = None
        if agn_kw is not None :
            self.agn = AGN( lmin = self.csp.l.min(),
                            lmax = self.csp.l.max(),
                            **agn_kw )
        
        if lstep is not None :
            self.lgrid = self.get_wavelenght_grid(lstep)
        else : 
            self.lgrid = numpy.arange(self.csp.l.size)

        # Build the redshift-dependent constant for lum->flux conversion
        if isinstance( cosmo, str ) :
            zz, DL = numpy.loadtxt( os.path.join( DL_DIR, f'{cosmo:s}.LumDist.txt' ),
                                    unpack = True )
        elif isinstance( cosmo, MM ) :
            zz, DL = cosmo['redshift'], cosmo['luminosity_distance']
        else :
            raise ValueError( 'Argument `cosmo` should be either a string '
                              'or a formatted dictionary.' )
        zz = numpy.ascontiguousarray(zz)
        DL = numpy.ascontiguousarray(
            1.e+26 * (1 + zz) * Lsun /
            ( 4 * numpy.pi * clight['A/s'] * DL**2 * Mpc_to_cm**2 )
        )
        self._fDL = lin_interp( zz, DL )
        
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
    
    def set_parameters ( self, age = None, redshift = None, sfh_kw = None, ism_kw = {}, agn_kw = None ) :
        """divided in nested dictionaries or not?"""
        if age is not None :
            self.age = age
        if redshift is not None :
            self.redshift = redshift
            self._z = 1. / ( 1 + self.redshift )

        if sfh_kw is not None :
            self.sfh.set_parameters(**sfh_kw)
            ism_kw.update( { 'Zgas'  : self.sfh.core.Zgas(self.age), 
                             'Mgas'  : self.sfh.core.Mgas(self.age),
                             'Mdust' : self.sfh.core.Mdust(self.age) } )
        if len( ism_kw ) > 0 :
            self.ism.set_parameters(**ism_kw)
        # if age is not None or sfh_kw is not None :
        #     self.csp.set_parameters( self.age, self.sfh )

        if agn_kw is not None :
            try :
                self.agn.set_parameters(**agn_kw)
            except AttributeError :
                raise AttributeError( 'Passing AGN-parameters to a GXY-class built '
                                      'without an AGN component is not allowed.' )
            
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
        EDD = Lsun * trap_int( 
            self.wl(), ( LattMC - Ltot ) #* self.wl()
        ) 
        EMC = Lsun * trap_int( 
            self.wl(), ( Lunatt - LattMC ) #* self.wl()
        )
        _ = self.ism.dd.temperature( EDD )
        _ = self.ism.mc.temperature( EMC )
         
        # emission from ISM
        Ltot += self.ism.mc.emission( self.wl() )
        Ltot += self.ism.dd.emission( self.wl() )
        
        if self.agn is not None :
            Ltot += self.agn.emission( self.wl(), (EDD+EMC)/Lsun )
        
        return Ltot

    def get_SED ( self ) :
        """ Returns the flux at given distance in units of milli-Jansky [mJy].
        lambda_R^2 * L_tot(lambda_R) * (1+z)/(4 * pi * c * D_L^2 )
        """
        return self.get_emission() * self.wl()**2 * self._fDL( self.redshift )
    

class PhotoGXY ( GXY ) :

    def __init__ ( self, age, redshift, lstep = None, cosmo = 'Planck18', 
                   sfh_kw = {}, csp_kw = {}, ism_kw = {} ) :
        super().__init__( age, redshift,
                          lstep=lstep, cosmo=cosmo,
                          sfh_kw=sfh_kw,
                          csp_kw=csp_kw,
                          ism_kw=ism_kw )
        self.pms = None

    def build_photometric_system ( self, *args, **kwargs ) :

        self.pms = PMS( *args, **kwargs)
        return;

    def photoSED ( self ) :

        if self.pms is None :
            raise Exception( "Photometric system has not been set. "
                             "Use function build_photometric_system() before." )
        return self.pms.get_fluxes( self.wl() * (1+self.redshift), self.get_SED() )
        


class SpectralGXY ( GXY ) :
    pass

