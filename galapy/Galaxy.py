""" Galaxy module
"""

# External imports
import numpy
import os
from collections.abc import MutableMapping as MM

# Internal imports
from galapy.StarFormationHistory import SFH
from galapy.CompositeStellarPopulation import CSP
from galapy.InterStellarMedium import ISM
from galapy.ActiveGalacticNucleus import AGN
from galapy.XRayBinaries import XRB
from galapy.NebularFreeFree import NFF
from galapy.PhotometricSystem import PMS
from galapy.Synchrotron import SNSYN
from galapy.Cosmology import CSM
from galapy.internal.utils import trap_int
from galapy.internal.interp import lin_interp
from galapy.internal.constants import Lsun, sunL, clight, Mpc_to_cm, hP
import galapy.internal.globs as GP_GBL
from galapy.internal.data import DataFile

class GXY () :
    """ Wrapping everything up
    * It has some own parameters (e.g. the age and redshift)
    * It shall give the possibility of selecting which components to consider
      (even though sfh, csp and maybe ism should be mandatory)
    * Should figure out a smart way of passing arguments to constructor
    """
    
    def __init__ ( self, age = None, redshift = None,
                   cosmo = 'Planck18', lstep = None,
                   do_Xray = False, do_Radio = False, do_AGN = False,
                   sfh = None, csp = None, ism = None,
                   agn = None, nff = None, syn = None ) :

        # Define a dictionary for parameters storage
        self.params = {}

        if age is None :
            self.age = 1.e+6
        else :
            self.age = age
        self.params[ 'age' ] = self.age

        if redshift is None :
            self.redshift = 0.
        else :
            self.redshift = redshift
        self.params[ 'redshift' ] = self.redshift

        # Build the redshift-dependent cosmological quantities:
        self.cosmo = CSM( cosmo )

        # Compute Universe Age at given redshift
        self.UA = self.cosmo.age( self.redshift )

        # Check that age and redshift are compatible
        if self.age > self.UA :
            raise RuntimeError( f"The age of the galaxy (t={self.age:e} years) "
                                f"cannot be larger than the age of the Universe which, "
                                f"at the given redshift ( z = {self.redshift:.3f} ),"
                                f" is {self.UA:e} years." )

        if sfh is None :
            sfh = {}
        if csp is None :
            csp = {}
        if ism is None :
            ism = {}
            
        self.sfh = SFH( **sfh )
        self.params[ 'sfh' ] = self.sfh.params
        ism.update( { 'Zgas'  : self.sfh.core.Zgas(self.age), 
                      'Mgas'  : self.sfh.core.Mgas(self.age),
                      'Mdust' : self.sfh.core.Mdust(self.age) } )
        
        # tell the CSP constructor to build with
        # SN support if Radio support is requested
        csp[ 'CCSN' ] = do_Radio 

        self.csp = CSP( **csp )
        self.csp.set_parameters( self.age, self.sfh )

        self.ism = ISM( **ism )
        self.params[ 'ism' ] = {}
        self.params[ 'ism' ].update( self.ism.mc.params )
        self.params[ 'ism' ].update( self.ism.dd.params )
        
        self.xrb = None
        if do_Xray :
            self.xrb = XRB( lmin = self.csp.l.min(),
                            lmax = self.csp.l.max(),
                            age = self.age,
                            psi = self.sfh( self.age ),
                            Mstar = self.sfh.Mstar( self.age ),
                            Zstar = self.sfh.Zstar( self.age ) )
            self.params[ 'xrb' ] = self.xrb.params
            
        self.agn = None
        if do_AGN :
            if agn is None :
                agn = {}
            self.agn = AGN( lmin = self.csp.l.min(),
                            lmax = self.csp.l.max(),
                            Xray = Xray,
                            **agn )
            self.params[ 'agn' ] = self.agn.params

        self.nff   = None
        self.snsyn = None
        if do_Radio :
            # Build the Nebular-Free support only if
            # it is not already included in the SSP library
            if 'br22.NT' not in self.csp.ssp_lib :
                if nff is None :
                    nff = {}
                self.nff = NFF( self.csp.l, **nff )
                self.Q_H_fact = 1. / clight['A/s'] / hP['erg*s']
                self.w912 = self.csp.l <= 912.
                self.params[ 'nff' ] = self.nff.params

            if syn is None :
                syn = {}
            syn[ 'RCCSN' ] = self.csp.core.RCCSN()
            self.snsyn = SNSYN( self.csp.l, **syn )
            self.params[ 'syn' ] = self.snsyn.params
                
        if lstep is not None :
            self.lgrid = self.get_wavelenght_grid(lstep)
        else : 
            self.lgrid = numpy.arange(self.csp.l.size)

        # Build the redshift-dependent constant for lum->flux conversion
        zz = numpy.ascontiguousarray(self.cosmo.DL.get_x())
        TF = numpy.ascontiguousarray(
            1.e+26 * (1 + zz) * Lsun /
            ( 4 * numpy.pi * clight['A/s'] * self.cosmo.DL.get_y()**2 * Mpc_to_cm**2 )
        )
        self._to_flux = lin_interp( zz, TF )
        
    def wl ( self, obs = False ) :
        """ Wavelenght grid with mask applied
        """
        if obs :
            return ( 1 + self.redshift ) * self.csp.l[ self.lgrid ]
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
        
    def set_parameters ( self, age = None, redshift = None,
                         sfh = None, ism = None,
                         agn = None, nff = None,
                         syn = None ) :
        """divided in nested dictionaries or not?"""

        # if age of sfh change, reset the csp's timetuple
        reset_csp = False
        reset_ism = False
        
        if redshift is not None :
            self.redshift = redshift
            self._z = 1. / ( 1 + self.redshift )
            self.UA = self.cosmo.age( self.redshift )
            
        if age is not None :
            # if the provided age is larger than the age of the Universe,
            # set the age to the age of the Universe
            self.age = numpy.min( [ age, self.UA ] )
            reset_csp = True
            reset_ism = True

        if sfh is not None :
            self.sfh.set_parameters(**sfh)
            reset_csp = True
            reset_ism = True

        if reset_csp :
            self.csp.set_parameters( self.age, self.sfh )
            if self.snsyn is not None :
                rccsn = self.csp.core.RCCSN()
                if syn is None :
                    syn = {}
                    syn['RCCSN'] = rccsn

        if ism is not None :
            if reset_ism :
                ism.update( { 'Zgas'  : self.sfh.core.Zgas(self.age), 
                              'Mgas'  : self.sfh.core.Mgas(self.age),
                              'Mdust' : self.sfh.core.Mdust(self.age) } )
                self.ism.set_parameters(**ism)
        elif reset_ism :
            ism = { 'Zgas'  : self.sfh.core.Zgas(self.age), 
                    'Mgas'  : self.sfh.core.Mgas(self.age),
                    'Mdust' : self.sfh.core.Mdust(self.age) }
            self.ism.set_parameters(**ism)

        # This one here should also be set only when either age or sfh changes:
        if self.xrb is not None :
            self.xrb.set_parameters( age = self.age,
                                     psi = self.sfh( self.age ),
                                     Mstar = self.sfh.Mstar( self.age ),
                                     Zstar = self.sfh.Zstar( self.age ) )

        if agn is not None :
            try :
                self.agn.set_parameters(**agn)
            except AttributeError :
                raise AttributeError( 'Passing AGN-parameters to a GXY-class built '
                                      'without an AGN component is not allowed.' )

        if nff is not None :
            try :
                nff.update( { 'Zgas' : ism[ 'Zgas' ] } )
                self.nff.set_parameters( **nff )
            except AttributeError :
                raise AttributeError( 'Passing NFF-parameters to a GXY-class built '
                                      'without an NFF component is not allowed.' )

        if syn is not None:
            try :
                self.snsyn.set_parameters( **syn )
            except AttributeError :
                raise AttributeError( 'Passing SYN-parameters to a GXY-class built '
                                      'without an SNSYN component is not allowed.' )
            
        return;
    
    def Lstellar ( self ) :
        
        # set derived stellar properties
        # self.csp.set_parameters( self.age, self.sfh )
        
        # emission from stars (using directly the core function for performance)  
        return self.csp.core.emission( self.lgrid )
    
    def get_emission ( self, store_attenuation = False, **kwargs ) :
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
        # self.csp.set_parameters( self.age, self.sfh )
        
        # attenuation from ISM
        attTotMC, attTot = self.ism.total_attenuation( self.wl(), self.csp.t )
        
        # emission from stars (using directly the core function for performance)  
        Lunatt = self.csp.core.emission( self.lgrid )
        LattMC = self.csp.core.emission( self.lgrid, attTotMC )
        Ltot   = self.csp.core.emission( self.lgrid, attTot )

        if store_attenuation or self.nff is not None :
            wn0 = Lunatt > 0.
            self.Aavg = numpy.ones_like( self.wl() )
            self.Aavg[wn0] = Ltot[wn0]/Lunatt[wn0]
        
        # set the resulting temperature of ISM
        EDD = Lsun * trap_int( 
            self.wl(), ( LattMC - Ltot )
        ) 
        EMC = Lsun * trap_int( 
            self.wl(), ( Lunatt - LattMC )
        )
        _ = self.ism.dd.temperature( EDD )
        _ = self.ism.mc.temperature( EMC )
         
        # emission from ISM
        Ltot += self.ism.mc.emission( self.wl() )
        Ltot += self.ism.dd.emission( self.wl() )

        # OPTIONAL COMPONENTS:

        # emission from AGN
        if self.agn is not None :
            Ltot += self.agn.emission( self.wl(), (EDD+EMC) * sunL )

        # emission from X-Ray binaries
        if self.xrb is not None :
            Ltot += self.xrb.emission( self.wl() )

        # Free-free emission from Nebulae
        if self.nff is not None :
            Q_H = trap_int( self.wl()[ self.w912 ],
                            Lunatt[ self.w912[ self.lgrid ] ] *
                            self.wl()[ self.w912 ] * self.Q_H_fact )
            Ltot += self.Aavg * self.nff.emission( Q_H, il = self.lgrid )

        # Synchrotron emission from SN-shocks
        if self.snsyn is not None :
            Ltot += self.snsyn.emission( self.lgrid ) * sunL
        
        return Ltot

    def get_SED ( self ) :
        """ Returns the flux at given distance in units of milli-Jansky [mJy].
        lambda_R^2 * L_tot(lambda_R) * (1+z)/(4 * pi * c * D_L^2 )
        """
        return self.get_emission() * self.wl()**2 * self._to_flux( self.redshift )
    

class PhotoGXY ( GXY ) :

    def __init__ ( self, *args, pms = None, **kwargs ) :

        super().__init__( *args, **kwargs )
        self.pms = pms

    def build_photometric_system ( self, *args, **kwargs ) :

        self.pms = PMS( *args, **kwargs)
        return;

    def photoSED ( self ) :

        if self.pms is None :
            raise Exception( "Photometric system has not been set. "
                             "Call function build_photometric_system() before." )
        return self.pms.get_fluxes( self.wl( obs = True ), self.get_SED() )
        


class SpectralGXY ( GXY ) :
    pass

