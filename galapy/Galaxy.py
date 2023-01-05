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

# gxy_tunables = ( 'age', 'redshift' )
# def gxy_build_params () :
#     pass

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
        # _ = {
        #     'lstep'    : lstep,
        #     'cosmo'    : cosmo,
        #     'do_Xray'  : do_Xray,
        #     'do_Radio' : do_Radio,
        #     'do_AGN'   : do_AGN,
        #     'csp'      : csp,
        # }

        # Define a dictionary to store components
        self.components = {
            'stellar' : None,
            'extinct' : None,
            'MC' : None,
            'DD' : None,
        }

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
        # NOTE THAT a more consistent implementation would be
        # to also set the parameter to false if synchrotron
        # is already included in the SSPs, left for future development:
        # csp[ 'CCSN' ] = do_Radio and 'br22.NT' not in csp['ssp_lib']

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
            self.components['XRB'] = None
            # This here below not really necessary as there is
            # no freedom in setting the X-ray Binaries parameters:
            # self.params[ 'xrb' ] = self.xrb.params
            
        self.agn = None
        if do_AGN :
            if agn is None :
                agn = {}
            self.agn = AGN( lmin = self.csp.l.min(),
                            lmax = self.csp.l.max(),
                            do_Xray = do_Xray,
                            **agn )
            self.params[ 'agn' ] = self.agn.params
            self.components['AGN'] = None

        self.nff   = None
        self.snsyn = None
        # Build the Nebular-Free support only if
        # it is not already included in the SSP library
        if do_Radio and 'br22.NTL' not in self.csp.ssp_lib :
            if nff is None :
                nff = {}
            self.nff = NFF( self.csp.l, **nff )
            self.Q_H_fact = 1. / clight['A/s'] / hP['erg*s']
            self.w912 = self.csp.l <= 912.
            self.params[ 'nff' ] = self.nff.params
            self.components['nebular'] = None

            # Build the Synchrotron support only if
            # it is not already included in the SSP library
            if 'br22.NT' not in self.csp.ssp_lib :
                if syn is None :
                    syn = {}
                syn[ 'RCCSN' ] = self.csp.core.RCCSN()
                self.snsyn = SNSYN( self.csp.l, **syn )
                self.params[ 'syn' ] = self.snsyn.params
                self.components['synchrotron'] = None
                
        if lstep is not None :
            self.lgrid = self.get_wavelenght_grid(lstep)
        else : 
            self.lgrid = numpy.arange(self.csp.l.size)
        
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
            self.params['redshift'] = self.redshift
            
        if age is not None :
            # if the provided age is larger than the age of the Universe,
            # set the age to the age of the Universe
            self.age = numpy.min( [ age, self.UA ] )
            reset_csp = True
            reset_ism = True
            self.params['age'] = self.age

        if sfh is not None :
            self.sfh.set_parameters(**sfh)
            reset_csp = True
            reset_ism = True
            self.params['sfh'].update(self.sfh.params)

        if reset_csp :
            self.csp.set_parameters( self.age, self.sfh )
            if self.snsyn is not None :
                if syn is None :
                    syn = {}
                syn['RCCSN'] = self.csp.core.RCCSN()

        ############################################################################
        # ISM:
        # --------------
        # None  | Reset 
        # ------+-------
        # True  | True
        # True  | False
        # False | True
        # False | False
        # ______|_______
        # --------------
        # This should cover all the possibilities:
        if ism is not None :
            if reset_ism :
                ism.update( { 'Zgas'  : self.sfh.core.Zgas(self.age), 
                              'Mgas'  : self.sfh.core.Mgas(self.age),
                              'Mdust' : self.sfh.core.Mdust(self.age) } )
            self.ism.set_parameters(**ism)
            self.params['ism'].update(self.ism.mc.params)
            self.params['ism'].update(self.ism.dd.params)
        else :
            if reset_ism :
                ism = { 'Zgas'  : self.sfh.core.Zgas(self.age), 
                        'Mgas'  : self.sfh.core.Mgas(self.age),
                        'Mdust' : self.sfh.core.Mdust(self.age) }
                self.ism.set_parameters(**ism)
                self.params['ism'].update(self.ism.mc.params)
                self.params['ism'].update(self.ism.dd.params)

        ############################################################################
        # This one here should also be set only when either age or sfh changes:
        if self.xrb is not None :
            self.xrb.set_parameters( age = self.age,
                                     psi = self.sfh( self.age ),
                                     Mstar = self.sfh.core.Mstar( self.age ),
                                     Zstar = self.sfh.core.Zstar( self.age ) )

        ############################################################################
        # AGN varies only if some new parameters are passed
        if agn is not None :
            try :
                self.agn.set_parameters(**agn)
                self.params['agn'].update(self.agn.params)
            except AttributeError :
                raise AttributeError( 'Passing AGN-parameters to a GXY-class built '
                                      'without an AGN component is not allowed.' )

        ############################################################################
        # NFF varies both whether some new parameters are passed or if some
        # of the other components above which they do depend changes
        if self.nff is not None :
            if nff is None :
                nff = {}
            nff.update( { 'Zgas' : self.sfh.core.Zgas(self.age) } )
            self.nff.set_parameters( **nff )
            self.params['nff'].update(self.nff.params)
            # try :
            #     nff.update( { 'Zgas' : ism[ 'Zgas' ] } )
            #     self.nff.set_parameters( **nff )
            #     self.params['nff'].update(self.nff.params)
            # except AttributeError :
            #     raise AttributeError( 'Passing NFF-parameters to a GXY-class built '
            #                           'without an NFF component is not allowed.' )

        ############################################################################
        # If the object exists and CCSN varied for some reason, csp has already
        # taken care of modifying the syn-dictionary, otherwise this enters only
        # when some new parameter has been set.
        if syn is not None:
            try :
                self.snsyn.set_parameters( **syn )
                self.params['syn'].update(self.snsyn.params)
            except AttributeError :
                raise AttributeError( 'Passing SYN-parameters to a GXY-class built '
                                      'without an SNSYN component is not allowed.' )
            
        return;
    
    def Lstellar ( self ) :
        
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
        # This here below should be removed as it is now authomatic with set_parameters()
        # set derived stellar properties
        # self.csp.set_parameters( self.age, self.sfh )
        
        # attenuation from ISM
        attTotMC, attTot = self.ism.total_attenuation( self.wl(), self.csp.t )
        
        # emission from stars (using directly the core function for performance)  
        Lunatt = self.csp.core.emission( self.lgrid )
        LattMC = self.csp.core.emission( self.lgrid, attTotMC )
        Ltot   = self.csp.core.emission( self.lgrid, attTot )
        self.components['stellar'] = Lunatt
        self.components['extinct'] = numpy.array(Ltot)

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
        self.components['MC'] = self.ism.mc.emission( self.wl() )
        self.components['DD'] = self.ism.dd.emission( self.wl() )
        Ltot += self.components['MC'] + self.components['DD']

        # OPTIONAL COMPONENTS:

        # emission from AGN
        if self.agn is not None :
            self.components['AGN'] = self.agn.emission( self.wl(), (EDD+EMC) * sunL )
            Ltot += self.components['AGN']

        # emission from X-Ray binaries
        if self.xrb is not None :
            self.components['XRB'] = self.xrb.emission( self.wl() )
            Ltot += self.components['XRB']

        # Free-free emission from Nebulae
        if self.nff is not None :
            Q_H = trap_int( self.wl()[ self.w912 ],
                            Lunatt[ self.w912[ self.lgrid ] ] *
                            self.wl()[ self.w912 ] * self.Q_H_fact )
            self.components['nebular'] = self.Aavg * self.nff.emission( Q_H, il = self.lgrid )
            Ltot += self.components['nebular']

        # Synchrotron emission from SN-shocks
        if self.snsyn is not None :
            self.components['synchrotron'] = self.snsyn.emission( self.lgrid ) * sunL
            Ltot += self.components['synchrotron']
        
        return Ltot

    def get_avgAtt ( self ) :
        """ Returns the average attenuation in absolute magnitudes 
        """
        if self.components['stellar'] is None or self.components['extinct'] is None :
            raise RuntimeError(
                "Cannot compute average attenuation if the components dictionary "
                "has not been computed. First run function GXY.get_emission()"
            )
        wn0 = self.components['stellar'] > 0
        AA = numpy.zeros_like( self.wl() )
        AA[wn0] = -2.5 * numpy.log10(
            self.components[ 'extinct' ][ wn0 ] / self.components[ 'stellar' ][ wn0 ]
        )
        return AA

    def get_SED ( self ) :
        """ Returns the flux at given distance in units of milli-Jansky [mJy].
        lambda_R^2 * L_tot(lambda_R) * (1+z)/(4 * pi * c * D_L^2 )
        """
        return self.cosmo.to_flux( self.redshift, self.wl(), self.get_emission() )

    def components_to_flux ( self ) :
        return {
            k : self.cosmo.to_flux( self.redshift, self.wl(), c ) 
            for k, c in self.components.items()
        }
    

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

