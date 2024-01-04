"""
Implements galaxy classes wrapping up the interplay between the different components 
contribution to the overall emission.
By instantiating an object of class ``GXY`` or derived it is possible to set-up the
different components, compute global and component-specific emission, compute the 
derived quantities and modify the parameter setting.
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
from galapy.InterGalacticMedium import IGM

from galapy.internal.utils import trap_int, find_nearest
from galapy.internal.interp import lin_interp
from galapy.internal.constants import Lsun, sunL, clight, Mpc_to_cm, hP
import galapy.internal.globs as GP_GBL
from galapy.internal.data import DataFile
from galapy.internal.abc import Model

########################################################################################

gxy_params_defaults = {

    ##########
    # Galaxy #
    ##########

    'age'      : ( 'Age of the galaxy',
                   [6., 11.], True, 'age' ),
    'redshift' : ( 'Redshift of the galaxy',
                   [0., 10.], False, 'redshift' ),

    ##########################
    # Star Formation History #
    ##########################
    
    'sfh.tau_quench' : ( 'Age of the abrupt quenching',
                         [6., 11.], True, '\\tau_\\mathrm{quench}' ),

    # In-Situ model
    'sfh.psi_max' : ( 'Normalisation',
                      [0., 4.], True, '\\psi_\\mathrm{max}' ),
    'sfh.tau_star' : ( 'Characteristic timescale',
                       [6., 11.], True, '\\tau_\\star' ),

    # Constant model
    'sfh.psi'   : ( 'Value of the constant SFR',
                    [0., 4.], True, '\\psi_0' ),
    'sfh.Mdust' : ( 'Total dust mass in galaxy at the given age',
                    [6., 14.], True, 'M_\\mathrm{dust}' ),
    'sfh.Zgxy'   : ( 'Metallicity of all phases in galaxy at the given age',
                     [0., 1.], False, 'Z_\\mathrm{gxy}' ),
    
    # Delayed-Exponential model
    'sfh.psi_norm' : ( 'Normalisation',
                       [0., 4.], True, '\\psi_\\mathrm{norm}' ),
    'sfh.k_shape'  : ( 'Shape parameter of the early evolution',
                       [0., 5.], False, '\\kappa' ),
    'sfh.tau_star' : ( 'Characteristic timescale',
                       [6., 11.], True, '\\tau_\\star' ),
    'sfh.Mdust' : ( 'Total dust mass in galaxy at the given age',
                    [6., 14.], True, 'M_\\mathrm{dust}' ),
    'sfh.Zgxy'   : ( 'Metallicity of all phases in galaxy at the given age',
                     [0., 1.], False, 'Z_\\mathrm{gxy}' ),
    
    # Log-Normal model
    'sfh.psi_norm'   : ( 'Normalisation',
                         [0., 4.], True, '\\psi_\\mathrm{norm}' ),
    'sfh.sigma_star' : ( 'Characteristic width',
                         [0., 5.], False, '\\sigma_\\star' ),
    'sfh.tau_star'   : ( 'Peak age',
                         [6., 11.], True, '\\tau_\\star' ),
    'sfh.Mdust' : ( 'Total dust mass in galaxy at the given age',
                    [6., 14.], True, 'M_\\mathrm{dust}' ),
    'sfh.Zgxy'   : ( 'Metallicity of all phases in galaxy at the given age',
                     [0., 1.], False, 'Z_\\mathrm{gxy}' ),

    ########################
    # Inter-Stellar Medium #
    ########################

    'ism.f_MC' : ( 'Fraction of dust in the MC phase',
                   [0., 1.], False, 'f_\\mathrm{MC}' ),

    ##################
    # Molecular Clouds
    
    'ism.norm_MC' : ( 'Normalisation of the MC extinction in the visible band',
                      [-1., 4.], True, '\\mathcal{C}_{V~\\mathrm{MC}}' ),
    'ism.N_MC'    : ( 'Number of MCs in the galaxy',
                      [0., 5.], True, 'N_\\mathrm{MC}' ),
    'ism.R_MC'    : ( 'Average radius of a MC',
                      [0., 5.], True, 'R_\\mathrm{MC}' ),
    'ism.tau_esc' : ( 'Time required by stars to start escaping their MC',
                      [4., 8.], True, '\\tau_\\mathrm{esc}' ),
    'ism.dMClow'  : ( 'Molecular clouds extinction index at wavelength < 100 mum (10^6 Ang)',
                      [0., 5.], False, '\\delta_{\\mathrm{MC}~l}' ),
    'ism.dMCupp'  : ( 'Molecular clouds extinction index at wavelength > 100 mum (10^6 Ang)',
                      [0., 5.], False, '\\delta_{\\mathrm{MC}~u}' ),
    
    ##############
    # Diffuse Dust
    
    'ism.norm_DD' : ( 'Normalisation of the DD extinction in the visible band',
                      [-1., 4.], True, '\\mathcal{C}_{V~\\mathrm{DD}}' ),
    'ism.Rdust'   : ( 'Radius of the diffuse dust region embedding stars and MCs',
                      [0., 5.], True, 'R_\\mathrm{DD}' ),
    'ism.f_PAH'   : ( 'Fraction of the total DD luminosity radiated by PAH',
                      [0., 1.], False, 'f_\\mathrm{PAH}' ),
    'ism.dDDlow'  : ( 'Diffuse dust extinction index at wavelength < 100 mum (10^6 Ang)',
                      [0., 5.], False, '\\delta_{\\mathrm{DD}~l}' ),
    'ism.dDDupp'  : ( 'Diffuse dust extinction index at wavelength > 100 mum (10^6 Ang)',
                      [0., 5.], False, '\\delta_{\\mathrm{DD}~u}' ),

    ###############
    # Synchrotorn #
    ###############

    'syn.alpha_syn'   : ( 'SN synchrotron spectral index',
                          [0., 5.], False, '\\alpha_\\mathrm{syn}' ), 
    'syn.nu_self_syn' : ( 'Self-absorption frequency of the SN synchrotron',
                          [0., 1.], False, '\\nu_\\mathrm{syn}^\\mathrm{self}'),

    #####################
    # Nebular Free-Free #
    #####################

    'nff.Zgas' : ( 'Metallicity of the ionised gas of nebular regions',
                   [0., 1.], False, 'Z_\\mathrm{gas}' ),
    'nff.Zi'   : ( 'Average atomic number of ions', [0., 10.], False, '\\mathcal{Z}^i'),

    ###########################
    # Active Galactic Nucleus #
    ###########################

    'agn.fAGN' : ( 'Fraction with respect to the total dust IR luminosity contributed by the AGN',
                   [-3., 3.], True, 'f_\\mathrm{AGN}' ),
    
    'agn.ct' : ( 'Torus half-aperture angle', 40, None, '\\Theta' ),
    'agn.al' : ( 'Density parameter (exponential part)', 0., None, '\\alpha' ),
    'agn.be' : ( 'Density parameter (power-law part)', -0.5, None, '\\beta' ),
    'agn.ta' : ( 'Optical depth at 9.7 mum', 6., None, '\\tau_{9.7}^\\mathrm{AGN}' ),
    'agn.rm' : ( 'Radial ratio of the torus', 60, None, 'R_\\mathrm{torus}^\\mathrm{AGN}' ),
    'agn.ia' : ( 'Inclination angle', 0.001, None, '\\Psi_\\mathrm{los}^\\mathrm{AGN}' ),

}

########################################################################################

class GXY ( Model ) :
    """Base galaxy class. All other galaxy models derive from this object.
    All the arguments are optional.
    Whether a particular component is included or not (e.g. radio, x-ray) the wavelength 
    grid will always span from 1 to :math:`10^{10}` angstroms.
    Note that by default, the only components that will be built are 

    1. star formation history (SFH): modelling the stellar mass growth     
    2. composite stellar populations (CSP): modelling unattenuated stellar emission
    3. inter-stellar medium (ISM): modelling the absorption/re-radiation due to 
       inter-stellar dust. This is divided in two components: molecular clouds (MC)
       and diffuse dust (DD)

    All the other components are not strictly necessay to generate an SED and, for the
    sake of performances, they should be built only if necessary. 
    
    Parameters
    ----------
    age : float
        (``None`` defaults to ``1.e+6``) age in years of the galaxy at the time it is observed
    redshift : float
        (``None`` defaults to ``0.0``) redshift of the galaxy
    cosmo : str or dict
        (default = ``'Planck18'``) cosmological model of choice, takes the same arguments 
        an object of type galapy.Cosmology.CSM would take. 
        If a string is passed, it should name one of the pre-computed cosmologies available in the
        database. Available cosmologies are 
        ``'WMAP7'``, ``'WMAP9'``, ``'Planck15'``, ``'Planck18'``.
        If a dictionary is passed, the class expects to find 3 key-value couples:
        key = ``'redshift'``, value = an array of redshift values;
        key = ``'luminosity_distance'``, value = an array of luminosity distances corresponding to 
        the redshift values of the first key-value couple;
        key = ``'age'``, value = an array of ages of the Universe corresponding to the redshift 
        values of the first key-value couple.
        All these arrays should have the same length.
    lstep :  scalar int or boolean mask
        (default = ``None``) Reduces the granularity of the wavelength grid.
        If the input is scalar, select one any `lstep` wavelengths.
        If the input is a boolean mask, only select the indexes of the
        array masked with lstep.
    do_Xray : bool
        (default = ``False``) build the components modelling X-Ray emission, 
        i.e. X-Ray binaries (low & high mass) and, if ``do_AGN=True``, the X-Ray
        part of the spectrum from the eventual AGN.
    do_Radio : bool
        (default = ``False``) build the components modelling Radio emission if necessary.
        If the chosen SSP library does not already include nebular and synchrotron emission,
        ``do_Radio = True`` will build objects of type NFF (nebular free-free) and SNSYN 
        (super-nova synchrotron). The eventual parameters of both components can be passed
        through arguments ``nff`` and ``syn``, respectively.
        Note that if ``do_Radio = False`` both ``nff`` and ``syn`` will be ignored.
    do_AGN : bool
        (default = ``False``) build AGN component. The AGN spectrum is modelled through 
        templates from Fritz et al. (2006) in the UV to IR bands and with an attenuated
        power-law in the X-ray band (only if also argument ``do_Xray`` is set to ``True``.
        Tunable parameters can be set by passing them to the argument ``agn``.
    do_IGM : bool
        (default = ``False``) build attenuation due to inter-galactic hydrogen with
        model from Inoue et al., (2014).
    sfh : dict
        (default = ``None``) arguments passed to the SFH object builder, if ``None`` is passed
        uses the default parameter of class ``galapy.StarFormationHistory.SFH``
    csp : dict
        (default = ``None``) arguments passed to the CSP object builder, if ``None`` is passed
        uses the default parameter of class ``galapy.CompositeStellarPopulation.CSP``
    ism : dict
        (default = ``None``) arguments passed to the ISM object builder, if ``None`` is passed
        uses the default parameter of class ``galapy.InterStellarMedium.ISM``
    agn : dict
        (default = ``None``) arguments passed to the AGN object builder, if ``None`` is passed
        uses the default parameter of class ``galapy.ActiveGalacticNucleus.AGN``. 
        If ``do_AGN=False`` this argument is ignored.
    nff : dict
        (default = ``None``) arguments passed to the NFF object builder, if ``None`` is passed
        uses the default parameter of class ``galapy.NebularFreeFree.NFF``. 
        If ``do_Radio=False`` this argument is ignored.
    syn : dict
        (default = ``None``) arguments passed to the SNSYN object builder, if ``None`` is passed
        uses the default parameter of class ``galapy.Synchrotron.SNSYN``. 
        If ``do_Radio=False`` this argument is ignored.
    """
    
    def __init__ ( self, age = None, redshift = None,
                   cosmo = 'Planck18', lstep = None,
                   do_Xray = False, do_Radio = False,
                   do_AGN = False, do_IGM = False,
                   sfh = None, csp = None, ism = None,
                   agn = None, nff = None, syn = None ) :

        # Define a dictionary for parameters storage
        self.params = {}

        # Define a dictionary to store components
        self.components = {
            'stellar' : None,
            'extinct' : None,
            'MC' : None,
            'DD' : None,
        }

        # If no age is passed set it to the minimum
        # possible age allowed by the time-integration grid
        if age is None :
            self.age = 1.e+6
        else :
            self.age = age
        self.params[ 'age' ] = self.age

        # If no redshift is passed set it to zero
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

        # Instantiate empty dictionaries for the
        # mandatory components when None is passed.
        if sfh is None :
            sfh = {}
        if csp is None :
            csp = {}
        if ism is None :
            ism = {}

        # Build the SFH model and set the parameters of ISM
        # to the values derived from SFH
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
        # csp[ 'CCSN' ] = do_Radio and 'parsec22.NT' not in csp['ssp_lib']

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
        if do_Radio and 'parsec22.NTL' not in self.csp.ssp_lib :
            if nff is None :
                nff = {}
            self.nff = NFF( self.csp.l, **nff )
            self.Q_H_fact = 1. / clight['A/s'] / hP['erg*s']
            self.w912 = self.csp.l <= 912.
            self.params[ 'nff' ] = self.nff.params
            self.components['nebular'] = None

            # Build the Synchrotron support only if
            # it is not already included in the SSP library
            if 'parsec22.NT' not in self.csp.ssp_lib :
                if syn is None :
                    syn = {}
                syn[ 'RCCSN' ] = self.csp.core.RCCSN()
                self.snsyn = SNSYN( self.csp.l, **syn )
                self.params[ 'syn' ] = self.snsyn.params
                self.components['synchrotron'] = None

        # If sub-gridding is required, build subgrid
        self.lstep = lstep
        if self.lstep is not None :
            self.lgrid = self.get_wavelength_grid(self.lstep)
        else : 
            self.lgrid = numpy.arange(self.csp.l.size)

        # If IGM transmission is required, build it
        # and compute transmission.
        self.igm = None
        self.igm_trans = numpy.ones_like(
            self.wl( obs = True )
        )
        if do_IGM :
            self.igm = IGM()
            self.igm_trans = self.igm.transmission(
                self.wl( obs = True ),
                self.redshift
            )

    def dump ( self ) :

        return dict(
            lstep    = self.lstep,
            do_Xray  = self.xrb is not None,
            do_Radio = self.nff is not None,
            do_AGN   = self.agn is not None,
            do_IGM   = self.igm is not None,
            cosmo    = { 'redshift' : self.cosmo.DL.get_x(),
                         'luminosity_distance' :  self.cosmo.DL.get_y(),
                         'age' :  self.cosmo.age.get_y() },
            csp = { 'ssp_lib' : self.csp.ssp_lib },
            **self.params
        )

    @classmethod
    def load ( cls, dictionary ) :
        return cls( **dictionary )
        
    def wl ( self, obs = False ) :
        """Returns the wavelength grid with the mask applied.
        
        Parameters
        ----------
        obs : bool
            (Optional, default = ``False``) if set to ``True``
            returns the observer's frame wavelength grid, otherwise
            the rest-frame grid is returned.
        
        Returns
        -------
        : 1d-array
        wavelength grid in observer's frame (``obs = True``) or in 
        rest frame (``obs = False``)
        """
        if obs :
            return ( 1 + self.redshift ) * self.csp.l[ self.lgrid ]
        return self.csp.l[self.lgrid]
    
    def get_wavelength_grid ( self, lstep ) :
        """ Reduces the granularity of the wavelength grid [optimization]
        
        Parameters
        ----------
        lstep : scalar int or boolean mask
            If the input is scalar, select one any `lstep` wavelengths.
            If the input is a boolean mask, only select the indexes of the
            array masked with lstep.
            
        Returns
        -------
        : integer array
            A list of indices of the wavelength grid
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
        """Preferred method for changing the value of free-parameters. 
        
        Parameters
        ----------
        age : float
            age of the galaxy
        redshift : float
            redshift of the galaxy
        sfh : dict
            (default = ``None``) arguments passed to SFH.set_parameters method, 
            if ``None`` is passed maintains the current parameterisation 
        csp : dict
            (default = ``None``) arguments passed to the CSP.set_parameters method, 
            if ``None`` is passed maintains the current parameterisation 
        ism : dict
            (default = ``None``) arguments passed to the ISM.set_parameters method, 
            if ``None`` is passed maintains the current parameterisation 
        agn : dict
            (default = ``None``) arguments passed to the AGN.set_parameters method, 
            if ``None`` is passed maintains the current parameterisation  
            If the GXY object has been built without AGN component this argument is ignored.
        nff : dict
            (default = ``None``) arguments passed to the NFF.set_parameters method, 
            if ``None`` is passed maintains the current parameterisation  
            If the GXY object has been built without AGN component this argument is ignored.
        syn : dict
            (default = ``None``) arguments passed to the SNSYN.set_parameters method, 
            if ``None`` is passed maintains the current parameterisation  
            If the GXY object has been built without AGN component this argument is ignored.

        Note
        ----
        The method is built to optimise the number of computations performed. 
        Do not pass arguments that would not change the current parameterisation.
        **e.g.1** passing ``sfh = {}`` is less performant than sticking to the 
        default ``sfh = None``.
        **e.g.2** passing the same value at each call is a waste of computational time: 
        ``redshift = 2.0`` passed at each call will considerably slow down execution.
        **Bottom line**: pass an argument only when necessary. 
        """

        # if age of sfh change, reset the csp's timetuple
        reset_csp = False
        reset_ism = False
        
        if redshift is not None :
            self.redshift = redshift
            self._z = 1. / ( 1 + self.redshift )
            self.UA = self.cosmo.age( self.redshift )
            self.params['redshift'] = self.redshift
            if self.igm is not None :
                self.igm_trans = self.igm.transmission(
                    self.wl( obs = True ),
                    self.redshift
                )
            
        if age is not None :
            # if the provided age is larger than the age of the Universe,
            # set the age to the age of the Universe
            self.age = age
            reset_csp = True
            reset_ism = True
            self.params['age'] = self.age

        if self.age > self.UA :
            raise RuntimeError( "Trying to set an age larger than the age "
                                "of the Universe at current redshift." )

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
                                     Mstar = self.sfh.Mstar( self.age ),
                                     Zstar = self.sfh.Zstar( self.age ) )

        ############################################################################
        # AGN varies only if some new parameters are passed
        if agn is not None :
            try :
                self.agn.set_parameters(**agn)
                self.params['agn'].update(self.agn.params)
            except AttributeError :
                raise AttributeError( 'Passing AGN-parameters to a GXY-class built '
                                      'without an AGN component is not allowed.' )
            except :
                raise

        ############################################################################
        # NFF varies both whether some new parameters are passed or if some
        # of the other components above which they do depend changes
        if self.nff is not None :
            if nff is None :
                nff = {}
            nff.update( { 'Zgas' : self.sfh.core.Zgas(self.age) } )
            self.nff.set_parameters( **nff )
            self.params['nff'].update(self.nff.params)

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
        """Unattenuated stellar emission in solar luminosities. 
        Approximates the integral
        
        .. math::
           
           L_\\lambda^\\text{CSP}(\\tau') = 
           \\int_0^{\\tau'}\\text{d}\\tau 
           L_\\lambda^\\text{SSP}\\bigl[\\tau, Z_\\ast(\\tau'-\\tau)\\bigr]\\psi(\\tau'-\\tau)

        where :math:`\\tau'` is the age of the galaxy, 
        :math:`L_\\lambda^\\text{SSP}[\\tau, Z\\ast]`
        is the luminosity of the Simple Stellar Population at given time :math:`\\tau` 
        and at given stellar metallicity :math:`Z_\\ast`, :math:`\\psi(\\tau)` is the 
        Star Formation History at time :math:`\\tau`
        
        Returns
        -------
        : 1-d array
        stellar emission on the default wavelength grid
        """
        
        # emission from stars (using directly the core function for performance)  
        return self.csp.core.emission( self.lgrid )

    def _passive_emission ( self, ftau = None ) :
        import warnings
        warnings.filterwarnings("error")
        if ftau is None :
            ftau = numpy.ones( self.csp.shape[:-1] )
        else :
            ftau = ftau.reshape( self.csp.shape[:-1] )
        
        itmin = find_nearest(self.csp.t, self.age-self.sfh.params['tau_quench']).item()
        itmax = max(find_nearest(self.csp.t, self.age).item(), itmin+1)
        itque = find_nearest(self.csp.t, self.sfh.params['tau_quench']).item()
    
        Mweight = trap_int(
            self.csp.t[:itque,numpy.newaxis], 
            (self.sfh( self.csp.t[:itque] ) * ftau[:,:itque]).T
        )
        if itmax-itmin > 1 :
            Lavg = self.csp.core.emission(self.lgrid, ftau) #.astype('float64')
            LocWeight = trap_int( 
                self.csp.t[itmin:itmax, numpy.newaxis], 
                (self.sfh( self.csp.t[itmin:itmax] ) * ftau[:,itmin:itmax]).T 
            )
            # Note that the '>1.e-20' is a NOT ELEGANT SOLUTION for the overflow bug
            # previously present in the following line
            # Lavg[wLW] /= LocWeight[wLW] -> RuntimeWarning: overflow encountered in divide
            # It is generated by Lavg being 'float32' and LocWeight 'float64' but 
            # if Lavg is casted to 'float32' another overflow arises in the loglikelihood computation
            # so this is not a perfect solution and should be investigated.
            wLW = (LocWeight > 1.e-20)
            Lavg[wLW] /= LocWeight[wLW]
            return Lavg * Mweight 
        return self.csp.L[:, itmax, self.csp._timetuple[2][0]] * Mweight
    
    def get_emission ( self, store_attenuation = False, **kwargs ) :
        """Computes the overall emission coming from a galaxy with given parameterisation.
        The resulting shape of the SED depends on the components the galaxy object has
        been built with. This authomatically deals with the interplay among the different 
        active components.

        Parameters
        ----------
        store_attenuation : bool
            (Optional, default = ``False``) if set to ``True``, stores the total,
            wavelength dependent, attenuation due to ISM in an internal variable (``Aavg``)

        **kwargs : dictionary, optional
            arguments passed to function ``set_parameters()``
        
        Returns
        -------
        : 1d-array
            the emission on the selected wavelength grid 
            in units of solar luminosities (:math:`[L_\odot]`)
        
        Note
        ----
        Even though only the overall emission is returned, the contribution of each
        component is stored in the internal dictionary attribute ``components``
        """
        if len( kwargs ) > 0 :
            self.set_parameters( **kwargs )
        
        # attenuation from ISM
        attTotMC, attTot = self.ism.total_attenuation( self.wl(), self.csp.t )
        
        # emission from stars (using directly the core function for performance)
        Lunatt, LattMC, Ltot = self.csp.core._kernel_emission( self.lgrid, attTotMC, attTot )
        # if self.params['age'] <= self.params['sfh']['tau_quench'] :
        #     Lunatt, LattMC, Ltot = self.csp.core._kernel_emission( self.lgrid, attTotMC, attTot )
        # else :
        #     Lunatt = self._passive_emission()
        #     LattMC = self._passive_emission( attTotMC )
        #     Ltot   = self._passive_emission( attTot )
        self.components['stellar'] = Lunatt
        self.components['extinct'] = numpy.array(Ltot) # copy the value

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
        """Returns the average attenuation in absolute magnitudes 
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
        """Returns the flux at given distance in units of milli-Jansky [mJy].

        .. math::
            
            S(\\lambda_O) = \\lambda_R^2 \cdot 
            \\dfrac{L_\\text{tot}(\\lambda_R) (1 + z)}{4\\;\\pi\\;c\\;D_L^2(z)} \cdot
            e^{-\\tau_\\text{IGM}(z)}

        where :math:`\\lambda_R` and :math:`\\lambda_O` are the rest-frame and observer's frame
        wavelength, respectively, :math:`L_\\text{tot}` is the total luminosity (as returned by
        function ``get_emission()``), :math:`z` is the redshift, :math:`c` the light-speed,
        :math:`D_L(z)` is the luminosity distance at observed redshift and, if it is included in 
        the model, :math:`e^{-\\tau_\\text{IGM}(z)}` is the IGM transmission at observed redshift.
        
        Returns
        -------
        : 1d-array
        SED flux in milli-Jansky
        """
        return self.cosmo.to_flux(
            self.redshift, self.wl(), self.get_emission()
        ) * self.igm_trans

    def components_to_flux ( self ) :
        """Utility function converting emissions in the internal ``components`` 
        dictionary to fluxes.
        
        Returns
        -------
        : dict
        copy of the internal dictionary ``components`` with the emission conveted 
        to fluxes in milli-Jansky.
        """
        return {
            k : self.cosmo.to_flux( self.redshift, self.wl(), c ) 
            for k, c in self.components.items()
        }
    

class PhotoGXY ( GXY ) :
    """Galaxy class including photometric system.
    This is a class derived from ``galapy.Galaxy`` which implements authomatic computation 
    of fluxes convolved with bandpass transmission filters.
    
    Parameters
    ----------
    pms : PMS instance, optional
      an instance of type ``galapy.PhotometricSystem.PMS``, if not passed,
      it should be built using function ``build_photometric_system``
    *args : tuple, optional
      arguments to be passed to the constructor of the base class ``GXY``
    **kwargs : dictionary, optional
      keyword arguments to be passed to the constructor of the base class ``GXY``
    
    See Also
    --------
    galapy.PhotometricSystem.PMS : class implementing the photometric system
    build_photometric_system : method for building a photometric system internally
    GXY : base class
    """

    def __init__ ( self, *args, pms = None, **kwargs ) :

        super().__init__( *args, **kwargs )
        self.pms = pms

    def dump ( self ) :
        return dict(
            pms_kwargs = self.pms.dump(),
            **super().dump()
        )

    @classmethod
    def load ( cls, dictionary ) :
        # ensure deep-copy
        dictionary = dict( dictionary )
        pms = PMS.load( dictionary.pop('pms_kwargs') )
        return cls( pms = pms, **dictionary )
    
    def build_photometric_system ( self, *args, **kwargs ) :
        """Forwards ``args`` and ``kwargs`` to the constructor of the photometric system.
        
        Parameters
        ----------
        *args : sequence
            positional arguments of the PMS contructor
        **kwargs : dictionary
            keyword arguments of the PMS constructor
        
        See Also
        --------
        galapy.PhotometricSystem.PMS : constructor of the photometric system class
        """

        self.pms = PMS( *args, **kwargs)
        return;

    def photoSED ( self ) :
        """Computes and returns the photometric bandpass fluxes in milliJansky.
        """

        if self.pms is None :
            raise Exception( "Photometric system has not been set. "
                             "Call function build_photometric_system() before." )
        return numpy.asarray(
            self.pms.get_fluxes( self.wl( obs = True ), self.get_SED() )
        )
        

class SpectralGXY ( GXY ) :
    """Not implemented yet
    """
    pass

