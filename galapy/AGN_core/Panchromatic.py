#######################################################################################
# External imports

import numpy
from scipy.ndimage import gaussian_filter
from scipy.special import erf, expit, gamma as Gamma, zeta
from scipy.constants import h,c,k
# [hc/k] = Ang * K
hc_k = 1.e+10 * h * c / k  

#######################################################################################
# Internal imports

from galapy.internal.utils import find_nearest, trap_int
from galapy.internal.interp import lin_interp
# from galapy.internal.abc import Model
from galapy.internal.constants import clight, Lsun, Ang_to_keV

#######################################################################################
# Support functions

def sigmoid ( x, k, tol=1.e-7 ) :
    """Sigmoid filter, returns zero if value lower than tolerance"""
    ret = expit(k*x)
    return ret * (ret>tol)
def inv_sigmoid ( x, k, tol=1.e-7 ) :
    """Inverse Sigmoid filter, returns zero if value higher than tolerance"""
    ret = 1.0 - expit(k*x)
    return ret * (ret>tol)
def gaussian ( x, k ) :
    """Gaussian filter, returns zero if value lower than tolerance"""
    return 0.5*(1.+erf(x/k))
def safe_exp(x):
    """
    Compute a numerically stable exponential using the logistic sigmoid (expit).

    Parameters
    ----------
    x : array_like
        Input value or array.

    Returns
    -------
    out : ndarray or scalar
        The exponential of `x`, computed safely without overflow.
    """
    x = numpy.asarray(x)
    
    # exp(x) = sigmoid(x) / sigmoid(-x)
    # where sigmoid(x) = expit(x) is numerically stable
    # To avoid division by zero, we add a very small epsilon
    eps = numpy.finfo(x.dtype).tiny if numpy.issubdtype(x.dtype, numpy.floating) else 1e-300
    # (alternatively we can fix the small epsilon to the following minimum
    # but the above expression costs only 5 micro-sec)
    # eps = 1.e-300 
    sigma_pos = expit(x)
    sigma_neg = 1.0 - sigma_pos  # equivalent to expit(-x) but cheaper and stable
    
    return sigma_pos / numpy.maximum(sigma_neg, eps)

#######################################################################################
# A generic class implementing a normalized piecewise power-law 

class piecewise_powerlaw () :
    
    def __init__ ( self, lgrid, 
                   expnt = [1.2, 0.8, -0.5, -3.0], 
                   bound = [12.5, 500., 1250., 5.e+4, numpy.inf],
                   sigma = .3, tol = 1.e-7
                 ) :
        self.lgrid  = lgrid
        self.smooth = 1./sigma
        self.llgrid = numpy.log(lgrid, dtype=float)
        self.lbound = numpy.log(bound, dtype=float)
        self.expnt  = numpy.array(expnt, dtype=float)
        
        # compute weights
        self.weight = numpy.array([
            sigmoid( self.llgrid - self.lbound[i], self.smooth, tol ) *
            inv_sigmoid( self.llgrid - self.lbound[i+1], self.smooth, tol )
            for i in range(self.lbound.size-1) 
        ])
        
        # compute coefficients
        self.coeff = numpy.zeros_like(self.expnt)
        indices = find_nearest( self.llgrid, self.lbound[:-1] )
        for ii in range(self.expnt.size-2, -1, -1) :
            jj = indices[ii+1]
            self.coeff[ii] = (
                self.coeff[ii+1] +
                (self.expnt[ii+1]-self.expnt[ii]) *
                self.llgrid[jj]
            )
            
        # compute emission
        self._emission = (
            numpy.exp(
                self.coeff[:,numpy.newaxis] +
                self.expnt[:,numpy.newaxis] * 
                self.llgrid 
            ) * self.weight
        ).sum(axis=0)
        
        # normalize emission
        # self._emission /= trap_int(self.lgrid, self._emission/self.lgrid)
        # the integral is in ln(lambda) because we consider nu*L(nu) -> L(lambda)/lambda
        self._emission /= trap_int(self.llgrid, self._emission) 
        
#######################################################################################
# Support functions

def get_view_fact ( theta ) :
    costheta = numpy.cos(theta)
    return 6./7. * costheta * ( 1 + 2 * costheta )
        
#######################################################################################
# Disk emission

class agn_disk ( piecewise_powerlaw ) :
    
    def __init__ (
            self, *args,
            Lbol = 1.e+44, theta_view = 0.0, delta = -0.36,
            **kwargs ) :
        kwargs.update(
            expnt = [1.2, 0.8, -0.5+delta, -3.0], 
            bound = [12.5, 500., 1250., 5.e+4, numpy.inf]
        )
        super().__init__(*args, **kwargs) 
        log_pw = -99 * numpy.ones_like(self._emission)
        # safe log:
        numpy.log(
            self._emission, out = log_pw,
            where=(self._emission>0.0)
        )
        # interpolator:
        self._fcall = lin_interp( 
            self.llgrid, 
            log_pw
        )
        # set free parameters
        _ = self.set_params( Lbol = Lbol, theta_view = theta_view ) 
        
    def set_params ( self, Lbol = None, theta_view = None, delta = None ) :
        changed = False
        if Lbol is not None :
            self.Lbol = Lbol
            changed |= True
        if theta_view is not None :
            self.theta_view = theta_view
            costheta = numpy.cos(self.theta_view)
            self.view_fact = get_view_fact( self.theta_view )
            changed |= True
        if delta is not None :
            pass # PLACEHOLDER <----------------------------------------------- !!!
        return changed
        
    def __call__ ( self, ll, **kwargs ) :
        _ = self.set_params(**kwargs)
        ret1 = self._fcall(numpy.log(ll))
        ret2 = numpy.zeros_like(ret1)
        # return zero values when the function is not defined:
        return self.Lbol * self.view_fact * numpy.where(
            ret1>-99, numpy.exp(ret1), ret2
        )

#######################################################################################
# AGN attenuation

_f37 = 0.42857142857142855
_f47 = 0.57142857142857140

class agn_attenuation () :

    ll_pr84 = numpy.log10([
        1275.,1330.,1385.,1435.,1490.,1545.,1595.,
        1647.,1700.,1755.,1810.,1860.,1910.,2000.,
        2115.,2220.,2335.,2445.,2550.,2665.,2778.,
        2890.,2995.,3105.,3704.,4255.,5291., 5292.
    ])
    AV_pr84 = numpy.array([
        13.54,12.52,11.51,10.80, 9.84, 9.28, 9.06, 
         8.49, 8.01, 7.71, 7.17, 6.90, 6.76, 6.38,
         5.85, 5.30, 4.53, 4.24, 3.91, 3.49, 3.15, 
         3.00, 2.65, 2.29, 1.67, 1.00, 0.00, 0.00
    ])
    
    def __init__ ( self, Delta = 0.7, EBV = 0.02 ) :
        self._fSMC_pr84 = lin_interp(
            type(self).ll_pr84, type(self).AV_pr84
        )
        _ = self.set_params( Delta = Delta, EBV = EBV )
        self.Av = self._Av( 10.**type(self).ll_pr84 )

    def set_params ( self, Delta = None, EBV = None ) :
        changed = False
        if Delta is not None :
            self.Delta = Delta
            sinD = numpy.sin( self.Delta )
            self.fact_trigonometric1 = _f47 + _f37 * numpy.cos(self.Delta)**2 - _f47 * sinD*sinD*sinD 
            self.fact_trigonometric2 = _f37 * sinD*sinD + _f47 * sinD*sinD*sinD
            changed |= True
        if EBV is not None :
            self.EBV = EBV
            changed |= True
        return changed

    def __str__ ( self ) :
        return f'{type(self).__name__}(Delta = {self.Delta:.4f}, EBV = {self.EBV:.4f})'

    def __repr__ ( self ) :
        return f'{type(self).__name__}(Delta = {self.Delta}, EBV = {self.EBV})'

    def _Av ( self, ll ) :
        return 10**(
            -0.4 * self.EBV * self._fSMC_pr84( numpy.log(ll)*0.43429448190325176 )
        )
            
    def attenuation ( self, ll, theta_view = 0.0, **kwargs ) :
        _ = self.set_params( **kwargs )
        self.Av = self._Av( ll )
        if theta_view < 0.5*numpy.pi - self.Delta :
            return self.Av
        return numpy.zeros_like( ll )
        

#######################################################################################
# AGN dust emission: torus + polar

class agn_dust () :
    
    def __init__ ( self, lgrid, TH = 1500., 
                   expnt = [0.8, -0.8, -1.5], 
                   bound = [5.e+4, 2.e+5, 4.e+5, numpy.inf], 
                   **kwargs ) :
        
        # set lambda grid
        self.lgrid = lgrid # shallow copy
        
        # piecewise power-law
        _pw = piecewise_powerlaw(
            self.lgrid, expnt = expnt[:-1], bound = bound[:-1], **kwargs
        )
        # normalized at 1/2 of total contribution
        _logpw = -99.0 * numpy.ones_like( _pw._emission )
        # safe log
        numpy.log(
            _pw._emission, out = _logpw, where = ( _pw._emission > 0.0 )
        )
        self._pw_fcall = lin_interp( _pw.llgrid, _logpw )
        self._idx = find_nearest( self.lgrid, bound[-2] )
        self._pw_fact = _pw._emission[self._idx]
        
        # modified black-body
        self.gamma = -( expnt[-1] - 3 )
        _bb_lnum = -self.gamma * _pw.llgrid
        # normalized at 1/2 of total contribution
        _bb_lnorm = numpy.log(
            0.5 / ( (1./hc_k)**self.gamma * Gamma(self.gamma+1) * zeta(self.gamma+1) )
        )
        self._bb_fact_num = numpy.exp( _bb_lnorm + _bb_lnum[self._idx] )
        self._bb_fcall = lin_interp(
            _pw.llgrid,
            _bb_lnorm + _bb_lnum 
        )
        self._bb_weight = lin_interp( self.lgrid, sigmoid( 
            _pw.llgrid - _pw.lbound[-1], _pw.smooth, 
            tol = kwargs.get('tol', 1.e-7) 
        ))
        _ = self.set_params(TH)
        
    def set_params ( self, TH = None ) :
        changed = False
        if TH is not None : 
            self.TH = TH
            self._continuity = self._pw_fact / (
                self._bb_fact_num *
                self._bb_den( self.lgrid[self._idx] )
            )
            changed |= True
        return changed
            
    def _bb_den ( self, ll ) :
        return ( self.TH**self.gamma ) / ( safe_exp(
            hc_k / ( self.TH * ll )
        ) - 1.0 )
            
    def __call__ ( self, ll, **kwargs ) :
        _ = self.set_params(**kwargs)
        logll = numpy.log(ll)
        return numpy.exp(self._pw_fcall(logll)) + self._continuity * numpy.exp(
            self._bb_fcall(logll)
        ) * self._bb_den(ll) * self._bb_weight(ll)

#######################################################################################
# Radio Jets

class agn_jets () :
    
    # shape 
    aflat = 0.1
    asteep = 0.8

    # fixed parameters
    glorentz = 10.
    Rext = 10.
    
    def __init__ ( 
        self, lgrid, 
        RL = 1.e+4, theta_view = 0.0, 
        LAD2500 = 1.e+44, 
        sigma = 0.3, tol = 1.e-18
    ) :
        self.lgrid = lgrid
        self.llgrid = numpy.log(lgrid)
        self.smooth = lambda x : gaussian_filter(
            x, 0.3 / ((self.llgrid.max()-self.llgrid.min())/self.llgrid.size) * numpy.log(10.)
        )
        self._ngrid5GHz = clight['A/s']/(self.lgrid*5.e+9)
        self._view_fact30 = get_view_fact(numpy.radians(30.0))
        self._glorentz2 = numpy.sqrt(type(self).glorentz*type(self).glorentz-1)
        self._weight = ((1.e+6<self.lgrid)&(self.lgrid<1.e+10)).astype(float)
        self._fflat = (
            2.0*inv_sigmoid(numpy.exp(
                (-type(self).aflat-2.5)*( numpy.log(clight['A/s']) - self.llgrid - numpy.log(2.e+8))
            ), 1.0, tol )
        )
        self._fsteep = (
            2.0*inv_sigmoid(numpy.exp(
                (-type(self).asteep-2.5)*( numpy.log(clight['A/s']) - self.llgrid - numpy.log(2.e+8))
            ), 1.0, tol )
        )
        _ = self.set_params( RL = RL, theta_view = theta_view, LAD2500 = LAD2500 )
        # 
        self._nfact = 1./(1.+numpy.sqrt(0.25*self._ngrid5GHz))
        self._fcall_flat = lin_interp(
            self.lgrid, self.smooth(
                self._weight * self._fflat * self._ngrid5GHz**(1-type(self).aflat) * self._nfact
            )
        )
        self._fcall_steep = lin_interp(
            self.lgrid, self.smooth(
                self._weight * self._fsteep * type(self).Rext * self._ngrid5GHz**(1-type(self).asteep) * self._nfact
            )
        )
            
    def set_params ( self, RL = None, theta_view = None, LAD2500 = None ) :
        changed = False
        if RL is not None :
            self.RL = RL
            changed |= True
        if theta_view is not None :
            self.theta_view = theta_view
            self.flat_fact = 1./(
                type(self).glorentz-self._glorentz2*numpy.cos(self.theta_view)
            )**(2+type(self).aflat)
            self._L5GHz = 1./(self.flat_fact+type(self).Rext)/(1+numpy.sqrt(0.25))
            changed |= True
        if LAD2500 is not None :
            self.LAD2500 = LAD2500*self._view_fact30
            changed |= True
        return changed
    
    def __call__ ( self, ll, **kwargs ) :
        _ = self.set_params(**kwargs)
        return self.RL * self.LAD2500 * self._L5GHz * (
            self.flat_fact * self._fcall_flat(ll) +
            self._fcall_steep(ll)
        ) * (2500*(5.e+9/clight['A/s']))

#######################################################################################
# AGN X-ray corona

class agn_xray () :
    
    Ecut = 3.e+2
    gamma = 1.8
    
    def __init__ ( 
        self, lgrid, 
        Lbol = 1.e+44, theta_view = 0.0,
        norm = 'bolometric', isotropic = True, 
        sigma = 0.3 
    ) :
        
        self.lgrid = lgrid
        self.llgrid = numpy.log(self.lgrid)
        self.El = Ang_to_keV(self.lgrid)
        self.norm = norm
        self.iso = isotropic
        self.smooth = lambda x : gaussian_filter(
            x, 0.3 / ((self.llgrid.max()-self.llgrid.min())/self.llgrid.size) * numpy.log(10.)
        )

        if self.El.max() < 2. :
            raise RuntimeError( "Cannot build the X-ray spectrum for "
                                "a wavelength grid starting at lambda > "
                                "6 Angstrom ~ 2 keV! "
                                "Set a smaller `lmin` value." )
        
        # compute lambda * L(lambda)
        lL = self.El**(-type(self).gamma+2) * safe_exp(-self.El / type(self).Ecut )
        
        # find interval for hard-X normalization
        wE = ( 2. <= self.El ) & ( self.El <= 10. )

        # compute emission law normalization
        lLnorm = 1. / trap_int( self.llgrid[wE], lL[wE] )
        
        # Compute the normalized power-law:
        # - high energy cut-off at 300 keV
        # - spectral index fixed to 1.8
        # - low energy cut-off at 50 Angstrom (~0.25 keV)
        ret = numpy.zeros_like( self.lgrid )
        wL = ( self.lgrid <= 5.e+1 ) 
        ret[wL] = lL[wL] * lLnorm

        # store normalized interpolator object for X-ray emission
        self._fcall = lin_interp( self.lgrid, self.smooth( ret ) )
        
        _ = self.set_params( Lbol = Lbol, theta_view = theta_view )
        
    def get_normalization ( self, norm = None ) :
        
        if self.norm == 'bolometric' :
            return 10.96 * ( 1. + ( numpy.log10( self.Lbol/Lsun ) / 11.48 )**17.79 )
           
    def set_params ( self, Lbol = None, theta_view = None ) :
        changed = False
        if Lbol is not None :
            self.Lbol = Lbol
            self.Lnorm = self.Lbol / self.get_normalization()
            changed |= True
        if theta_view is not None :
            pass
        return changed
    
    def __call__ ( self, ll, **kwargs ) :
        return self.Lnorm * self._fcall( ll )

#######################################################################################
# AGN components wrapping model

def agn_build_params ( **kwargs ) :
    """ Standard function for building the parameters dictionary of class Panchromatic()
    """

    return {
        'Lbol'       : kwargs.get('Lbol',       1.e+44),
        'theta_view' : kwargs.get('theta_view', 0.0   ),
        'delta'      : kwargs.get('delta',     -0.36  ),
        'TH'         : kwargs.get('TH',      1500.    ),
        'Delta'      : kwargs.get('Delta',      0.7   ),
        'EBV'        : kwargs.get('EBV',         .02  ),
        'RL'         : kwargs.get('RL',         1.e+4 ),
    }

class Panchromatic () :

    def __init__ (
            self, lgrid,
            Lbol = 1.e+44, theta_view = 0.0, delta = -0.36,
            TH = 1500.,
            Delta = 0.7, EBV = .02,
            RL = 1.e+4,
            do_Xray = True
    ) :
        self.lgrid = numpy.array( lgrid )
        self.disk = agn_disk( self.lgrid )
        self.attn = agn_attenuation()
        self.dust = agn_dust( self.lgrid, sigma=.2 )
        self.jets = agn_jets( self.lgrid, sigma=.3 )
        self.xray = None
        self.components = {
            'disk_int' : numpy.zeros_like( self.lgrid ),
            'disk_obs' : numpy.zeros_like( self.lgrid ),
            'dust_obs' : numpy.zeros_like( self.lgrid ),
            'jets_obs' : numpy.zeros_like( self.lgrid ),
        }
        if do_Xray :
            self.xray = agn_xray( self.lgrid )
            self.components['xray_obs'] = numpy.zeros_like( self.lgrid ),
        self.params = agn_build_params(
            Lbol = Lbol, theta_view = theta_view, delta = delta,
            TH = TH,
            Delta = Delta, EBV = EBV,
            RL = RL 
        )
        self.set_parameters( **self.params )

    def set_parameters (
            self,
            Lbol = None, theta_view = None, delta = None,
            TH = None,
            Delta = None, EBV = None,
            RL = None,
            **kwargs
    ) :
        LAD2500 = None
        if self.disk.set_params( Lbol = Lbol, theta_view = theta_view, delta = delta ) :
            LAD2500 = self.disk(2500.) / self.disk.view_fact
        self.dust.set_params( TH = TH )
        self.attn.set_params( Delta = Delta, EBV = EBV )
        self.jets.set_params( RL = RL, theta_view = theta_view, LAD2500 = LAD2500 )
        if self.xray is not None :
            self.xray.set_params( Lbol = Lbol )
        self.params.update( **{
            k : v
            for k, v in locals().items()
            if k in set(self.params.keys())
            and v is not None
        } )
        return;

    def __call__ ( self, ll, **kwargs ) :
        self.set_parameters( **kwargs )
        self.components['disk_int'] = self.disk( ll )
        self.components['disk_obs'] = self.attn.attenuation(
            ll, theta_view = self.disk.theta_view ) * self.components['disk_int']
        self.components['dust_obs'] = self.dust( ll ) * trap_int(
            ll, self.components['disk_int'] / ( self.disk.view_fact * ll ) * (
                ( 1 - self.attn.Av) * self.attn.fact_trigonometric1 +
                self.attn.fact_trigonometric2 )
        )
        self.components['jets_obs'] = self.jets( ll )
        if self.xray is not None :
            self.components['xray_obs'] = self.xray( ll )
            return (
                self.components['disk_obs'] +
                self.components['dust_obs'] +
                self.components['jets_obs'] +
                self.components['xray_obs']
            ) / ll
        return (
            self.components['disk_obs'] +
            self.components['dust_obs'] +
            self.components['jets_obs']
        ) / ll

    def emission ( self, *args, **kwargs ) :
        return self.__call__( *args, **kwargs )
        
#######################################################################################

