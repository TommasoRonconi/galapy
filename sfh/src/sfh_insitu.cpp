#include <sfh_insitu.h>

// =============================================================================

double sed::sfh_insitu::_sfr ( const double tau ) const noexcept {

  double xx = tau * _paramsrc[ 3 ];
  return _paramsrc[ 0 ] *						\
    ( std::exp( - sed::sfh_insitu::_i_s * xx ) -			\
      std::exp( - xx * ( _paramsrc[ 1 ] - sed::R_chabrier_inst ) ) );
  
}

// =============================================================================

double sed::sfh_insitu::_dust_to_gas ( const double tau ) const noexcept {

  double xx = tau * _i_s * _paramsrc[ 3 ];
  double fact = _s * ( _paramsrc[ 1 ] - sed::R_chabrier_inst ) - 1;
  double eps = _kSN + ( _eps_acc * _s * _yD ) / ( fact + _s * _kSN );
  double factx = fact * xx;
  double sex = _s * eps * xx;
  return _s * _s * _s * _eps_acc * _yD * _yZ /
    ( fact * ( fact + _s * _kSN ) * ( fact + _s * eps ) ) *
    ( 1 - factx / ( std::exp( factx ) - 1 ) *
      ( 1 + fact / ( _s * eps ) * ( 1 - ( 1 - std::exp( - sex ) ) / sex ) ) );

}

// =============================================================================

double sed::sfh_insitu::get_Mgas ( const double tau ) const noexcept {

  return ( *this )( tau ) * _paramsrc[ 2 ];

}

// =============================================================================

double sed::sfh_insitu::get_Mdust ( const double tau ) const noexcept {

  return _dust_to_gas( tau ) * get_Mgas( tau );

}

// =============================================================================

double sed::sfh_insitu::get_Zgas ( const double tau ) const noexcept {

  if ( tau >= 1.e+5 ) {
    double xx = tau * _i_s * _paramsrc[ 3 ];
    double fact = _s * ( _paramsrc[ 1 ] - sed::R_chabrier_inst ) - 1;
    return _s * _yZ / fact * ( 1 - fact * xx / ( std::exp( fact * xx ) - 1 ) );
  }
  else
    return 0.;

}

// =============================================================================

double sed::sfh_insitu::get_Zstar ( const double tau ) const noexcept {

  if ( tau >= 1.e+5 ) {
    double xx = tau * _i_s * _paramsrc[ 3 ];
    double gamm =_paramsrc[ 1 ] - sed::R_chabrier_inst;
    double fact1 = _s * gamm;
    double fact2 = fact1 - 1;
    double exp1 = std::exp( - xx ), exp2 = std::exp( - fact1 * xx );
  
    return _yZ / gamm *
      ( 1 - fact1 / fact2 *
	( exp1 - exp2 * ( 1 + fact2 * xx ) ) /
	( fact2 + exp2 - fact1 * exp1 ) );
  }
  else
    return 0.;
}

// =============================================================================
