#include <ism.h>

// =============================================================================
// =============================================================================
// Non-member functions
// =============================================================================

double sed::total_attenuation ( const double wavelength,
				const double att_dd_ll,
				const double att_mc_ll,
				const double eta_mc_tt ) {

  return -2.5 * sed::cnst::nl_10 *
    std::log( att_dd_ll * ( 1 - eta_mc_tt * ( 1 - att_mc_ll ) ) );

}

// =============================================================================
// =============================================================================
// Base ISM class functions
// =============================================================================

double sed::ism::temperature ( const double Etot ) noexcept {

  if ( Etot > 0. ) {
    std::size_t NL = 200, NT = 200;
    auto Lvec = utl::log_vector( NL,
				 std::pow( 10, _logLmin ),
				 std::pow( 10, _logLmax ) );
    auto Tvec = utl::log_vector( NT,
				 std::pow( 10, _logTmin ),
				 std::pow( 10, _logTmax ) );
    std::vector< double > Et ( NT, 0. );

    // REGULARIZATION OF THE INTEGRAL:
    // dividing everything by Etot guarantees number approaching 0
    // without scattering too many orders of magnitude around the origin
    double totE = 1. / Etot;

    // Subtrahend is defined as
    // totE * Etot / NL = Etot / ( NL * Etot ) = 1 / NL
    double subFact = 1. / NL;
  
    double kL, dL;
    for ( std::size_t iL = 1; iL < NL; ++iL ) {
      kL = _paramsrc[ 1 ] * ( 1 - attenuation( Lvec[ iL ] ) );
      dL = ( Lvec[ iL ] - Lvec[ iL - 1 ] ); // dLambda in [ cm ] units
      for ( std::size_t iT = 0; iT < NT; ++iT )
	Et[ iT ] += totE * dL * kL *
	  utl::blackbody( 1.e-8 * Lvec[ iL ], Tvec[ iT ] )
	  - subFact;
    }
    // LINEAR INTERPOLATION:
    // [ maybe here I could use an interpolator instead ]
    auto idx = utl::find_low( 0., Et );
    _Temp = utl::line_from_2points( 0.,
				    Et[ idx ], Tvec[ idx ],
				    Et[ idx + 1 ], Tvec[ idx + 1 ] );

    // setting lambda wien's displacement
    _lambda_wien = 3.e+7 / _Temp;
    
    return _Temp;
  }
  else return 0.;

}

// =============================================================================
// =============================================================================
// DIFFUSE Phase
// =============================================================================

double sed::diffuse::_A_Vband ( const double * const param ) const noexcept {

  return 2.e-2 *                  // 2 * ( 1 / 10^8 ) * ( 10^{-3} )^-2
    ( 1 - param[ 0 ] ) *          // ( 1 - f_MC )
    param[ 1 ] *                  // N_V^diff
    param[ 2 ] /                  // M_dust
    ( param[ 3 ] * param[ 3 ] );  // ( R_dust )^{-2}

}

// =============================================================================

double sed::diffuse::_fact_greybody ( const double * const param ) const noexcept {

  return
    16 * sed::cnst::pi * sed::cnst::pi / 3 * // 16 pi^2 / 3
    1.e-8 *                                  // [ cm^-1 ] -> [ Ang^-1 ]
    ( sed::cnst::pc * sed::cnst::pc *        // [ pc ] -> [ cm ]
      param[ 3 ] * param[ 3 ] );             // ( R_dust )^{2}

}

// =============================================================================

double sed::diffuse::emission ( const double lambda ) const noexcept {

  return sed::cnst::solL *
    ( this->_emission( lambda, _Temp, _paramsrc[ 1 ] ) + // dust contribute
      _paramsrc[ 2 ] * _Labs * _fpah( lambda ) );        // PAH contribute

}

// =============================================================================
// =============================================================================
// CLOUD Phase
// =============================================================================

double sed::cloud::_A_Vband ( const double * const param ) const noexcept {

  return 2.56e-4 *                         // ( 16 * 16 / 10^6 )
    param[ 1 ] *	                   // N_V^MC
    param[ 4 ] * sed::cnst::solZ *         // Z_gas / Z_sol
    param[ 0 ] * param[ 6 ] /              // M_MC = f_MC * M_gas
    ( param[ 2 ] + 1.e-7 ) /               //        / ( N_MC + 1.e-7 ) 
    ( param[ 3 ] * param[ 3 ] );           // ( R_MC )^{-2}

}

// =============================================================================

double sed::cloud::_fact_greybody ( const double * const param ) const noexcept {

  return
    16 * sed::cnst::pi * sed::cnst::pi / 3 * // 16 pi^2 / 3
    1.e-8 *                                  // [ cm^-1 ] -> [ Ang^-1 ]
    sed::cnst::pc * sed::cnst::pc *          // [ pc ] -> [ cm ]
    param[ 3 ] * param[ 3 ] *                // R_MC^2
    param[ 2 ];                              // N_MC
    

}

// =============================================================================

double sed::cloud::emission ( const double lambda ) const noexcept {

  return sed::cnst::solL * this->_emission( lambda, _Temp, _paramsrc[ 1 ] );

}

// =============================================================================

double sed::cloud::time_attenuation ( const double lambda,
				      const double tau ) const noexcept {
  
  return 1 - this->eta( tau ) * ( 1. - this->attenuation( lambda ) );
  
}

// =============================================================================

double sed::cloud::average_attenuation ( const double lambda,
					 const double tau ) const noexcept {

  // place-holder operation to shut-up 'warning: unused parameter':
  return lambda * tau;

}

// =============================================================================
