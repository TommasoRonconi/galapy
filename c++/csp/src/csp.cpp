#include <csp.h>

// =============================================================================

sed::csp::csp ( const sed::csp::vect_double lambda,
		const sed::csp::vect_double tau,
		const sed::csp::vect_double Z,
		const sed::csp::vect_double LltZ,
		const int do_CCSN_rate ) noexcept : 
  _lambda{ lambda }, _Nlambda{ lambda.size() },
  _tau{ tau }, _Ntau{ tau.size() },
  _Z{ Z }, _NZ{ Z.size() },
  _LltZ{ LltZ }, _NL{ LltZ.size() }
{

  // if boolean true compute table for RCCSN:
  if ( do_CCSN_rate ) {
    const utl::interpolator< utl::lin_interp > fR0 { sed::ZZ_CCSN, sed::R0_CCSN };
    const utl::interpolator< utl::lin_interp > fR1 { sed::ZZ_CCSN, sed::R1_CCSN };
    _RCCSNtZ = sed::csp::vect_double( _Ntau * _NZ );
    for ( std::size_t it = 0; it < _Ntau; ++it )
      for ( std::size_t iz = 0; iz < _NZ; ++iz )
	_RCCSNtZ[ it * _NZ + iz ] =
	  1.e-9 * fR0( _Z[ iz ] ) * std::pow( 1.e-6 * _tau[ it ], -fR1( _Z[ iz ] ) );
  }
  
}

// =============================================================================

void sed::csp::set_params ( const sed::csp::vect_double & dM,
			    const sed::csp::vect_double & Zstar,
			    const std::vector< std::size_t > & iz_low,
			    const std::size_t it_last ) noexcept {

  _dM = dM; _Zstar = Zstar;
  _iz_low = iz_low; _it_last = it_last;
  
  return;

}

// =============================================================================

double sed::csp::emission ( const std::size_t il,
			    const double * const Tfact ) const noexcept {

  double Lout = 0.;
  for ( std::size_t it = 1; it < _it_last; ++it )
      Lout += Tfact[ it ] * _sfh_lum_timeintegrand( il, it );
  return Lout;

}

// =============================================================================

std::vector< double > sed::csp::kernel_emission ( const std::size_t il,
						  const double * const MCfact,
						  const double * const DDfact )
  const noexcept {
  
  std::vector< double > Lout = { 0., 0., 0.};
  double dL;
  for ( std::size_t it = 1; it < _it_last; ++it ) {
    dL = _sfh_lum_timeintegrand( il, it );
    Lout[0] += dL;
    Lout[1] += MCfact[ it ] * dL;
    Lout[2] += DDfact[ it ] * dL;
  }
  return Lout;

}

// =============================================================================
