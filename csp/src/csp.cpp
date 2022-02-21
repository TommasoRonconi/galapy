#include <csp.h>

// =============================================================================

void sed::csp::set_params ( const sed::csp::vect_double & psi,
			    const sed::csp::vect_double & Zstar,
			    const std::vector< std::size_t > & iz_low,
			    const std::size_t it_last ) noexcept {

  _psi = psi; _Zstar = Zstar;
  _iz_low = iz_low; _it_last = it_last;
  
  return;

}

// =============================================================================

double sed::csp::emission ( const std::size_t il,
			    const double * const Tfact ) const noexcept {

  double Lout = 0.;
  for ( std::size_t it = 1; it <= _it_last; ++it )
      Lout += Tfact[ it ] * _sfh_lum_timeintegrand( il, it );
  return Lout;

}

// =============================================================================
