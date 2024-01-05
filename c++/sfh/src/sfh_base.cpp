#include <sfh_base.h>

// =============================================================================
// Definitions of sed::sfh_base class functions
// =============================================================================

double sed::sfh_base::operator() ( const double tt ) const noexcept {

  return _sfr( tt ) *							\
    ( 1 - utl::heaviside< double, double >( tt, _tau_quench ) );

}

// =============================================================================

void sed::sfh_base::model ( const double * const tau,
			    double * const psi,
			    const std::size_t size ) const noexcept {
  
  for ( size_t ii = 0; ii < size; ++ii )
    psi[ ii ] = ( *this )( tau[ ii ] );
  
  return;

}

// =============================================================================

double sed::sfh_base::get_Mstar ( const double tau,
				  const std::size_t npoints ) const noexcept {

  std::vector< double > tvec = utl::log_vector( npoints, 1., tau );
  double Mstar = 0;
  for ( std::size_t it = 1; it < npoints; ++it ) {
    double tloc = 0.5 * ( tvec[ it ] + tvec[ it - 1 ] );
    Mstar +=
      ( tvec[ it ] - tvec[ it - 1 ] ) *
      ( 1 - sed::R_chabrier( tau - tloc ) ) *
      ( *this )( tloc );
  }
  return Mstar;
  
}

// =============================================================================

void sed::sfh_base::eval ( const double * const tau,
			   double * const psi,
			   const std::size_t size,
			   const double * const param ) noexcept {
  
  if ( param ) set_params( param );
  model( tau, psi, size );
  return;
      
}

// =============================================================================

// This is the version for PyBind11
void sed::sfh_base::time_grid ( const double age,
				const std::vector< double > & tgrid,
				const std::vector< double > & Zgrid,
				std::vector< double > & out_dMgrid,
				std::vector< double > & out_Zgrid,
				std::vector< std::size_t > & out_Zidx,
				std::size_t & out_last_idx ) {
  
  // 1) compute star formation history
  // 2) compute star-metal enrichment history
  // 3) find metallicity interval for interpolation  
  std::fill( out_dMgrid.begin(), out_dMgrid.end(), 0. );
  std::fill( out_Zgrid.begin(),  out_Zgrid.end(),  0. );
  std::fill( out_Zidx.begin(),   out_Zidx.end(),   0  ); 
  out_last_idx = 0;

  double time = age - tgrid[ out_last_idx ], next_time;
  double dM = 0.0;
  std::fill( _dM.begin(), _dM.end(), 0. ); // reset _dM vector
  std::size_t jj = 0;
  while ( time > 0. ) {
    next_time = age - tgrid[ out_last_idx + 1 ];
    while ( time > next_time ) {
      // the control factor '*( time > 0.0 )' guarantees
      // no negative values of time are summed up in the process
      _dM[++jj] = ( *this )( time ) * 1.e+5 * ( time > 0.0 ); // dt=10^5 is fixed
      dM += _dM[jj];
      time -= 1.e+5; // dt=10^5 is fixed
    }
    /*(1)*/ out_dMgrid[ out_last_idx ] = dM;
    /*(2)*/ out_Zgrid[ out_last_idx ]  = ( *this ).get_Zstar( time );
    /*(3)*/ out_Zidx[ out_last_idx ]   = utl::find_low( out_Zgrid[ out_last_idx ], Zgrid );
    
    ++out_last_idx; // increment index
    time = next_time; // update time
  }
  
  return;

}

// =============================================================================

// This is the version for CPython
// void sed::sfh_base::time_grid ( const double age,
// 				const double * const tgrid,
// 				const std::size_t tgrid_size,
// 				const double * const Zgrid,
// 				const std::size_t Zgrid_size,
// 			        double * const * const out_psigrid,
// 			        double * const * const out_Zgrid,
// 			        std::size_t * const * const out_Zidx,
// 				std::size_t * const out_last_idx ) {
  
//   // 1) compute star formation history
//   // 2) compute star-metal enrichment history
//   // 3) find metallicity interval for interpolation  
//   std::fill( *out_psigrid, *out_psigrid + tgrid_size, 0. );
//   std::fill( *out_Zgrid,   *out_Zgrid   + tgrid_size, 0. );
//   std::fill( *out_Zidx,    *out_Zidx    + tgrid_size, 0  ); 
//   *out_last_idx = 0;

//   double time = age - tgrid[ *out_last_idx ];
//   while ( time > 0. ) {
//     /*(1)*/ (*out_psigrid)[ *out_last_idx ]  = ( *this )( time );
//     /*(2)*/ (*out_Zgrid)[ *out_last_idx ]    = ( *this ).get_Zstar( time );
//     /*(3)*/ (*out_Zidx)[ *out_last_idx ] = utl::find_low( (*out_Zgrid)[ *out_last_idx ],
// 							  Zgrid, Zgrid_size );
//     ++( *out_last_idx ); // increment index
//     time = age - tgrid[ *out_last_idx ]; // update time
//   }
  
//   return;

// }

// =============================================================================
