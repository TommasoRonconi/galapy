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

  std::vector< double > tvec = utl::lin_vector( npoints, 1., tau );
  double Mstar = 0;
  for ( std::size_t it = 1; it < npoints; ++it ) {
    double tloc = 0.5 * ( tvec[ it ] + tvec[ it - 1 ] );
    Mstar +=
      ( tvec[ it ] - tvec[ it - 1 ] ) *
      ( 1 - sed::R_chabrier( tau - tloc ) ) *
      ( *this )( tloc );
  }
  // for ( std::size_t it = 1; it < npoints; ++it ) 
  //   Mstar +=
  //     ( tvec[ it ] - tvec[ it - 1 ] ) *
  //     ( 1 - sed::R_chabrier( tau - tvec[ it ] ) )*
  //     ( *this )( tvec[ it ] );
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

// void sed::sfh_base::time_grid ( const double age,
// 				const std::vector< double > & tgrid,
// 				const std::vector< double > & Zgrid ) {

//   // 1) compute star formation history
//   // 2) compute star-metal enrichment history
//   // 3) find metallicity interval for interpolation
//   _psi_grid  = std::vector< double >( tgrid.size(), 0. );
//   _Z_grid    = std::vector< double >( tgrid.size(), 0. );
//   _Zidx_grid = std::vector< std::size_t >( tgrid.size(), 0 );
//   _last_idx  = 0;

//   double time = age - tgrid[ _last_idx ];
//   while ( time > 0. ) {
//     /*(1)*/ _psi_grid[ _last_idx ]  = ( *this )( time );
//     /*(2)*/ _Z_grid[ _last_idx ]    = ( *this ).get_Zstar( time );
//     /*(3)*/ _Zidx_grid[ _last_idx ] = utl::find_low( _Z_grid[ _last_idx ], Zgrid );
//     ++_last_idx; // increment index
//     time = age - tgrid[ _last_idx ]; // update time
//   }
  
//   return;

// }

// // =============================================================================

// void sed::sfh_base::time_grid ( const double age,
// 				const double * const tgrid,
// 				const std::size_t tgrid_size,
// 				const double * const Zgrid,
// 				const std::size_t Zgrid_size ) {

//   // 1) compute star formation history
//   // 2) compute star-metal enrichment history
//   // 3) find metallicity interval for interpolation
//   _psi_grid  = std::vector< double >( tgrid_size, 0. );
//   _Z_grid    = std::vector< double >( tgrid_size, 0. );
//   _Zidx_grid = std::vector< std::size_t >( tgrid_size, 0 );
//   _last_idx  = 0;

//   double time = age - tgrid[ _last_idx ];
//   while ( time > 0. ) {
//     /*(1)*/ _psi_grid[ _last_idx ]  = ( *this )( time );
//     /*(2)*/ _Z_grid[ _last_idx ]    = ( *this ).get_Zstar( time );
//     /*(3)*/ _Zidx_grid[ _last_idx ] = utl::find_low( _Z_grid[ _last_idx ],
// 						     Zgrid, Zgrid_size );
//     ++_last_idx; // increment index
//     time = age - tgrid[ _last_idx ]; // update time
//   }
  
//   return;

// }

// =============================================================================

void sed::sfh_base::time_grid ( const double age,
				const double * const tgrid,
				const std::size_t tgrid_size,
				const double * const Zgrid,
				const std::size_t Zgrid_size,
			        double * const * const out_psigrid,
			        double * const * const out_Zgrid,
			        std::size_t * const * const out_Zidx,
				std::size_t * const out_last_idx ) {

  // 1) compute star formation history
  // 2) compute star-metal enrichment history
  // 3) find metallicity interval for interpolation  
  std::fill( *out_psigrid, *out_psigrid + tgrid_size, 0. );
  std::fill( *out_Zgrid,   *out_Zgrid   + tgrid_size, 0. );
  std::fill( *out_Zidx,    *out_Zidx    + tgrid_size, 0  ); 
  *out_last_idx = 0;

  double time = age - tgrid[ *out_last_idx ];
  while ( time > 0. ) {
    /*(1)*/ (*out_psigrid)[ *out_last_idx ]  = ( *this )( time );
    /*(2)*/ (*out_Zgrid)[ *out_last_idx ]    = ( *this ).get_Zstar( time );
    /*(3)*/ (*out_Zidx)[ *out_last_idx ] = utl::find_low( (*out_Zgrid)[ *out_last_idx ],
							  Zgrid, Zgrid_size );
    ++( *out_last_idx ); // increment index
    time = age - tgrid[ *out_last_idx ]; // update time
  }
  
  return;

}

// =============================================================================
