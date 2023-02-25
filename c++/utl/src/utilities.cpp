#include <utilities.h>

//==============================================================================================

std::vector<double> utl::lin_vector ( const size_t nbin, const double min, const double max ) {

  double bin = ( max - min ) / ( nbin - 1 );
  std::vector< double > vec ( nbin );
  for (size_t ii = 0; ii < nbin; ++ii ) vec[ ii ] = min + bin * ii;
  // vec.emplace_back( min + bin * ii );
  
  return vec;

}

//==============================================================================================

std::vector<double> utl::log_vector ( const size_t nbin, const double min, const double max ) {

  double bin = ( std::log( max ) - std::log( min ) )/( nbin - 1 );
  std::vector< double > vec ( nbin );
  for ( size_t ii = 0; ii < nbin; ++ii ) vec[ ii ] = std::exp( std::log( min ) + bin * ii );
  // vec.emplace_back( std::exp( std::log( min ) + bin * ii ) );
  
  return vec;

}

//==============================================================================================
