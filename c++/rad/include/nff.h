/**
 *  @file rad/include/nff.h
 *
 *  @brief Class implementing the Nebular Free-Free emission (i.e. Bremsstrahlung)
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __NFF_H__
#define __NFF_H__

// external includes
// #include <vector>

// internal includes
#include <utilities.h>
#include <serialize.h>

namespace sed {

  class nff : Serializable {

  private :

    // =============================================================================
    // Private variables
    
    // average # of protons in iones (set to 1 for a pure Hydrogen plasma)
    double _lnZi = std::log( 1 );
    double _Zgas = 0.01;

    //
    std::vector< double > _paramsrc;

    // units array
    std::vector< double > _nu;

    // term not parameterization dependent in gaunt factor:
    std::vector< double > _fact_gff;
    
    // for converting emission units [erg s^-1 Hz^-1] -> [erg s^-1 Ang^-1]:
    std::vector< double > _Hz2Ang;
    
    // =============================================================================
    // Private functions
    
  public :

    // =============================================================================
    // C.tor/D.tor
    
    nff () = default;
    
    nff ( const std::vector< double > lambda ) {
      
      _nu.reserve( lambda.size() );
      _Hz2Ang.reserve( lambda.size() );
      _fact_gff.reserve( lambda.size() );
      for ( auto && _l : lambda ) {
	_nu.push_back( 1.e+8 * sed::cnst::cc / _l );
	_Hz2Ang.push_back( _nu.back() / _l );
	_fact_gff.push_back( std::log( 1.e-9 * _nu.back() ) );
      }
      
    }
    
    virtual ~nff () = default;
    
    // =============================================================================
    // Public functions

    // note that the value of lnZgas is defined as
    // log( Zgas / 0.02 )
    inline double lTe ( const double lZgas ) const noexcept {

      return sed::cnst::ln_10 * ( 3.89 - lZgas * ( 0.0205 * lZgas + 0.4802 ) );

    }

    inline double gff ( const std::size_t il, const double lnTe4 ) const noexcept {

      return
	std::log( std::exp( 5.96 - sed::cnst::ip * 1.732050808 *
			    ( _fact_gff[ il ] + _lnZi - 1.5 * lnTe4 ) )
		  + sed::cnst::e_1 );
      
    }
    
    // In:
    // param[ 0 ] = Z_gas
    // param[ 1 ] = Z_i
    // Out:
    // idx_0 = lnTe4
    // idx_1 = 1./Te
    void set_params ( const double * const params ) noexcept {

      _Zgas = params[ 0 ];
      _lnZi = std::log( params[ 1 ] );
      double lnTe = lTe( std::log10( 50 * _Zgas ) );
      _paramsrc = { lnTe - 4 * sed::cnst::ln_10,
		    std::exp( -lnTe ) };

    }

    double emission ( const std::size_t il, const double Q_H ) const noexcept {

      return 1.8e-27 * Q_H * gff( il, _paramsrc[ 0 ] ) * _Hz2Ang[ il ] *
	std::exp( 0.3 * _paramsrc[ 0 ] -
		  sed::cnst::_hPBk * _nu[ il ] * _paramsrc[ 1 ] );
  
    }

    // =============================================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	SerialPOD< double >::serialize_size( _lnZi ) +
	SerialPOD< double >::serialize_size( _Zgas ) +
	SerialVecPOD< double >::serialize_size( _paramsrc ) +
	SerialVecPOD< double >::serialize_size( _nu ) +
	SerialVecPOD< double >::serialize_size( _fact_gff ) +
	SerialVecPOD< double >::serialize_size( _Hz2Ang );

    }

    virtual char * serialize ( char * data ) const {
      
      data = SerialPOD< double >::serialize( data, _lnZi );
      data = SerialPOD< double >::serialize( data, _Zgas );
      data = SerialVecPOD< double >::serialize( data, _paramsrc );
      data = SerialVecPOD< double >::serialize( data, _nu );
      data = SerialVecPOD< double >::serialize( data, _fact_gff );
      data = SerialVecPOD< double >::serialize( data, _Hz2Ang );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {
      
      data = SerialPOD< double >::deserialize( data, _lnZi );
      data = SerialPOD< double >::deserialize( data, _Zgas );
      data = SerialVecPOD< double >::deserialize( data, _paramsrc );
      data = SerialVecPOD< double >::deserialize( data, _nu );
      data = SerialVecPOD< double >::deserialize( data, _fact_gff );
      data = SerialVecPOD< double >::deserialize( data, _Hz2Ang );
      return data;
      
    }

    // =============================================================================
    
  }; // endclass nff

} // endnamespace sed

#endif //__NFF_H__
