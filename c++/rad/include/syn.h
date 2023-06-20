/**
 *  @file rad/include/syn.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __SYN_H__
#define __SYN_H__

// STL includes
// #include <vector>

// internal includes
#include <utilities.h>
#include <serialize.h>

namespace sed {

  class syn : Serializable {

  private :

    // ==================================================
    // private variables of the class

    std::vector< double > _1_of_nu;    // 1 / frequency[ GHz ]
    std::vector< double > _f_emission; // constant term of the emission law
    std::vector< double > _f_optdepth; // constant term of the optical depth

    // vector with the running-values of the parameters
    std::vector< double > _paramsrc;

    double _alpha_syn = 0.75;     // [ adimensional ]
    double _nu_self = 0.2;        // [ GHz ]

    // constant conversion factor wavelength [ Angstroms ] -> inverse of the frequency [ 1./GHz ]
    double _Ang2GHz =
      1.e-8 *          // cm to Angstrom
      1.e+9 /          // s^-1 to GHz 
      sed::cnst::cc;   // speed of light [cm * s^-1];
    
    // ==================================================
    // private functions of the class

    double _opt_depth ( const std::size_t il, const double fact ) const noexcept {

      return fact * _paramsrc[ 1 ] * _f_optdepth[ il ];
      
    }

  public :

    syn () = default;

    virtual ~syn () = default;

    syn ( const std::vector< double > & lambda ) noexcept {

      _paramsrc = std::vector< double >( 2 );
      set_params();
      _1_of_nu    = std::vector< double >( lambda.size() );
      _f_emission = std::vector< double >( lambda.size() );
      _f_optdepth = std::vector< double >( lambda.size() );
      int irc = 0;
      for ( auto && _l : lambda ) {
        _1_of_nu[irc] = _Ang2GHz * _l;
	_f_emission[irc] = 1. / ( 1. + 1. / std::sqrt( 20. * _1_of_nu[irc] ) );
	_f_optdepth[irc] = _1_of_nu[irc] * _1_of_nu[irc] * std::sqrt( _1_of_nu[irc] );
	++irc;
      }

    }
    
    // In:
    // param[ 0 ] = alpha_syn
    // param[ 1 ] = nu_self
    // Out:
    // idx_0 = alpha_syn
    // idx_1 = nu_self^( alpha_syn + 2.5 )
    void set_params ( const double * const params = nullptr ) noexcept {

      if ( params ) {
	_alpha_syn = params[ 0 ];
	_nu_self   = params[ 1 ];
      }
      _paramsrc[ 0 ] = _alpha_syn;
      _paramsrc[ 1 ] = std::pow( _nu_self, _alpha_syn + 2.5 );

    }
    
    double opt_depth ( const std::size_t il ) const noexcept {
      
      double l_tothe_asyn = std::pow( _1_of_nu[ il ], _paramsrc[ 0 ] );
      return _opt_depth( il, l_tothe_asyn );

    }

    double energy ( const std::size_t il, const double fact ) const noexcept {

      double l_tothe_asyn = std::pow( _1_of_nu[ il ], _paramsrc[ 0 ] );
      double tau = _opt_depth( il, l_tothe_asyn );
      return fact *
	l_tothe_asyn *
	_f_emission[ il ] *
	( 1. - std::exp( -tau ) ) / tau;

    }

    // =============================================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      return
	SerialPOD< double >::serialize_size( _alpha_syn ) +
	SerialPOD< double >::serialize_size( _nu_self ) +
	SerialVecPOD< double >::serialize_size( _paramsrc ) +
	SerialVecPOD< double >::serialize_size( _1_of_nu ) +
	SerialVecPOD< double >::serialize_size( _f_emission ) +
	SerialVecPOD< double >::serialize_size( _f_optdepth );

    }

    virtual char * serialize ( char * data ) const {
      
      data = SerialPOD< double >::serialize( data, _alpha_syn );
      data = SerialPOD< double >::serialize( data, _nu_self );
      data = SerialVecPOD< double >::serialize( data, _paramsrc );
      data = SerialVecPOD< double >::serialize( data, _1_of_nu );
      data = SerialVecPOD< double >::serialize( data, _f_emission );
      data = SerialVecPOD< double >::serialize( data, _f_optdepth );
      return data;

    }

    virtual const char * deserialize ( const char * data ) {
      
      data = SerialPOD< double >::deserialize( data, _alpha_syn );
      data = SerialPOD< double >::deserialize( data, _nu_self );
      data = SerialVecPOD< double >::deserialize( data, _paramsrc );
      data = SerialVecPOD< double >::deserialize( data, _1_of_nu );
      data = SerialVecPOD< double >::deserialize( data, _f_emission );
      data = SerialVecPOD< double >::deserialize( data, _f_optdepth );
      return data;
      
    }

    // =============================================================================
    
  }; // endclass syn

} // endnamespace sed

#endif //__SYN_H__
