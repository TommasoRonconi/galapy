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

namespace sed {

  class syn {

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

    // constant conversion factor wavelenght [ Angstroms ] -> inverse of the frequency [ 1./GHz ]
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

    ~syn () = default;

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
    
  }; // endclass syn

} // endnamespace sed

#endif //__SYN_H__