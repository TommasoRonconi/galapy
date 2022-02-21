/**
 *  @file csp/include/csp.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __CSP_H__
#define __CSP_H__

// external includes
#include <vector>

// internal includes
#include <utilities.h>

namespace sed {

  class csp {

  private :

    using vect_double = std::vector< double >;

    // =============================================================================
    // Private variables

    // SSP variables
    vect_double _lambda;
    std::size_t _Nlambda;
    vect_double _tau;
    std::size_t _Ntau;
    vect_double _Z;
    std::size_t _NZ;
    vect_double _LltZ;
    std::size_t _NL;

    // Parameters
    vect_double _psi, _Zstar;
    std::vector< std::size_t > _iz_low;
    std::size_t _it_last;

    // =============================================================================
    // Private functions

    inline std::size_t
    _Lidx ( const std::size_t il,
	    const std::size_t it,
	    const std::size_t iz )
      const noexcept { return il * _NZ * _Ntau + it * _NZ + iz; }
    
    inline double
    _sfh_lum_timeintegrand ( const std::size_t il,
			     const std::size_t it )
      const noexcept { return				\
	( _tau[ it ] - _tau[ it - 1 ] ) *		\
	_psi[ it ] *					\
	utl::line_from_2points( _Zstar[ it ],
				_Z[ _iz_low[ it ] ],
				luminosity( il, it, _iz_low[ it ] ),
				_Z[ _iz_low[ it ] + 1 ],
				luminosity( il, it, _iz_low[ it ] + 1 ) ); }
    
  public :

    // =============================================================================
    // C.tor/D.tor
    
    csp () = default;

    csp ( const vect_double lambda,
	  const vect_double tau,
	  const vect_double Z,
	  const vect_double LltZ ) :
      _lambda{ lambda }, _Nlambda{ lambda.size() },
      _tau{ tau }, _Ntau{ tau.size() },
      _Z{ Z }, _NZ{ Z.size() },
      _LltZ{ LltZ }, _NL{ LltZ.size() } {}

    ~csp () = default;
    
    // =============================================================================
    // Public functions

    void set_params ( const vect_double & psi,
		      const vect_double & Zstar,
		      const std::vector< std::size_t > & _iz_low,
		      const std::size_t _it_last ) noexcept;

    vect_double get_lambda () const noexcept { return _lambda; }
    vect_double get_tau    () const noexcept { return _tau; }
    vect_double get_Z      () const noexcept { return _Z; }
    
    std::size_t lambda_size () const noexcept { return _Nlambda; }
    std::size_t tau_size    () const noexcept { return _Ntau; }
    std::size_t Z_size      () const noexcept { return _NZ; }
    std::size_t size        () const noexcept { return _NL; }

    inline double
    luminosity ( const std::size_t il,
		 const std::size_t it,
		 const std::size_t iz )
      const noexcept { return _LltZ[ _Lidx( il, it, iz ) ]; }

    double emission ( const std::size_t il,
		      const double * const Tfact ) const noexcept;
    
  }; // endclass csp

} // endnamespace sed

#endif //__CSP_H__
