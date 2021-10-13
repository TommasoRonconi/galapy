/**
 *  @file include/utilities.h
 *
 *  @brief Utilities header containing support functions 
 *
 *  This file declares a set of support functions useful for
 *  several purposes.
 *  It implements an interface for constants, conversion factors
 *  integration methods and random selections
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */


#ifndef __UTILITIES__
#define __UTILITIES__

// STL includes
#include <vector>
#include <cmath>

namespace sed {

  namespace cnst {
    
    /**
     * @name some static constants
     *
     * @{
     */

    /// value of \f$\pi\f$
    static const double pi = 3.1415926535897932;

    /// value of \f$\pi^{-1}\f$
    static const double ip = 1/3.1415926535897932;

    /// value of \f$e^{ 1 }\f$
    static const double e_1 = std::exp( 1 );

    /// value of \f$e^{ - 10 }\f$
    static const double pxe_10 = std::exp( -10 );

    /// value of \f$\ln( 10 )\f$
    static const double ln_10 = std::log( 10 );

    /// value of \f$[\ln( 10 )]^{-1}\f$
    static const double nl_10 = 1. / ln_10;

    /// value of the speed of light \f$c\f$ in \f$[ cm / s ]\f$
    static const double cc = 2.99792458e+10;

    /// value of 1 parsec in centimeters \f$[ cm ]\f$
    static const double pc = 3.08567758130573e+18;

    /// value of 1 year in seconds \f$[ s ]\f$
    static const double yr = 3600 * 24 * 365.256363004;

    /// value of the newtonian gravitational constant \f$[ m^3 kg^{ -1 } s^{ -2 } ]\f$
    static const double Gn = 6.6738480e-11;

    /// inverse of the newtonian gravitational constant \f$[ m^{ -3 } kg s^2 ]\f$
    static const double nG = 1.49838594e+10;

    /// mass of the sun in grams
    static const double Msol = 1.98855e+33;

    /// inverse of the sun mass in grams
    static const double solM = 1./Msol;

    /// solar metallicity mass fraction
    static const double Zsol = 0.0134;

    /// inverse of the solar metallicity mass fraction
    static const double solZ = 1./Zsol;

    /// solar luminosity in \f$[ erg s^{-1} ]\f$
    static const double Lsol = 3.828e33;

    /// inverse of the solar luminosity \f$[ erg^{-1} s ]\f$
    static const double solL = 1./Lsol;

    /// Planck constant \f$[ erg \cdot s ]\f$
    static const double hP = 6.626196e-27;

    /// Boltzmann constant \f$[ erg \cdot K^{-1} ]\f$
    static const double kB = 1.380622e-16;

    /// Used in the Black-Body law \f$2 h_P \cdot c^2 [ erg cm^2 s^{-1} ]\f$
    static const double _2hPc2 = 2 * hP * cc * cc;

    /// Used in the Black-Body law \f$h_P / k_B [ cm \cdot K ]\f$
    static const double _hPccBk = hP * cc / kB;
  
    /// @} End of static constant

  } // endnamespace cnst
  
} //endnamespace sed

namespace utl {

  /// type-safe C++ implementation of the sign (a.k.a. signum) function
  /// from first comment:
  /// stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
  template < typename T >
  int sgn ( T val ) {
    return ( T( 0 ) < val ) - ( val < T( 0 ) );
  }

  template < typename T, typename V >
  T heaviside ( V x, V delay ) {
    return 0.5 * ( 1 + sgn( x - delay ) );
  }

  /**
   * @brief Function that returns a linearly spaced vector in the range [min,max]
   *
   * @param nbin number of sub-intervals in which to divide the interval
   * 
   * @param min minimum of the interval
   *
   * @param max maximum of the interval
   *
   * @return a vector with nbin linearly spaced numbers between min and max
   */
  template < typename T >
  std::vector< T > lin_vector ( const size_t nbin, const T min, const T max ) {
    
    T bin = ( max - min ) / ( nbin - 1 );
    std::vector< T > vec ( nbin );
    for ( std::size_t ii = 0; ii < nbin; ++ii )
      vec[ ii ] = min + bin * ii;
  
    return vec;

  }

  /**
   * @brief Function that returns a logarithmically spaced vector in the range [min, max]
   *	                                                                     
   * @param nbin number of sub-intervals in which to divide the interval
   *	                                                                     
   * @param min minimum of the interval
   *	                                                                     
   * @param max maximum of the interval
   *	                                                                     
   * @return a vector with nbin logarithmically spaced numbers between min and max<br>
   *         <b>Note:</b> the interval [ln(min), ln(max)] is divided in nbin linearly
   *         spaced numbers
   */
  template < typename T >
  std::vector< T > log_vector ( const size_t nbin, const T min, const T max ) {

    T bin = ( std::log( max ) - std::log( min ) )/( nbin - 1 );
    std::vector< T > vec ( nbin );
    for ( std::size_t ii = 0; ii < nbin; ++ii )
      vec[ ii ] = std::exp( std::log( min ) + bin * ii );
  
    return vec;

  }

  /**
   * @brief Naive (a.k.a. brute force) implementation of search within sorted vector of the 
   *        lower index for interpolation/extrapolation between two values
   *
   * @param val value to search for
   * @param vec sorted vector
   * @param start_indx starting point in vector (default = 0, first element of vector)
   * 
   * @result std::size_t index of vector vec.
   *         if val is smaller than the first element, returns the first index (0)
   *         if val is larger than the last element, returns the index before 
   *         the last one ( vec.size()-2 )
   * @warning the algorithm assumes:
   *          - vec.size >= 2,
   *          - vec is sorted in ascending order
   *          - start_indx <= vec.size() 
   *          no exception nor error is raised, 
   *          if the above is not respected, behaviour is undefined
   */  
  inline std::size_t
  find_low ( const double val,
	     const std::vector< double > & vec,
	     const std::size_t start_indx = 0 ) noexcept {
    
    std::size_t indx = start_indx;
    while ( indx < vec.size() - 1 ) {
      if ( val < vec[ indx + 1 ] )
	return indx;
      ++indx;
    }
    return indx - 1;

  }
  
  inline double
  line_from_2points ( const double xx,
		      const double x1, const double y1,
		      const double x2, const double y2 ) noexcept
  { return ( y2 - y1 ) / ( x2 - x1 ) * ( xx - x1 ) + y1; }

  /**
   * @brief Black-Body Emission Law
   *
   * @param ll wavelenght [cm]
   * @param TT temperature [K]
   * 
   * @result
   */
  inline double blackbody ( const double ll, const double TT ) noexcept {

    return sed::cnst::_2hPc2 / ( ll * ll * ll * ll * ll *
				 ( std::exp( sed::cnst::_hPccBk / ( ll * TT ) ) - 1 ) ); 
    
  }

} //endnamespace utl

#endif //__UTILITIES__
