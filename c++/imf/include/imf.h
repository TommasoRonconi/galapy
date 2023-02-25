/**
 *  @file imf/include/imf.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __IMF_H__
#define __IMF_H__

// external includes
#include <cmath>

// internal includes
// #include <utilities.h> // [ external global includes ] 

namespace sed {

  /// Recycled fraction of gas from stellar evolution (for Chabrier IMF) in
  /// the instantaneous recycling approximation
  static const double R_chabrier_inst = 0.45;
  
  /** @brief Recycled fraction of gas from stellar evolution (for Chabrier IMF)
   *
   *  @param tau time must be expressed in units of yr
   * 
   *  @return recycled fraction at given galaxy age 
   */
  inline double R_chabrier ( const double tau ) { return 0.05 * std::log( 1 + 2.5e-6 * tau ); }

} // endnamespace sed

#endif //__IMF_H__
