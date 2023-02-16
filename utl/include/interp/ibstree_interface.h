#ifndef __IBSTREE_INTERFACE__
#define __IBSTREE_INTERFACE__

/// STL includes
#include <memory>

/// internal includes
#include "base_interface.h"
#include "interval_tree.h"

struct IntAcc : Serializable{

  IntAcc () = default;
  virtual ~IntAcc () = default;

  virtual double eval ( const double xx ) const noexcept = 0;
  virtual double integrate ( const double aa, const double bb ) const noexcept = 0;

}; // endstruct IntAcc

struct LinIntAcc : public IntAcc {

  double m;
  double q;
  double integral;

  LinIntAcc () = default;

  LinIntAcc ( const double x1, const double x2,
	      const double y1, const double y2 ) {

    m = ( y2 - y1 ) / ( x2 - x1 );
    q = y1 - m * x1;
    integral = integrate( x1, x2 );
    
  }
  virtual ~LinIntAcc () = default;

  inline virtual double eval ( const double xx ) const noexcept override {

    return m * xx + q;

  }

  inline virtual double integrate ( const double aa, const double bb ) const noexcept override {

    return
      0.5 * m * bb * bb + q * bb -
      0.5 * m * aa * aa - q * aa;

  } 

  // =============================================================================
  // Serialize Object:

  virtual std::size_t serialize_size () const {

    return
      SerialPOD< double >::serialize_size( m ) +
      SerialPOD< double >::serialize_size( q ) +
      SerialPOD< double >::serialize_size( integral );

  }

  virtual char * serialize ( char * data ) const {

    data = SerialPOD< double >::serialize( data, m );
    data = SerialPOD< double >::serialize( data, q );
    data = SerialPOD< double >::serialize( data, integral );
    return data;

  }

  virtual const char * deserialize ( const char * data ) {

    data = SerialPOD< double >::deserialize( data, m );
    data = SerialPOD< double >::deserialize( data, q );
    data = SerialPOD< double >::deserialize( data, integral );
    return data;

  }

  // =============================================================================
  
}; // endstruct LinIntAcc

namespace utl {

  class lin_interp : public base_interface {

  private:

    ibstree< double, LinIntAcc > _T {};
    void _alloc () {

      for ( std::size_t ii = 0; ii < _thinness - 1; ++ii )
	_T.insert( interval< double >{ _xv[ ii ], _xv[ ii +1 ] },
		   LinIntAcc{ _xv[ ii ], _xv[ ii +1 ],
			       _fv[ ii ], _fv[ ii +1 ] } );
      _T.balance();

    }     
    
  public:
    
    lin_interp () = default;
    lin_interp ( std::function< double ( double ) > func,
		 const double x_min, const double x_max,
		 const std::size_t thinness,
		 const std::string interp_type = "linear" )
      : base_interface{ x_min, x_max, thinness } {

      _xv = lin_vector( thinness, x_min, x_max );
      for ( auto && _x : _xv )
	_fv.emplace_back( func( _x ) );

      // [ BST is an overkill when the X-domain is regularly spaced ]
     
    }

    lin_interp( const std::vector< double > & xv,
		const std::vector< double > & fv,
		const std::string interp_type = "linear" )
      : base_interface{ xv.front(), xv.back(), xv.size() } {

      _xv = xv; _fv = fv;
      _alloc();

    }

    
    /// move constructor
    lin_interp ( lin_interp && ii )
      : base_interface{ ii }, _T{ std::move( _T ) } {}
    
    /// copy constructor
    lin_interp ( const lin_interp & ii )
      : base_interface{ ii } {
      
      _xv = ii.get_xv(); _fv = ii.get_fv();
      _alloc();
      
    }

    /// destructor
    virtual ~lin_interp () = default;

    /// move assignment
    lin_interp & operator= ( lin_interp && ii ) noexcept = default;


    /// copy-assignment operator
    lin_interp & operator= ( lin_interp other ) {

      std::swap( _T, other._T );
      other.swap( *this );
      
      return * this;
	
    }    

    double eval ( const double xx ) const noexcept override {

      return _T.find( xx )->value().eval( xx );

    }

    double integrate ( const double aa, const double bb ) const noexcept override {
      
      // find iterators to intervals containing
      // the lower and upper integral limits
      auto it = _T.find( aa );
      auto stop = _T.find( bb );
      // If the limits belong to the same interval
      // perform integration and return
      if ( it == stop ) 
	return it->value().integrate( aa, bb );
      // initialize return value to integral in first interval [aa, x_j)
      double integral = it->value().integrate( aa, it->key().upp() );

      // traverse the tree adding up integral in each interval visited
      while ( ++it != stop ) {
	integral += it->value().integral;
      }

      // sum last contribute in interval [x_k, bb)
      integral += it->value().integrate( it->key().low(), bb );

      // return value
      return integral;

    }

    // =============================================================================
    // Overload arithmetic operators

    /// overload of operator += for same type add
    virtual lin_interp & operator+= ( const lin_interp & rhs ) {
	
      if ( _thinness != rhs._thinness)
	throw utl_err::size_invalid {
	  "Error in addition: right hand side has different size from left hand side!"
	    };

      for ( size_t ii = 0; ii < _thinness; ++ii ) _fv[ ii ] += rhs._fv[ ii ];

      _alloc();
	
      return *this;
	
    }

    /// overload of operator -= for same type subtract
    virtual lin_interp & operator-= ( const lin_interp & rhs ) {

      if ( _thinness != rhs._thinness)
	throw utl_err::size_invalid {
	  "Error in subtraction: right hand side has different size from left hand side!"
	    };

      for ( size_t ii = 0; ii < _thinness; ++ii ) _fv[ ii ] -= rhs._fv[ ii ];
      
      _alloc();
	
      return *this;
	
    }

    /// overload of operator *= for same type mult
    virtual lin_interp & operator*= ( const lin_interp & rhs ) {
	
      if ( _thinness != rhs._thinness)
	throw utl_err::size_invalid {
	  "Error in multiplication: right hand side has different size from left hand side!"
	    };

      for ( size_t ii = 0; ii < _thinness; ++ii ) _fv[ ii ] *= rhs._fv[ ii ];
      
      _alloc();
	
      return *this;
	
    }

    /// overload of operator /= for same type div
    virtual lin_interp & operator/= ( const lin_interp & rhs ) {
	
      if ( _thinness != rhs._thinness)
	throw utl_err::size_invalid {
	  "Error in division: right hand side has different size from left hand side!"
	    };

      for ( size_t ii = 0; ii < _thinness; ++ii ) _fv[ ii ] /= rhs._fv[ ii ];
      
      _alloc();
	
      return *this;
	
    }

    /// overload of operator += for adding a scalar
    virtual lin_interp & operator+= ( const double & rhs ) {

      for ( size_t ii = 0; ii < _thinness; ++ii ) _fv[ ii ] += rhs;
      
      _alloc();
	
      return *this;
      
    }

    /// overload of operator *= for multiplying by a scalar
    virtual lin_interp & operator*= ( const double & rhs ) {

      for ( size_t ii = 0; ii < _thinness; ++ii ) _fv[ ii ] *= rhs;
      
      _alloc();
	
      return *this;
      
    }

    /// overload of operator *= for dividing by a scalar
    friend lin_interp operator/ ( const double & lhs,
				  lin_interp rhs ) {

      for ( size_t ii = 0; ii < rhs._thinness; ++ii ) rhs._fv[ ii ] *= lhs / rhs._fv[ ii ];
      
      rhs._alloc();
	
      return rhs;
            
    }

    // =============================================================================
    // Serialize Object:

    virtual std::size_t serialize_size () const {

      if ( _T.planted() ) 
	return
	  base_interface::serialize_size() +
	  ( _thinness - 1 ) * ( _T.ctop()->key().serialize_size() +
				_T.ctop()->value().serialize_size() );
      else
	return base_interface::serialize_size();

    }

    virtual char * serialize ( char * data ) const {

      data = base_interface::serialize( data );
      std::vector< const node< double, LinIntAcc > * > store;
      store.reserve( _thinness - 1 );
      _T.extract( store );
      for ( auto && _p : store ) {
	data = _p->key().serialize( data );
	data = _p->value().serialize( data );
      }
      return data;

    }

    virtual const char * deserialize ( const char * data ) {

      data = base_interface::deserialize( data );
      interval< double > IT {};
      LinIntAcc LIA {};
      for ( std::size_t ii = 0; ii < _thinness - 1; ++ii ) {
	data = IT.deserialize( data );
	data = LIA.deserialize( data );
	_T.insert( IT, LIA );
      }
      return data;

    }

    // =============================================================================

  }; // endclass lin_interp

} // endnamespace utl

#endif //__IBSTREE_INTERFACE__
