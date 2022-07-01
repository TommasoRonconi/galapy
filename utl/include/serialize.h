#ifndef __SERIALIZE_H__
#define __SERIALIZE_H__

// STL includes
#include <vector>
#include <algorithm> // std::copy

//=============================================================================================

// Abstract Base Class
class Serializable {
  
public :

  virtual size_t serialize_size () const = 0;
  virtual char * serialize ( char* dataOut ) const = 0;
  virtual const char * deserialize ( const char* dataIn ) = 0;

}; // endclass Serializable

//=============================================================================================

template < typename POD >
class SerialPOD {

public :
  
  static std::size_t serialize_size ( POD str ) { return sizeof( POD ); }
  
  static char * serialize ( char * target, const POD value ) {

    std::copy( &value, &value + 1, reinterpret_cast< POD * >( target ) );
    return target + serialize_size( value );
    
  }
  
  static const char * deserialize ( const char * source, POD & target ) {

    std::copy( reinterpret_cast< const POD * >( source ),
	       reinterpret_cast< const POD * >( source ) + 1,
	       &target );
    return source + serialize_size( target );
    
  }
  
}; // endclass SerialPOD

//=============================================================================================
//== these two specializations of the above template class might be required for char*-types ==
//== (Note that they require to include <cstring> and some minor mods to align with the rest) =
//=============================================================================================
// template<>
// size_t SerialPOD< char * >::serialize_size( char * str ) {
//   return sizeof( std::size_t ) + strlen( str );
// }
// template<>
// const char* SerialPOD< char * >::deserialize( const char * source, char * & target ) {
//   size_t length;
//   memcpy( &length, source, sizeof( std::size_t ) );
//   memcpy( &target, source + sizeof( std::size_t ), length );
//   return source + sizeof( std::size_t ) + length;
// }
//=============================================================================================
//=============================================================================================
//=============================================================================================

template < typename POD >
class SerialVecPOD {

public :

  static std::size_t serialize_size ( const std::vector< POD > & str ) {

    return sizeof( std::size_t ) + str.size() * sizeof( POD );

  }

  static char * serialize ( char * target, const std::vector< POD > & value ) {

    target = SerialPOD< std::size_t >::serialize( target, value.size() );
    std::copy( value.begin(), value.end(), reinterpret_cast< POD * >( target ) );
    return target + serialize_size( value );

  }

  static const char * deserialize ( const char * source, std::vector< POD > & target ) {

    std::size_t vec_len;
    source = SerialPOD< std::size_t >::deserialize( source, vec_len );
    target.resize( vec_len );
    std::copy( reinterpret_cast< const POD * >( source ),
    	       reinterpret_cast< const POD * >( source + vec_len * sizeof( POD ) ),
    	       target.begin() );
    return source + vec_len * sizeof( POD );

  }

}; // endclass SerialVecPOD

//=============================================================================================

#endif // __SERIALIZE_H__
