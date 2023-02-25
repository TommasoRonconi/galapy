/**
 *  @file ibstree/node.tpp
 *
 *  @brief Definition of templated functions of struct node
 *
 *  This file defines the interface of the templated functions of the struct node
 *  It is included at the end of the header file ibstree/node.h.
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */


// ===========================================================================


template< class T, class U >
node< T, U > * node< T, U >::insert ( const interval< T > key, const U value ) {

  node * n = nullptr;
  
  if( key < content.first ) {
    if( left ) 
      n = left->insert( key, value );
    else {
      left.reset( new node { key, value, this } );
      n = left.get();
    }
  }

  if ( key > content.first ) {
    if ( right ) 
      n = right->insert( key, value );
    else {
      right.reset( new node{ key, value, parent } );
      n = right.get();
    }
  }

  return n;
      
}


// ===========================================================================

template< class T, class U >
void node< T, U >::extract ( std::vector< const node< T, U > * > & store )
  const noexcept {

  store.emplace_back( this );
  if ( left )
    left->extract( store );
  if ( right )
    right->extract( store );
  return;
  
}

// ===========================================================================


template< class T, class U >
node< T, U > * node< T, U >::find ( const T key ) {
  
  if ( key == content.first ) 
    return this;

  if ( ( left ) && ( key < content.first ) ) 
    return left->find( key );
	
  if ( ( right ) && ( key > content.first ) ) 
    return right->find( key );

  // out-of-bounds not allowed:
  // return nullptr;

  // out-of-bounds allowed:
  return this;

}


// ===========================================================================
