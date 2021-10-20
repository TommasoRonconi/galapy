/**
 *  @file ibstree/ibstree.tpp
 *
 *  @brief Definition of templated functions of class ibstree
 *
 *  This file defines the interface of the templated functions of the class ibstree
 *  It is included at the end of the header file ibstree/ibstree.h.
 *  Implementing the definition of the class member functions in a separate
 *  file allows us to maintain an ordered structure, while using this strategy
 *  instead of declaring each possible setting of the templated types keeps the
 *  generality of the class.
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

// ===========================================================================


template < class ot, class ou >
std::ostream& operator<< (std::ostream& os, const ibstree< ot, ou >& t) {
  
  const_iterator< ot, ou > it = t.cbegin();
  
  if ( it.operator->() ) {
    const_iterator< ot, ou > stop = t.cend();
    while ( it != stop ) {
      os << it->key() << ":\t" << it->value() << "\n";
      ++it;
    }

  }
  else os << "Empty Tree!";

  return os;
  
}


// ===========================================================================


template < class T, class U >
ibstree< T, U >::ibstree ( const ibstree & T_other ) {

  // Lambda-function making recursive new insertions
  // from top to bottom of hierarchy, starting from some const_iterator
  std::function< void ( const_iterator ) > deep_copy =
    [ this, &deep_copy ] ( const_iterator it ) {
    
    insert( it->key(), it->value() );
    
    if ( it->left )
      deep_copy( const_iterator{ it->left.get() } );
    if ( it->right )
      deep_copy( const_iterator{ it->right.get() } );
    
    return;
    
  };

  deep_copy( const_iterator{ T_other.root.get() } );

}


// ===========================================================================


template < class T, class U >
typename ibstree< T, U >::iterator ibstree< T, U >::insert ( const interval<T> key, const U value ) {

  if ( root ) {
    iterator it { root->insert( key, value ) };
    if ( key < tail->key() ) tail = tail->left.get();
    return it;
  }
  else {
    root.reset( new node{ key, value } );
    tail = root.get();
    return iterator { root.get() };
  }

}

// ===========================================================================


template < class T, class U >
void ibstree< T, U >::balance() {

  // allocate the nodes in vector sorted by key
  std::vector<node*> nodes;
  iterator it = begin();
  for ( ; it != end(); ++it )
    nodes.push_back( new node { it->key(), it->value() } );
  
  // delete un-balanced tree
  clear();

  // if vector not-empty allocate new root and call recursive
  // function to populate new tree, when finished update tail
  if ( nodes.size() > 0 ) {
    root.reset( nodes[ 0.5 * nodes.size() ] );
    kernel_balance( iterator { root.get() }, nodes );
    tail = nodes[ 0 ];
  }
  else {
    tail = root.get();
  }  
 
}

// ===========================================================================


template < class T, class U >
void ibstree< T, U >::kernel_balance( iterator here, const std::vector<node*>& nodes ) {

  auto begin = nodes.begin();
  auto mid = nodes.begin() + 0.5 * nodes.size();
  auto last = nodes.end();
  
  std::vector< node* > left_half { begin, mid };
  if ( left_half.size() > 0 ) {
    here->left.reset( left_half[ 0.5 * left_half.size() ] );
    here->left->parent = here.operator->();
    kernel_balance( iterator{ here->left.get() }, left_half );
  }

  std::vector< node* > right_half { ++mid, last };
  if ( right_half.size() > 0 ) {
    here->right.reset( right_half[ 0.5 * right_half.size() ] );
    here->right->parent = here->parent;
    kernel_balance( iterator{ here->right.get() }, right_half );
  }

  return;
  
}

// ===========================================================================
