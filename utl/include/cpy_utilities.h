// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// NumPy headers
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarraytypes.h>
#include <numpy/ndarrayobject.h>
// STL headers
#include <iostream>
#include <vector>
#include <algorithm> // std::copy
#include <cstring>
  
// ==================================================================================
// =================================== C-LINKAGE ====================================
// ==================================================================================

extern "C" {
  
  // ==================================================================================

  static int return_1D_array ( PyObject * bufobj,
			       Py_buffer * view,
			       const char * format = "d" ) {
    
    /* Attempt to extract buffer information from it */ 
    if ( PyObject_GetBuffer( bufobj, view,
			     PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT ) == -1 ) return -1;
  
    if ( view->ndim != 1 ) {
      PyErr_SetString( PyExc_TypeError, "Expected a 1-dimensional array" );
      PyBuffer_Release( view );
      return -1;
    }
  
    /* Check the type of items in the array */
    if ( strcmp( view->format, format ) != 0 ) {
      char ret_string[50];
      strcpy( ret_string, "Wrong format, expecting " );
      strcat( ret_string, view->format );
      strcat( ret_string, ", received ");
      strcat( ret_string, format );
      PyErr_SetString( PyExc_TypeError, ret_string );
      PyBuffer_Release( view );
      return -1;
    }

    return 0;

  }

  // ==================================================================================

} // endextern "C"

// ==================================================================================
// =================================== TEMPLATES ====================================
// ==================================================================================

template< class T >
PyObject *  Py_call_fromScalarOrArray ( T * self, PyObject * args ) {

  PyObject * xx_buf;
  
  if ( !PyArg_ParseTuple( args, "O", &xx_buf ) ) {
    return NULL;
  }
  if ( PyObject_TypeCheck( xx_buf, &PyFloat_Type ) ) {
    double xx = PyFloat_AS_DOUBLE( xx_buf );
    return PyFloat_FromDouble( ( *( self->ptrObj ) )( xx ) );
  }
  else {
    Py_buffer xx;
    if ( return_1D_array( xx_buf, &xx, "d" ) == -1 ) {
      return NULL;
    }
    PyObject * pyList = PyList_New( xx.shape[ 0 ] );
    // std::vector< double > xv ( (double*)xx.buf, (double*)xx.buf + xx.shape[ 0 ] );
    // PyBuffer_Release( &xx );
    // for ( std::size_t ii = 0; ii < xv.size(); ++ii )
    // 	if ( PyList_SetItem( pyList, ii,
    // 			     PyFloat_FromDouble( ( *( self->ptrObj ) )( xv[ ii ] ) )
    // 			     ) != 0 ) return NULL;
    for ( long int ii = 0; ii < xx.shape[ 0 ]; ++ii )
      if ( PyList_SetItem( pyList, ii,
			   PyFloat_FromDouble( ( *( self->ptrObj ) )
					       ( *( (double*)xx.buf + ii ) ) )
			   ) != 0 ) return NULL;
    PyBuffer_Release( &xx );
    return pyList;
  }

}

// ==================================================================================

template< class T >
static int NPyArrayToCArray1D ( PyArrayObject * NPyArr, T ** CArr, std::size_t * size ) {

  if ( PyArray_NDIM( NPyArr ) != 1 ) {
    PyErr_SetString( PyExc_AttributeError, "Expecting a 1-dimensional numpy.array" );
    return -1;
  }

  *size = PyArray_DIM( NPyArr, 0 );
  *CArr = new T [ *size ];
  std::copy( reinterpret_cast< T * >( PyArray_DATA( NPyArr ) ),
	     reinterpret_cast< T * >( PyArray_DATA( NPyArr ) ) + (*size),
	     (*CArr) );
  return 0;
    
}

// ==================================================================================

template< class T >
static int NPyArrayToCxxVector1D ( PyArrayObject * NPyArr, std::vector< T > & CxxVec ) {

  if ( PyArray_NDIM( NPyArr ) != 1 ) {
    PyErr_SetString( PyExc_ValueError, "Expected a 1-dimensional numpy.array" );
    return -1;
  }

  std::size_t size = PyArray_DIM( NPyArr, 0 );
  CxxVec = std::vector< T >( reinterpret_cast< T * >( PyArray_DATA( NPyArr ) ),
			     reinterpret_cast< T * >( PyArray_DATA( NPyArr ) ) + size );

  return 0;
  
}

// ==================================================================================

template< class T, int typenum >
static PyObject * CxxVectorToNPyArray1D ( const std::vector< T > & CxxVec ) {

  PyObject * output = NULL;
  npy_intp size = CxxVec.size();
  T * data = new T [ CxxVec.size() ];
  std::copy( CxxVec.begin(), CxxVec.end(), data );

  output = PyArray_SimpleNewFromData( 1, &size, typenum,
				      reinterpret_cast< void * >( data ) );

  return output;
  
}

// ==================================================================================
