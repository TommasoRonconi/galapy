// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// STL headers
#include <iostream>
#include <cstring>
#include <vector>

static const char* PICKLE_VERSION_KEY = "_pickle_version";
static int PICKLE_VERSION = 1;

// ==================================================================================
// =================================== TEMPLATES ====================================
// ==================================================================================

/* Pickle the object */
template < class CPyT >
PyObject * CPy___getstate__ ( CPyT * self, PyObject * Py_UNUSED(ignored) ) {

  /* Serialize the object */
  long len = ( long ) self->ptrObj->serialize_size();
  char * data = new char[ len ];
  self->ptrObj->serialize( data );

  /* Convert the serialized object into bytes */
  PyObject * bytes = PyBytes_FromStringAndSize( data, len );
  if ( bytes == NULL ) return NULL;

  /* Build the picklable dictionary */
  PyObject * ret = Py_BuildValue( "{slsOsi}",
				  "bytesLen", len,
				  "ptrBytes", bytes,
				  PICKLE_VERSION_KEY, PICKLE_VERSION );
  
  // still not sure whether the bytes objects takes the ownership of data
  // but if we delete it pickling and unpickling returns a memory error
  // (nonetheless the python doc says that PyBytes_FromString*
  //  makes a copy of the string, if no-memory leaks, fine with the comment)
  // delete [] data;
    
  return ret;
    
}
  
// ==================================================================================

/* Un-pickle the object */
template < class CPyT, class T >
PyObject * CPy___setstate__( CPyT * self, PyObject * state ) {

  /* ---------------------------------- */
  /* Here we make various sanity checks */
  /* ---------------------------------- */

  /* Error check. */
  if ( !PyDict_CheckExact( state ) ) {
    PyErr_SetString( PyExc_ValueError, "Pickled object is not a dict." );
    return NULL;
  }
    
  /* Version check. */
  /* Borrowed reference but no need to increment as we create a C long
   * from it. */
  PyObject * temp = PyDict_GetItemString( state, PICKLE_VERSION_KEY );
  if ( temp == NULL ) {
    /* PyDict_GetItemString does not set any error state so we have to. */
    PyErr_Format( PyExc_KeyError, "No \"%s\" in pickled dict.",
  		  PICKLE_VERSION_KEY );
    return NULL;
  }
  int pickle_version = ( int ) PyLong_AsLong( temp );
  if ( pickle_version != PICKLE_VERSION ) {
    PyErr_Format( PyExc_ValueError,
  		  "Pickle version mismatch. Got version %d but expected version %d.",
  		  pickle_version, PICKLE_VERSION );
    return NULL;
  }

  /* ---------------------------------- */
  /* Here set the custom object members */
  /* ---------------------------------- */

  // First get the bytes lenght
  /* Borrowed reference but no need to incref as we create a C long from it. */
  PyObject * bLen = PyDict_GetItemString( state, "bytesLen" );
  if ( bLen == NULL ) {
    /* PyDict_GetItemString does not set any error state so we have to. */
    PyErr_Format( PyExc_KeyError, "No \"%s\" in pickled dict.",
		  "bytesLen" );
    return NULL;
  }
  long len = PyLong_AsLong( bLen );
    
  // For a Generic PyObject
  PyObject * bytes = PyDict_GetItemString( state, "ptrBytes" );
  if ( bytes == NULL ) {
    /* PyDict_GetItemString does not set any error state so we have to. */
    PyErr_SetString(PyExc_KeyError, "No \"ptrBytes\" in pickled dict.");
    return NULL;
  }
  self->ptrObj = new T {};
  char * data = new char [ len ];
  Py_ssize_t check;
  if ( PyBytes_AsStringAndSize( PyBytes_FromObject( bytes ),
				&data, &check ) == -1 ) return NULL;
  // Might want to use this one instead, in order to
  // check the bytes read are of the right lenght:
  // if ( PyBytes_AsStringAndSize( PyBytes_FromObject( bytes ),
  // 				  &data, &check ) == -1 ||
  // 	 ( long ) check != len ) return NULL;
  // NOTE that the '||' operator should first check the 1st,
  // if false, check the 2nd.
  self->ptrObj->deserialize( data );
  
  // still not sure whether the bytes objects takes the ownership of data
  // but if we delete it pickling and unpickling returns a memory error
  // (nonetheless the python doc says that PyBytes_FromString*
  //  makes a copy of the string, if no-memory leaks, fine with the comment)
  // delete [] data;

  Py_RETURN_NONE;
}

// ==================================================================================
// ==================================================================================
