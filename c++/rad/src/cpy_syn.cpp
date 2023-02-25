// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Internal headers
#include <cpy_utilities.h>
#include <cpy_serialize.h>
#include <syn.h>
// STL headers
#include <vector>
#include <cstring>
#include <iostream>

extern "C" {
  
  // ========================================================================================
  // ======================================== CPySYN ========================================
  // ========================================================================================

  // Base class SYN
  typedef struct {
    PyObject_HEAD
    sed::syn * ptrObj;
  } CPySYN;

  // ========================================================================================

  // Initialize object Synchrotron
  static int CPySYN_init ( CPySYN * self, PyObject *args ) {

    PyArrayObject * Lbuf = NULL;
    std::vector< double > Lvec;
        
    /* Get the passed Python object */
    if ( !PyArg_ParseTuple( args, "O!", &PyArray_Type, &Lbuf ) ) return -1;
  
    /* Convert NumPy-array to C++ vector */
    if ( NPyArrayToCxxVector1D< double >( Lbuf, Lvec ) == -1 ) return -1;

    /* Call the constructor */
    self->ptrObj = new sed::syn{ Lvec };
    
    return 0;
    
  }
  
  // ========================================================================================

  static void CPySYN_dealloc ( CPySYN * self ) {

    delete self->ptrObj;
    Py_TYPE( self )->tp_free( self );
    return;
    
  }
  
  // ========================================================================================

  /* Pickle the object */
  static PyObject * CPySYN___getstate__ ( CPySYN * self, PyObject * Py_UNUSED(ignored) ) {

    PyObject * ret = CPy___getstate__< CPySYN >( self, NULL );
    if ( !ret ) {
      PyErr_SetString( PyExc_TypeError,
		       "Unable to get state from CPySYN object" );
      return NULL;
    }

    return ret;
    
  }
  
  // ========================================================================================

  /* Un-Pickle the object */
  static PyObject * CPySYN___setstate__ ( CPySYN * self, PyObject * state ) {

    if ( !CPy___setstate__< CPySYN, sed::syn >( self, state ) ) {
      PyErr_SetString( PyExc_TypeError,
		       "Unable to set state from CPySYN object" );
      return NULL;
    }
    Py_RETURN_NONE;
    
  }
  
  // ========================================================================================
  
  static const char DocString_set_params[] =
    "Function for setting the parameters of the model.\n"
    "\nParameters"
    "\n----------"
    "\nparams : array-like\n"
    "\nReturns"
    "\n-------"
    "\n: None"; 
  static PyObject * CPySYN_set_params ( CPySYN * self, PyObject * args ) {

    PyArrayObject * NPyBuf = NULL;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!", &PyArray_Type, &NPyBuf ) ) return NULL;

    /* Convert NumPy-array to C-array */
    double * params;
    std::size_t len;
    if ( NPyArrayToCArray1D< double >( NPyBuf, &params, &len ) == -1 ) return NULL;

    /* Call C++ member function */
    self->ptrObj->set_params( params );

    // equivalent to return None
    Py_RETURN_NONE;

  }
  
  // ========================================================================================
  
  static const char DocString_opt_depth[] =
    "Computes the optical depth of the synchrotron emission at given "
    "wavelenght element in grid.\n"
    "\nParameters"
    "\n----------"
    "\nil : scalar int\n"
    "\tIndex in the wavelenght grid\n"
    "\nReturns"
    "\n-------"
    "\ntau : scalar float\n"
    "\tOptical depth.\n"; 
  static PyObject * CPySYN_opt_depth ( CPySYN * self, PyObject * args ) {

    int il;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "i", &il ) ) return NULL;

    /* Call C++ member function */
    return PyFloat_FromDouble( self->ptrObj->opt_depth( il ) );

  }
  
  // ========================================================================================

    static const char DocString_energy[] =
    "energy( self, il, fact )\n"
    "--\n\n"
    "Computes the synchrotron energy at given index in the wavelenght-grid.\n"
    "\nParameters"
    "\n----------"
    "\nil : array of int\n"
    "\tarray of indexes of the positions in the wavelenght-grid"
    " for which to compute the energy.\n"
    "\nfact : float\n"
    "\tNormalization of the synchrotron spectrum, depends on the radiation source\n"
    "\nReturns"
    "\n-------"
    "\nL_SYN : array or scalar float\n"
    "\tthe self-absorbed Synchrotron energy\n";
  static PyObject * CPySYN_energy ( CPySYN * self, PyObject * args ) {
    
    PyArrayObject * il_buf;
    double fact; 

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!d", &PyArray_Type, &il_buf, &fact ) ) return NULL;

    /* Convert Numpy-array to C-array */
    std::size_t * il_arr;
    std::size_t il_size;
    if ( NPyArrayToCArray1D< std::size_t >( (PyArrayObject*)il_buf, &il_arr, &il_size ) == -1 )
      return NULL;

    double * outarr = new double [ il_size ];
    for ( unsigned int ii = 0; ii < il_size; ++ii ) 
      outarr[ ii ] = self->ptrObj->energy( il_arr[ ii ], fact );

    /* Clear heap */
    delete [] il_arr;

    PyObject * ret = PyArray_SimpleNewFromData( 1, (npy_intp*)&il_size, NPY_DOUBLE,
						reinterpret_cast< void * >( outarr ) );
    PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
    return ret;
    
  }
  
  // ========================================================================================
  
  static PyMethodDef CPySYN_Methods[] =
    {
     { "set_params",
       (PyCFunction) CPySYN_set_params,
       METH_VARARGS,
       DocString_set_params },
     { "opt_depth",
       (PyCFunction) CPySYN_opt_depth,
       METH_VARARGS,
       DocString_opt_depth },
     { "energy",
       (PyCFunction) CPySYN_energy,
       METH_VARARGS,
       DocString_energy },
     { "__getstate__",
       (PyCFunction) CPySYN___getstate__,
       METH_NOARGS,
       "Pickle the Custom object" },
     { "__setstate__",
       (PyCFunction) CPySYN___setstate__,
       METH_O,
       "Un-pickle the Custom object" },
     {NULL, NULL, 0, NULL}
    };

  static PyTypeObject CPySYN_t = { PyVarObject_HEAD_INIT( NULL, 0 )
				   "galapy.SYN_core.CSYN"   /* tp_name */
  };
  
  // ========================================================================================

  static struct PyModuleDef syn_module = {
					  PyModuleDef_HEAD_INIT,
					  "galapy.SYN_core",
					  "Python wrap of c++ SYN component implementation.\n"
					  "Build an object of type syn as:\n"
					  ">>> import galapy.SYN_core as csyn\n"
					  ">>> import numpy as np\n"
					  ">>> ll = np.logspace(1,8,100)\n"
					  ">>> syn = csyn.CSYN( ll )\n",
					  -1,
					  NULL, NULL, NULL, NULL, NULL				  
  }; /* endPyModuleDef syn_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_SYN_core( void ) {

    /* -------------------- */
    /* Create CPySYN module */
    /* -------------------- */
    PyObject * m;
    m = PyModule_Create( &syn_module );
    if ( m == NULL )
        return NULL;

    /* -------------------------------------- */
    /* Adding new object CSYN to CPySYN module */
    /* -------------------------------------- */
    // Describe CSYN object
    CPySYN_t.tp_new       = PyType_GenericNew;
    CPySYN_t.tp_basicsize = sizeof( CPySYN );
    CPySYN_t.tp_dealloc   = (destructor) CPySYN_dealloc;
    CPySYN_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    CPySYN_t.tp_doc       = "Synchrotron object";
    CPySYN_t.tp_methods   = CPySYN_Methods;
    CPySYN_t.tp_init      = (initproc) CPySYN_init;
    if ( PyType_Ready( &CPySYN_t ) < 0 )
      return NULL;
    Py_INCREF( &CPySYN_t );
    // Add CSYN object to the module
    PyModule_AddObject( m, "CSYN", (PyObject *)&CPySYN_t );

    /* -------------------------------- */
    /* Initialize NumPy Array interface */
    /* -------------------------------- */
    import_array();

    /* ----------------- */
    /* return new module */
    /* ----------------- */
    return m;
    
  }
  
  // ========================================================================================

} // endextern "C"
