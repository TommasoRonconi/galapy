// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Internal headers
#include <cpy_utilities.h>
#include <cpy_serialize.h>
#include <transmission.h>
// STL headers
#include <vector>
#include <cstring>
#include <iostream>

extern "C" {

  // ========================================================================================
  // ======================================== CPyBPT ========================================
  // ========================================================================================

  typedef struct {
    PyObject_HEAD
    utl::transmission * ptrObj;
  } CPyBPT;
  
  // ========================================================================================
  
  // initialize CPyBPT Object
  static const char DocString_BPTinit[] =
    "BPT( self, tau_quench, model = 'insitu')\n"
    "--\n\n"
    "Class defining ...\n\n"
    "\nParameters"
    "\n----------"
    "\n : \n"
    "\t\n"
    "\n : \n"
    "\t\n";
  static int CPyBPT_init ( CPyBPT *self, PyObject *args, PyObject *kwds ) {
    
    PyArrayObject * NPyLL = NULL;
    PyArrayObject * NPyFL = NULL;
    std::vector< double > ll, fl;
    
    /* Get the passed Python objects */
    if ( !PyArg_ParseTuple( args, "O!O!", &PyArray_Type, &NPyLL, &PyArray_Type, &NPyFL ) )
      return -1;
    
    if ( PyArray_DIM( NPyLL, 0 ) != PyArray_DIM( NPyFL, 0 ) ) {
      PyErr_SetString( PyExc_ValueError,
		       "The input arrays must have the same size." );
      return -1;
    }
    
    /* Convert NumPy-arrays to C++ vectors */
    if ( NPyArrayToCxxVector1D< double >( NPyLL, ll ) == -1 ) return -1;
    if ( NPyArrayToCxxVector1D< double >( NPyFL, fl ) == -1 ) return -1;

    self->ptrObj = new utl::transmission{ ll, fl };

    return 0;
    
  }

  // ========================================================================================

  static void CPyBPT_dealloc ( CPyBPT * self ) {

    delete self->ptrObj;
    Py_TYPE( self )->tp_free( self );
    return;
    
  }

  // ========================================================================================

  /* Pickle the object */
  static PyObject * CPyBPT___getstate__ ( CPyBPT * self, PyObject * Py_UNUSED(ignored) ) {

    PyObject * ret = CPy___getstate__< CPyBPT >( self, NULL );
    if ( !ret ) {
      PyErr_SetString( PyExc_TypeError,
		       "Unable to get state from CPyBPT object" );
      return NULL;
    }

    return ret;
    
  }
  
  // ========================================================================================

  /* Un-Pickle the object */
  static PyObject * CPyBPT___setstate__ ( CPyBPT * self, PyObject * state ) {

    if ( !CPy___setstate__< CPyBPT, utl::transmission >( self, state ) ) {
      PyErr_SetString( PyExc_TypeError,
		       "Unable to set state from CPyBPT object" );
      return NULL;
    }
    Py_RETURN_NONE;
    
  }
  
  // ========================================================================================

  static PyObject * CPyBPT_call ( CPyBPT * self, PyObject * args ) {

    PyArrayObject * NPyBuf = NULL;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!", &PyArray_Type, &NPyBuf ) ) return NULL;
    
    npy_intp size = PyArray_DIM( NPyBuf, 0 );
    double * bpt = new double [size];
    for ( long int ii = 0; ii < size; ++ii )
      bpt[ ii ] = ( *( self->ptrObj ) )
	( *reinterpret_cast< double * >( PyArray_GETPTR1( NPyBuf, ii ) ) );

    PyObject * ret = PyArray_SimpleNewFromData( 1, &size, NPY_DOUBLE,
						reinterpret_cast< void * >( bpt ) );
    PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
    return ret;

  }

  // ========================================================================================
  
  // Get the bandpass flux
  static const char DocString_get_bandpass_flux[] =
    "get_bandpass_flux( self, ... )\n"
    "--\n\n"
    "Function computing ...\n\n"
    "\nParameters"
    "\n----------"
    "\n : \n"
    "\t\n"
    "\n : \n"
    "\t\n";
  static PyObject * CPyBPT_get_bandpass_flux ( CPyBPT * self, PyObject * args ) {

    PyArrayObject * NPyLL = NULL;
    PyArrayObject * NPyFL = NULL;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!O!",
			    &PyArray_Type, &NPyLL,
			    &PyArray_Type, &NPyFL ) ) return NULL;

    npy_intp size = PyArray_DIM( NPyLL, 0 );
    if ( size != PyArray_DIM( NPyFL, 0 ) ) {
      PyErr_SetString( PyExc_ValueError,
		       "The input arrays must have the same size." );
      return NULL;
    }
    
    return PyFloat_FromDouble( self->ptrObj->get_bandpass_flux
			      ( reinterpret_cast< double * >( PyArray_DATA( NPyLL ) ),
				reinterpret_cast< double * >( PyArray_DATA( NPyFL ) ),
				size )
			      );

  }

  // ========================================================================================
  
  // Get the bandpass flux
  static const char DocString_get_lmin[] =
    "get_lmin( self, ... )\n"
    "--\n\n"
    "Function computing ...\n\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n";
  static PyObject * CPyBPT_get_lmin ( CPyBPT * self ) {
    
    return PyFloat_FromDouble( self->ptrObj->lmin );

  }

  // ========================================================================================
  
  // Get the bandpass flux
  static const char DocString_get_lmax[] =
    "get_lmax( self, ... )\n"
    "--\n\n"
    "Function computing ...\n\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n";
  static PyObject * CPyBPT_get_lmax ( CPyBPT * self ) {
    
    return PyFloat_FromDouble( self->ptrObj->lmax );

  }

  // ========================================================================================
  
  // Get the bandpass flux
  static const char DocString_get_lpiv[] =
    "get_lpiv( self, ... )\n"
    "--\n\n"
    "Function computing ...\n\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n";
  static PyObject * CPyBPT_get_lpiv ( CPyBPT * self ) {
    
    return PyFloat_FromDouble( self->ptrObj->lpiv );

  }

  // ========================================================================================
  
  // Get the bandpass flux
  static const char DocString_get_norm[] =
    "get_norm( self, ... )\n"
    "--\n\n"
    "Function computing ...\n\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n";
  static PyObject * CPyBPT_get_norm ( CPyBPT * self ) {
    
    return PyFloat_FromDouble( self->ptrObj->norm );

  }

  // ========================================================================================
  
  // Get the X-Axis grid
  static const char DocString_get_xaxis[] =
    "get_xaxis( self, ... )\n"
    "--\n\n"
    "Function computing ...\n\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n";
  static PyObject * CPyBPT_get_xaxis ( CPyBPT * self ) {
    
    return CxxVectorToNPyArray1D< double, NPY_FLOAT64 >( self->ptrObj->get_xaxis() );

  }
  
  // Get the Y-Axis grid
  static const char DocString_get_yaxis[] =
    "get_yaxis( self, ... )\n"
    "--\n\n"
    "Function computing ...\n\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n";
  static PyObject * CPyBPT_get_yaxis ( CPyBPT * self ) {
    
    return CxxVectorToNPyArray1D< double, NPY_FLOAT64 >( self->ptrObj->get_yaxis() );

  }

  // ========================================================================================
  
  static PyMethodDef CPyBPT_Methods[] =
    {
     { "get_bandpass_flux",
       (PyCFunction) CPyBPT_get_bandpass_flux,
       METH_VARARGS,
       DocString_get_bandpass_flux },
     { "get_lmin",
       (PyCFunction) CPyBPT_get_lmin,
       METH_NOARGS,
       DocString_get_lmin },
     { "get_lmax",
       (PyCFunction) CPyBPT_get_lmax,
       METH_NOARGS,
       DocString_get_lmax },
     { "get_lpiv",
       (PyCFunction) CPyBPT_get_lpiv,
       METH_NOARGS,
       DocString_get_lpiv },
     { "get_norm",
       (PyCFunction) CPyBPT_get_norm,
       METH_NOARGS,
       DocString_get_norm },
     { "get_xaxis",
       (PyCFunction) CPyBPT_get_xaxis,
       METH_NOARGS,
       DocString_get_xaxis },
     { "get_yaxis",
       (PyCFunction) CPyBPT_get_yaxis,
       METH_NOARGS,
       DocString_get_yaxis },
     { "__getstate__",
       (PyCFunction) CPyBPT___getstate__,
       METH_NOARGS,
       "Pickle the Custom object" },
     { "__setstate__",
       (PyCFunction) CPyBPT___setstate__,
       METH_O,
       "Un-pickle the Custom object" },
     {NULL, NULL, 0, NULL}
    };

  static PyTypeObject CPyBPT_t = { PyVarObject_HEAD_INIT( NULL, 0 )
				   "galapy.BandpassTransmission.BPT"   /* tp_name */
  };
  
  // ========================================================================================
  // ===================================== BUILD MODULE =====================================
  // ========================================================================================

  static struct PyModuleDef bpt_module = {
					  PyModuleDef_HEAD_INIT,
					  "galapy.BandpassTransmission",
					  "C-Python implementation of "
					  "bandpass transmission (BPT).\n",
					  -1,
					  NULL, NULL, NULL, NULL, NULL				  
  }; /* endPyModuleDef bpt_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_BandpassTransmission( void ) {

    /* ---------------------------------- */
    /* Create BandpassTransmission module */
    /* ---------------------------------- */
    PyObject * m;
    m = PyModule_Create( &bpt_module );
    if ( m == NULL )
        return NULL;

    /* ---------------------------------------------------- */
    /* Adding new object BPT to BandpassTransmission module */
    /* ---------------------------------------------------- */
    // Describe Transmission object
    CPyBPT_t.tp_new       = PyType_GenericNew;
    CPyBPT_t.tp_basicsize = sizeof( CPyBPT );
    CPyBPT_t.tp_dealloc   = (destructor) CPyBPT_dealloc;
    CPyBPT_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    CPyBPT_t.tp_doc       = DocString_BPTinit;
    CPyBPT_t.tp_call      = (ternaryfunc) CPyBPT_call;
    CPyBPT_t.tp_methods   = CPyBPT_Methods;
    //~ CPyBPT_t.tp_members=Noddy_members;
    CPyBPT_t.tp_init      = (initproc) CPyBPT_init;
    if ( PyType_Ready( &CPyBPT_t ) < 0 )
      return NULL;
    Py_INCREF( &CPyBPT_t );
    // Add BPT object to the module
    PyModule_AddObject( m, "BPT", (PyObject *)&CPyBPT_t );

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

