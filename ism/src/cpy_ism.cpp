// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Internal headers
#include <cpy_utilities.h>
#include <ism.h>
// STL headers
#include <vector>
#include <cstring>
#include <iostream>

extern "C" {
  
  // ========================================================================================
  // ======================================== CPyISM ========================================
  // ========================================================================================

  // Base class ISM
  typedef struct {
    PyObject_HEAD
    sed::ism * ptrObj;
  } CPyISM;

  // Derived class Diffuse Dust
  typedef struct {
    CPyISM super;
    sed::diffuse * ptrObj;
  } CPyDD;

  // Derived class Molecular Cloud
  typedef struct {
    CPyISM super;
    sed::cloud * ptrObj;
  } CPyMC;

  // ========================================================================================

  // Initialize derived object Diffuse-Dust
  static int CPyDD_init ( CPyDD * self ) {

    self->super.ptrObj = self->ptrObj = new sed::diffuse{};
    return 0;
    
  }

  // Initialize derived object Molecular-Cloud
  static int CPyMC_init ( CPyMC * self ) {

    self->super.ptrObj = self->ptrObj = new sed::cloud{};
    return 0;
    
  }
  
  // ========================================================================================

  static void CPyISM_dealloc ( CPyISM * self ) {

    delete self->ptrObj;
    Py_TYPE( self )->tp_free( self );
    return;
    
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
  static PyObject * CPyISM_set_params ( CPyISM * self, PyObject * args ) {

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
  
  static const char DocString_set_temperature[] =
    "Function for setting the temperature of the ISM.\n"
    "\nParameters"
    "\n----------"
    "\ntemp : scalar float\n"
    "\ttemperature of the medium\n"
    "\nReturns"
    "\n-------"
    "\n: None"; 
  static PyObject * CPyISM_set_temperature ( CPyISM * self, PyObject * args ) {

    double temp;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "d", &temp ) ) return NULL;

    /* Call C++ member function */
    self->ptrObj->set_temperature( temp );

    // equivalent to return None
    Py_RETURN_NONE;

  }
  
  // ========================================================================================
  
  static const char DocString_temperature[] =
    "Computes the temperature of the ISM assuming input total energy.\n"
    "\nParameters"
    "\n----------"
    "\nEtot : scalar float\n"
    "\tTotal energy of the ISM\n"
    "\nReturns"
    "\n-------"
    "\nT : scalar float\n"
    "\tTemperature of the given state of ISM.\n"; 
  static PyObject * CPyISM_temperature ( CPyISM * self, PyObject * args ) {

    double Etot;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "d", &Etot ) ) return NULL;

    /* Call C++ member function */
    return PyFloat_FromDouble( self->ptrObj->temperature( Etot ) );

  }

  // ========================================================================================

  static const char DocString_extinction[] =
    "Computes the ISM extinction at given wavelenght.\n"
    "\nParameters"
    "\n----------"
    "\nlambda : array or scalar float\n"
    "\t\n"
    "\nReturns"
    "\n-------"
    "\n : array or scalar float\n"
    "\t\n";
  static PyObject * CPyISM_extinction ( CPyISM * self, PyObject * buf ) {
    
    /* If first argument is scalar return scalar, else return NumPy array */
    if ( PyObject_TypeCheck( buf, &PyFloat_Type ) )
      return PyFloat_FromDouble( self->ptrObj->extinction( PyFloat_AsDouble( buf ) ) );
    else {

      /* Convert Numpy-array to C-array */
      double * lambda;
      std::size_t size;
      if ( NPyArrayToCArray1D< double >( (PyArrayObject*)buf, &lambda, &size ) == -1 ) return NULL;

      /* Call C++ member function */
      double * outarr = new double [ size ];
      for ( unsigned int ii = 0; ii < size; ++ii )
	outarr[ ii ] = self->ptrObj->extinction( lambda[ ii ] );

      /* Clear heap */
      delete [] lambda;

      // return PyArray_SimpleNewFromData( 1, (npy_intp*)&size, NPY_DOUBLE,
      // 					reinterpret_cast< void * >( outarr ) );

      PyObject * ret = PyArray_SimpleNewFromData( 1, (npy_intp*)&size, NPY_DOUBLE,
						  reinterpret_cast< void * >( outarr ) );
      PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
      return ret;
      
    }
    
  }

  // ========================================================================================

  static const char DocString_attenuation[] =
    "Computes the ISM attenuation at given wavelenght.\n"
    "\nParameters"
    "\n----------"
    "\nlambda : array or scalar float\n"
    "\t\n"
    "\nReturns"
    "\n-------"
    "\n : array or scalar float\n"
    "\t\n";
  static PyObject * CPyISM_attenuation ( CPyISM * self, PyObject * buf ) {
    
    /* If first argument is scalar return scalar, else return NumPy array */
    if ( PyObject_TypeCheck( buf, &PyFloat_Type ) ) 
      return PyFloat_FromDouble( self->ptrObj->attenuation( PyFloat_AsDouble( buf ) ) );
    else {

      /* Convert Numpy-array to C-array */
      double * lambda;
      std::size_t size;
      if ( NPyArrayToCArray1D< double >( (PyArrayObject*)buf, &lambda, &size ) == -1 ) return NULL;

      /* Call C++ member function */
      double * outarr = new double [ size ];
      for ( unsigned int ii = 0; ii < size; ++ii )
	outarr[ ii ] = self->ptrObj->attenuation( lambda[ ii ] );

      /* Clear heap */
      delete [] lambda;
      
      // return PyArray_SimpleNewFromData( 1, (npy_intp*)&size, NPY_DOUBLE,
      // 					reinterpret_cast< void * >( outarr ) );

      PyObject * ret = PyArray_SimpleNewFromData( 1, (npy_intp*)&size, NPY_DOUBLE,
						  reinterpret_cast< void * >( outarr ) );
      PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
      return ret;
      
    }
    
  }

  // ========================================================================================

  static const char DocString_emission[] =
    "Computes the ISM emission at given wavelenght.\n"
    "\nParameters"
    "\n----------"
    "\nlambda : array or scalar float\n"
    "\t\n"
    "\nReturns"
    "\n-------"
    "\nL_ISM : array or scalar float\n"
    "\t\n";
  static PyObject * CPyISM_emission ( CPyISM * self, PyObject * buf ) {
    
    /* If first argument is scalar return scalar, else return NumPy array */
    if ( PyObject_TypeCheck( buf, &PyFloat_Type ) ) 
      return PyFloat_FromDouble( self->ptrObj->emission( PyFloat_AsDouble( buf ) ) );
    else {

      /* Convert Numpy-array to C-array */
      double * lambda;
      std::size_t size;
      if ( NPyArrayToCArray1D< double >( (PyArrayObject*)buf, &lambda, &size ) == -1 ) return NULL;

      /* Call C++ member function */
      double * outarr = new double [ size ];
      for ( unsigned int ii = 0; ii < size; ++ii )
	outarr[ ii ] = self->ptrObj->emission( lambda[ ii ] );

      /* Clear heap */
      delete [] lambda;
      
      // return PyArray_SimpleNewFromData( 1, (npy_intp*)&size, NPY_DOUBLE,
      // 					reinterpret_cast< void * >( outarr ) );

      PyObject * ret = PyArray_SimpleNewFromData( 1, (npy_intp*)&size, NPY_DOUBLE,
						  reinterpret_cast< void * >( outarr ) );
      PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
      return ret;
      
    }
    
  }
  
  // ========================================================================================

  static const char DocString_get_AV[] =
    "Returns the extinction value in the visible band.\n";
  static PyObject * CPyISM_get_AV ( CPyISM * self ) {

    return PyFloat_FromDouble( self->ptrObj->get_params()[ 0 ] );
    
  }  
  
  // ========================================================================================
  // Derived Type Diffuse-Dust

  static PyMethodDef CPyDD_Methods[] = {
					{ "set_params",
					  (PyCFunction) CPyISM_set_params,
					  METH_VARARGS,
					  DocString_set_params },
					{ "set_temperature",
					  (PyCFunction) CPyISM_set_temperature,
					  METH_VARARGS,
					  DocString_set_temperature },
					{ "temperature",
					  (PyCFunction) CPyISM_temperature,
					  METH_VARARGS,
					  DocString_temperature },
  					{ "extinction",
  					  (PyCFunction) CPyISM_extinction,
  					  METH_O,
  					  DocString_extinction },
					{ "attenuation",
					  (PyCFunction) CPyISM_attenuation,
					  METH_O,
					  DocString_attenuation },
  					{ "emission",
  					  (PyCFunction) CPyISM_emission,
  					  METH_O,
  					  DocString_emission },
					{ "A_V",
					  (PyCFunction) CPyISM_get_AV,
					  METH_NOARGS,
					  DocString_get_AV },
					// { "",
					//   (PyCFunction) CPyISM_,
					//   METH_VARARGS | METH_KEYWORDS,
					//   DocString_ },
					{NULL, NULL, 0, NULL}
  };

  static PyTypeObject CPyDD_t = { PyVarObject_HEAD_INIT( NULL, 0 )
				  "ISM_core.CDD"   /* tp_name */
  };
  
  // ========================================================================================
  // ======================================== CPyMC =========================================
  // ========================================================================================

  static const char DocString_eta[] =
    "\n"
    "\nParameters"
    "\n----------"
    "\ntau : array or scalar float\n"
    "\t\n"
    "\nReturns"
    "\n-------"
    "\n : array or scalar float\n"
    "\t\n";
  static PyObject * CPyMC_eta ( CPyMC * self, PyObject * buf ) {

    /* If first argument is scalar return scalar, else return NumPy array */
    if ( PyObject_TypeCheck( buf, &PyFloat_Type ) ) 
      return PyFloat_FromDouble( self->ptrObj->eta( PyFloat_AsDouble( buf ) ) );
    else {

      /* Convert Numpy-array to C-array */
      double * tau;
      std::size_t size;
      if ( NPyArrayToCArray1D< double >( (PyArrayObject*)buf, &tau, &size ) == -1 ) return NULL;

      /* Call C++ member function */
      double * outarr = new double [ size ];
      for ( unsigned int ii = 0; ii < size; ++ii )
      	outarr[ ii ] = self->ptrObj->eta( tau[ ii ] );

      /* Clear heap */
      delete [] tau;
      
      // return PyArray_SimpleNewFromData( 1, (npy_intp*)&size, NPY_DOUBLE,
      // 					reinterpret_cast< void * >( outarr ) );

      PyObject * ret = PyArray_SimpleNewFromData( 1, (npy_intp*)&size, NPY_DOUBLE,
						  reinterpret_cast< void * >( outarr ) );
      PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
      return ret;
      
    } 
    
  }

  // ========================================================================================
  // Derived Type Molecular-Clouds

  static PyMethodDef CPyMC_Methods[] = {
  					{ "set_params",
  					  (PyCFunction) CPyISM_set_params,
  					  METH_VARARGS,
  					  DocString_set_params },
  					{ "set_temperature",
  					  (PyCFunction) CPyISM_set_temperature,
  					  METH_VARARGS,
  					  DocString_set_temperature },
					{ "temperature",
					  (PyCFunction) CPyISM_temperature,
					  METH_VARARGS,
					  DocString_temperature },
  					{ "extinction",
  					  (PyCFunction) CPyISM_extinction,
  					  METH_O,
  					  DocString_extinction },
					{ "attenuation",
					  (PyCFunction) CPyISM_attenuation,
					  METH_O,
					  DocString_attenuation },
  					{ "emission",
  					  (PyCFunction) CPyISM_emission,
  					  METH_O,
  					  DocString_emission },
					{ "A_V",
					  (PyCFunction) CPyISM_get_AV,
					  METH_NOARGS,
					  DocString_get_AV },
  					{ "eta",
  					  (PyCFunction) CPyMC_eta,
  					  METH_O,
  					  DocString_eta },
  					// { "",
  					//   (PyCFunction) CPyISM_,
  					//   METH_VARARGS | METH_KEYWORDS,
  					//   DocString_ },
  					{NULL, NULL, 0, NULL}
  };

  static PyTypeObject CPyMC_t = { PyVarObject_HEAD_INIT( NULL, 0 )
  				  "ISM_core.CMC"   /* tp_name */
  };
  
  // ========================================================================================
  // =================================== Non-Type Methods ===================================
  // ========================================================================================
  
  static const char DocString_TotalAttenuation[] =
    "Function for ...\n"
    "\nParameters"
    "\n----------"
    "\n : \n"
    "\nReturns"
    "\n-------"
    "\n : "; 
  static PyObject * method_TotalAttenuation ( PyObject * self, PyObject * args ) {

    PyObject * llBuf, * attDDBuf, * attMCBuf, * etaMCBuf;
    
    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O|OOO", &llBuf, &attDDBuf, &attMCBuf, &etaMCBuf ) ) return NULL;

    /* If first argument is scalar return scalar, else return NumPy array */
    if ( PyObject_TypeCheck( llBuf, &PyFloat_Type) )
      return PyFloat_FromDouble( sed::total_attenuation( PyFloat_AsDouble( llBuf ),
							 PyFloat_AsDouble( attDDBuf ),
							 PyFloat_AsDouble( attMCBuf ),
							 PyFloat_AsDouble( etaMCBuf ) ) );
    else {

      if ( PyArray_Size( llBuf ) != PyArray_Size( attDDBuf ) ||
	   PyArray_Size( attDDBuf ) != PyArray_Size( attMCBuf ) ) {
	PyErr_SetString( PyExc_AttributeError,
			 "Error in input parameters: "
			 "wavelenght-depending arrays must have same size." );
	return NULL;
      }
      
      /* Convert NPy arrays into C-style arrays */
      double * llArr, * attDDArr, * attMCArr, * etaMCArr;
      std::size_t ll_size, tt_size;
      if ( NPyArrayToCArray1D< double >( (PyArrayObject*)llBuf,
					 &llArr, &ll_size ) == -1 ) return NULL;
      if ( NPyArrayToCArray1D< double >( (PyArrayObject*)attDDBuf,
					 &attDDArr, &ll_size ) == -1 ) return NULL;
      if ( NPyArrayToCArray1D< double >( (PyArrayObject*)attMCBuf,
					 &attMCArr, &ll_size ) == -1 ) return NULL;
      if ( NPyArrayToCArray1D< double >( (PyArrayObject*)etaMCBuf,
					 &etaMCArr, &tt_size ) == -1 ) return NULL;

      npy_intp dims = ll_size * tt_size;
      
      /* Call C++ member function */
      double * outarr = new double [ ll_size * tt_size ];
      for ( unsigned int il = 0; il < ll_size; ++il )
	for ( unsigned int it = 0; it < tt_size; ++it )
	  outarr[ il * tt_size + it ] =
	    sed::total_attenuation( llArr[ il ], attDDArr[ il ], attMCArr[ il ], etaMCArr[ it ] );

      delete [] llArr;
      delete [] attDDArr;
      delete [] attMCArr;
      delete [] etaMCArr;

      // return PyArray_SimpleNewFromData( 1, &dims, NPY_DOUBLE,
      // 					reinterpret_cast< void * >( outarr ) );

      PyObject * ret = PyArray_SimpleNewFromData( 1, (npy_intp*)&dims, NPY_DOUBLE,
						  reinterpret_cast< void * >( outarr ) );
      PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
      return ret;
      
    }

  }
    
  // ========================================================================================

  static PyMethodDef CPyISM_NT_Methods[] =
    {
     { "total_attenuation",
       (PyCFunction) method_TotalAttenuation,
       METH_VARARGS,
       DocString_TotalAttenuation },
     { NULL, NULL, 0, NULL }
    };
  
  // ========================================================================================
  // ===================================== BUILD MODULE =====================================
  // ========================================================================================

  static struct PyModuleDef ism_module = {
					  PyModuleDef_HEAD_INIT,
					  "ISM_core",
					  "Python wrap of c++ ISM component implementation.\n"
					  "Build an object of type ism as:\n"
					  ">>> import galapy.ISM_core as cism\n"
					  ">>> dd = cism.CDD()"
					  ">>> mc = cism.CMC()",
					  -1,
					  NULL, NULL, NULL, NULL, NULL				  
  }; /* endPyModuleDef ism_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_ISM_core( void ) {

    /* -------------------- */
    /* Create CPyISM module */
    /* -------------------- */
    PyObject * m;
    m = PyModule_Create( &ism_module );
    if ( m == NULL )
        return NULL;

    /* -------------------------------------- */
    /* Adding new object CDD to CPyISM module */
    /* -------------------------------------- */
    // Describe CDD object
    CPyDD_t.tp_new       = PyType_GenericNew;
    CPyDD_t.tp_basicsize = sizeof( CPyDD );
    CPyDD_t.tp_dealloc   = (destructor) CPyISM_dealloc;
    CPyDD_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    CPyDD_t.tp_doc       = "Diffuse-Dust objects";
    // CPyDD_t.tp_call      = (ternaryfunc) CPyISM_call;
    CPyDD_t.tp_methods   = CPyDD_Methods;
    CPyDD_t.tp_init      = (initproc) CPyDD_init;
    if ( PyType_Ready( &CPyDD_t ) < 0 )
      return NULL;
    Py_INCREF( &CPyDD_t );
    // Add CDD object to the module
    PyModule_AddObject( m, "CDD", (PyObject *)&CPyDD_t );

    /* -------------------------------------- */
    /* Adding new object CMC to CPyISM module */
    /* -------------------------------------- */
    // Describe CMC object
    CPyMC_t.tp_new       = PyType_GenericNew;
    CPyMC_t.tp_basicsize = sizeof( CPyMC );
    CPyMC_t.tp_dealloc   = (destructor) CPyISM_dealloc;
    CPyMC_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    CPyMC_t.tp_doc       = "Molecular-Cloud objects";
    // CPyMC_t.tp_call      = (ternaryfunc) CPyISM_call;
    CPyMC_t.tp_methods   = CPyMC_Methods;
    CPyMC_t.tp_init      = (initproc) CPyMC_init;
    if ( PyType_Ready( &CPyMC_t ) < 0 )
      return NULL;
    Py_INCREF( &CPyMC_t );
    // Add CMC object to the module
    PyModule_AddObject( m, "CMC", (PyObject *)&CPyMC_t );

    /* ---------------------------------------- */
    /* Adding Non-Type methods to CPyISM module */
    /* ---------------------------------------- */
    PyModule_AddFunctions( m, CPyISM_NT_Methods );

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
