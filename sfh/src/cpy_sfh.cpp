// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Internal headers
#include <cpy_utilities.h>
#include <sfh.h>
// STL headers
#include <vector>
#include <cstring>
#include <iostream>

PyObject * _Create_SFHmodelDict () {

  PyObject * dict = PyDict_New();
  if ( PyDict_SetItemString( dict, "insitu",     PyLong_FromLong( 1 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "constant",   PyLong_FromLong( 2 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "delayedexp", PyLong_FromLong( 3 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "lognormal",  PyLong_FromLong( 4 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "burst",      PyLong_FromLong( 5 ) ) < 0 ) return NULL;
  return dict;
  
}

extern "C" {

  // ========================================================================================
  // ======================================== CPySFH ========================================
  // ========================================================================================

  typedef struct {
    PyObject_HEAD
    sed::sfh_base * ptrObj;
  } CPySFH;

  // ========================================================================================

  static PyObject * sfh_model = _Create_SFHmodelDict();

  // ========================================================================================

  // initialize CPySFH Object
  static const char DocString_SFHinit[] =
    "Class defining the model of Star Formation History.\n\n"
    "The possible models to choose are\n\n"
    "#. 'insitu'\n"
    "#. 'constant'\n"
    "#. 'delayedexp'\n"
    "#. 'lognormal'\n"
    "#. 'burst'\n"
    "\nParameters"
    "\n----------"
    "\ntau_quench : float\n"
    "\tEventual abrupt quenching time for star formation.\n"
    "\tShould be expressed in years. Refers to the age of the galaxy.\n"
    "\nmodel : string\n"
    "\tOne among ('insitu', 'constant', 'delayedexp', 'lognormal', 'burst').\n"
    "\tDefault is 'insitu'.\n";
  static int CPySFH_init ( CPySFH *self, PyObject *args, PyObject *kwds ) {
    
    double tau_quench;
    char const * model = "insitu";
    
    // work-around to silence compile warning (const char* -> char* forbidden in C++):
    char ** kwlist = new char * [ 3 ] { (char*)"", (char*)"model", NULL };
    
    /* Get the passed Python object */
    if ( !PyArg_ParseTupleAndKeywords( args, kwds, "d|s", kwlist, &tau_quench, &model ) ) {
      return -1;
    }

    switch ( (int) PyLong_AsLong( PyDict_GetItemString( sfh_model, model ) ) ) {

    case 1 :
      self->ptrObj = new sed::sfh_insitu{ tau_quench };
      break;

    case 2 :
      self->ptrObj = new sed::sfh_constant{ tau_quench };
      break;

    case 3 :
      self->ptrObj = new sed::sfh_delayedexp{ tau_quench };
      break;

    case 4 :
      self->ptrObj = new sed::sfh_lognorm{ tau_quench };
      break;

    // case 5 :
    //   self->ptrObj = new sed::sfh_burst{};
    //   break;

    default :
      PyErr_SetString( PyExc_TypeError,
    		       "SFH model not valid. "
    		       "Valid models are: "
    		       "'insitu', 'constant', "
    		       "'delayedexp', 'lognormal', "
    		       "'burst'" );
      return -1;

    } // endswitch ( model )

    /* Clear heap */
    delete [] kwlist;

    return 0;
    
  }

  // ========================================================================================

  static void CPySFH_dealloc ( CPySFH * self ) {

    delete self->ptrObj;
    Py_TYPE( self )->tp_free( self );
    return;
    
  }

  // ========================================================================================

  static PyObject * CPySFH_call ( CPySFH * self, PyObject * args ) {

    return Py_call_fromScalarOrArray< CPySFH >( self, args );

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
  static PyObject * CPySFH_set_params ( CPySFH * self, PyObject * args ) {

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
  
  static const char DocString_Mstar[] =
    "Function computing the galaxy stellar mass at given time.\n"
    "It approximates the integral:\n"
    "\n.. math::\n\n\tM_\\ast(\\tau') = \\int_0^{\\tau'}\\text{d}\\tau "
    "\\bigl[1 - \\mathcal{R}_\\text{IMF}(\\tau)\\bigr]\\psi(\\tau)\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tgalactic age\n"
    "\nnpoints : int\n"
    "\tthinness for approximated integral computation\n"
    "\nReturns"
    "\n-------"
    "\nM* : float\n"
    "\tthe stellar mass of the galaxy at time `tau`\n"; 
  static PyObject * CPySFH_Mstar ( CPySFH * self, PyObject * args, PyObject * kwds ) {

    double tau;
    Py_ssize_t npoints = 100;
    
    // work-around to silence compile warning (const char* -> char* forbidden in C++):
    char ** kwlist = new char * [ 3 ] { (char*)"", (char*)"npoints", NULL };

    // Parse arguments
    if ( !PyArg_ParseTupleAndKeywords( args, kwds, "d|n", kwlist, &tau, &npoints ) ) {
      return NULL;
    }

    /* Clear heap */
    delete [] kwlist;

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Mstar( tau, (std::size_t)npoints ) );

  }

  // ========================================================================================
  
  static const char DocString_Mdust[] =
    "Function ..\n"
    "\nParameters"
    "\n----------"
    "\n : \n"
    "\t\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n"; 
  static PyObject * CPySFH_Mdust ( CPySFH * self, PyObject * args ) {

    double tau;

    // Parse arguments
    if ( !PyArg_ParseTuple( args, "d", &tau ) ) {
      return NULL;
    }        

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Mdust( tau ) );

  }

  // ========================================================================================
  
  static const char DocString_Mgas[] =
    "Function ..\n"
    "\nParameters"
    "\n----------"
    "\n : \n"
    "\t\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n"; 
  static PyObject * CPySFH_Mgas ( CPySFH * self, PyObject * args ) {

    double tau;

    // Parse arguments
    if ( !PyArg_ParseTuple( args, "d", &tau ) ) {
      return NULL;
    }        

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Mgas( tau ) );

  }

  // ========================================================================================
  
  static const char DocString_Zgas[] =
    "Function ..\n"
    "\nParameters"
    "\n----------"
    "\n : \n"
    "\t\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n"; 
  static PyObject * CPySFH_Zgas ( CPySFH * self, PyObject * args ) {

    double tau;

    // Parse arguments
    if ( !PyArg_ParseTuple( args, "d", &tau ) ) {
      return NULL;
    }        

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Zgas( tau ) );

  }

  // ========================================================================================
  
  static const char DocString_Zstar[] =
    "Function ..\n"
    "\nParameters"
    "\n----------"
    "\n : \n"
    "\t\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n"; 
  static PyObject * CPySFH_Zstar ( CPySFH * self, PyObject * args ) {

    double tau;

    // Parse arguments
    if ( !PyArg_ParseTuple( args, "d", &tau ) ) {
      return NULL;
    }        

    /* Call member function and return PyObject*/ 
    return PyFloat_FromDouble( self->ptrObj->get_Zstar( tau ) );

  }

  // ========================================================================================
  
  static const char DocString_time_grid[] =
    "Function computing ...\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tgalactic age\n"
    "\nnpoints : int\n"
    "\nReturns"
    "\n-------"
    "\n : \n"
    "\t\n"; 
  static PyObject * CPySFH_time_grid ( CPySFH * self, PyObject * args, PyObject * kwds ) {

    double age;
    PyArrayObject * tgrid_NPyBuf = NULL, * Zgrid_NPyBuf = NULL;
    PyArrayObject * NPyPsiGrid = NULL, * NPyZGrid = NULL, * NPyZidxGrid = NULL;
    std::vector< double > tgrid_CxxVec, Zgrid_CxxVec;
    
    if ( !PyArg_ParseTuple( args, "dO!O!",
			    &age,
			    &PyArray_Type, &tgrid_NPyBuf,
			    &PyArray_Type, &Zgrid_NPyBuf ) ) return NULL;

    /* Convert NumPy-array to C-array */
    if ( NPyArrayToCxxVector1D< double >( tgrid_NPyBuf, tgrid_CxxVec ) == -1 ) return NULL;
    if ( NPyArrayToCxxVector1D< double >( Zgrid_NPyBuf, Zgrid_CxxVec ) == -1 ) return NULL;

    /* Call member functions and convert to NumPy-arrays */
    self->ptrObj->time_grid( age, tgrid_CxxVec, Zgrid_CxxVec );
    NPyPsiGrid  =
      ( PyArrayObject * )
      CxxVectorToNPyArray1D< double, NPY_DOUBLE >( self->ptrObj->get_psi_grid() );
    NPyZGrid    =
      ( PyArrayObject * )
      CxxVectorToNPyArray1D< double, NPY_DOUBLE >( self->ptrObj->get_Z_grid() );
    NPyZidxGrid =
      ( PyArrayObject * )
      CxxVectorToNPyArray1D< std::size_t, NPY_UINT64 >( self->ptrObj->get_Zidx_grid() );
    
    // return tuple
    return PyTuple_Pack( 4, NPyPsiGrid, NPyZGrid, NPyZidxGrid,
  			 PyLong_FromSize_t( self->ptrObj->get_last_grid_idx() ) );
    
  }

  // ========================================================================================
  
  static PyMethodDef CPySFH_Methods[] =
    {
     { "set_params",
       (PyCFunction) CPySFH_set_params,
       METH_VARARGS,
       DocString_set_params },
     { "Mstar",
       (PyCFunction) CPySFH_Mstar,
       METH_VARARGS | METH_KEYWORDS,
       DocString_Mstar },
     { "Mdust",
       (PyCFunction) CPySFH_Mdust,
       METH_VARARGS,
       DocString_Mdust },
     { "Mgas",
       (PyCFunction) CPySFH_Mgas,
       METH_VARARGS,
       DocString_Mgas },
     { "Zgas",
       (PyCFunction) CPySFH_Zgas,
       METH_VARARGS,
       DocString_Zgas },
     { "Zstar",
       (PyCFunction) CPySFH_Zstar,
       METH_VARARGS,
       DocString_Zstar },
     { "time_grid",
       (PyCFunction) CPySFH_time_grid,
       METH_VARARGS,
       DocString_time_grid },
     {NULL, NULL, 0, NULL}
    };

  static PyTypeObject CPySFH_t = { PyVarObject_HEAD_INIT( NULL, 0 )
				   "CPySFH.CSFH"   /* tp_name */
  };
  
  // ========================================================================================
  // ===================================== BUILD MODULE =====================================
  // ========================================================================================

  static struct PyModuleDef sfh_module = {
					  PyModuleDef_HEAD_INIT,
					  "CPySFH",
					  "Python wrap of c++ SFH component implementation.\n"
					  "Build an object of type sfh as:\n"
					  ">>> import galapy.internal.CPySFH as csfh\n"
					  ">>> sfh = csfh.CSFH()",
					  -1,
					  NULL, NULL, NULL, NULL, NULL				  
  }; /* endPyModuleDef sfh_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_CPySFH( void ) {

    /* -------------------- */
    /* Create CPySFH module */
    /* -------------------- */
    PyObject * m;
    m = PyModule_Create( &sfh_module );
    if ( m == NULL )
        return NULL;

    /* --------------------------------------- */
    /* Adding new object CSFH to CPySFH module */
    /* --------------------------------------- */
    // Describe CSFH object
    CPySFH_t.tp_new       = PyType_GenericNew;
    CPySFH_t.tp_basicsize = sizeof( CPySFH );
    CPySFH_t.tp_dealloc   = (destructor) CPySFH_dealloc;
    CPySFH_t.tp_flags     = Py_TPFLAGS_DEFAULT;
    CPySFH_t.tp_doc       = DocString_SFHinit;
    CPySFH_t.tp_call      = (ternaryfunc) CPySFH_call;
    CPySFH_t.tp_methods   = CPySFH_Methods;
    //~ CPySFH_t.tp_members=Noddy_members;
    CPySFH_t.tp_init      = (initproc) CPySFH_init;
    if ( PyType_Ready( &CPySFH_t ) < 0 )
      return NULL;
    Py_INCREF( &CPySFH_t );
    // Add CSFH object to the module
    PyModule_AddObject( m, "CSFH", (PyObject *)&CPySFH_t );

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

