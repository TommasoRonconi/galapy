// Python headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Internal headers
#include <cpy_utilities.h>
#include <cpy_serialize.h>
#include <sfh.h>
// STL headers
#include <vector>
#include <cstring>
#include <iostream>

// static const char* PICKLE_VERSION_KEY = "_pickle_version";
// static int PICKLE_VERSION = 1;

// ==========================================================================================

PyObject * _Create_SFHmodelDict () {

  PyObject * dict = PyDict_New();
  if ( PyDict_SetItemString( dict, "insitu",     PyLong_FromLong( 1 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "constant",   PyLong_FromLong( 2 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "delayedexp", PyLong_FromLong( 3 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "lognormal",  PyLong_FromLong( 4 ) ) < 0 ) return NULL;
  if ( PyDict_SetItemString( dict, "burst",      PyLong_FromLong( 5 ) ) < 0 ) return NULL;
  return dict;
  
}

// ==========================================================================================

sed::sfh_base * set_sfh_model ( const int modelID ) {

  sed::sfh_base * ptrObj = NULL;
  
  switch ( modelID ) {
    
  case 1 :
    ptrObj = new sed::sfh_insitu{};
    break;

  case 2 :
    ptrObj = new sed::sfh_constant{};
    break;

  case 3 :
    ptrObj = new sed::sfh_delayedexp{};
    break;

  case 4 :
    ptrObj = new sed::sfh_lognorm{};
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
    return NULL;

  } // endswitch ( modelID )
  
  return ptrObj;
  
}

// ==========================================================================================

extern "C" {

  // ========================================================================================
  // ======================================== CPySFH ========================================
  // ========================================================================================

  typedef struct {
    PyObject_HEAD
    sed::sfh_base * ptrObj;
    int model = 1;
  } CPySFH;

  // ========================================================================================

  static PyObject * sfh_model = _Create_SFHmodelDict();

  // ========================================================================================

  // initialize CPySFH Object
  static const char DocString_CSFHinit[] =
    "CSFH( self, tau_quench, model = 'insitu')\n"
    "--\n\n"
    "Class defining the model of Star Formation History (SFH).\n\n"
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

    /* Get the model ID and set the ptr object accordingly */
    self->model = (int) PyLong_AsLong( PyDict_GetItemString( sfh_model, model ) );
    self->ptrObj = set_sfh_model( self->model );
    if ( !self->ptrObj ) { 
      PyErr_SetString( PyExc_TypeError,
    		       "SFH model not valid. "
    		       "Valid models are: "
    		       "'insitu', 'constant', "
    		       "'delayedexp', 'lognormal', "
    		       "'burst'" );
      return -1;
    }
    self->ptrObj->set_tau_quench( tau_quench );

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

    PyArrayObject * NPyBuf = NULL;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "O!", &PyArray_Type, &NPyBuf ) ) return NULL;
    
    npy_intp size = PyArray_DIM( NPyBuf, 0 );
    double * psi = new double [size];
    self->ptrObj->model( reinterpret_cast< double* >( PyArray_DATA( NPyBuf ) ), psi, size );

    PyObject * ret = PyArray_SimpleNewFromData( 1, &size, NPY_DOUBLE,
						reinterpret_cast< void * >( psi ) );
    PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
    return ret;

  }

  // ========================================================================================

  /* Pickle the object */
  static PyObject * CPySFH___getstate__ ( CPySFH * self, PyObject * Py_UNUSED(ignored) ) {

    PyObject * ret = CPy___getstate__< CPySFH >( self, NULL );
    if ( !ret ) {
      PyErr_SetString( PyExc_TypeError,
		       "Unable to get state from CPySFH object" );
      return NULL;
    }
    if ( PyDict_SetItemString
	 ( ret, "model", PyLong_FromLong( ( long ) self->model ) ) < 0 )
      return NULL;

    return ret;
    
  }
  
  // // ========================================================================================

  /* Un-pickle the object */
  static PyObject * CPySFH___setstate__( CPySFH * self, PyObject * state ) {
    
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
    
    // Read the bytes and check they are not NULL
    PyObject * bytes = PyDict_GetItemString( state, "ptrBytes" );
    if ( bytes == NULL ) {
      /* PyDict_GetItemString does not set any error state so we have to. */
      PyErr_SetString(PyExc_KeyError, "No \"ptrBytes\" in pickled dict.");
      return NULL;
    }
    
    // Read the model ID:
    /* Borrowed reference but no need to incref as we create a C long from it. */
    PyObject *modelID = PyDict_GetItemString( state, "model" );
    if ( modelID == NULL ) {
      /* PyDict_GetItemString does not set any error state so we have to. */
      PyErr_SetString( PyExc_KeyError, "No \"modelID\" in pickled dict." );
      return NULL;
    }
    self->model = ( int ) PyLong_AsLong( modelID );

    // Allocate memory for the SFH object
    self->ptrObj = set_sfh_model( self->model );
    if ( !self->ptrObj ) { 
      PyErr_SetString( PyExc_TypeError,
    		       "SFH model not valid. "
    		       "Valid models are: "
    		       "'insitu', 'constant', "
    		       "'delayedexp', 'lognormal', "
    		       "'burst'" );
      return NULL;
    }

    // Finally deserialized the un-pickled data
    char * data = new char [ len ];
    Py_ssize_t check;
    if ( PyBytes_AsStringAndSize( PyBytes_FromObject( bytes ),
				  &data, &check ) == -1 ) return NULL;
    // if ( PyBytes_AsStringAndSize( PyBytes_FromObject( bytes ),
    // 				  &data, &check ) == -1 ||
    // 	 ( long ) check != len ) return NULL;
    self->ptrObj->deserialize( data );
  
    // still not sure whether the bytes objects takes the ownership of data
    // but if we delete it pickling and unpickling returns a memory error
    // (nonetheless the python doc says that PyBytes_FromString*
    //  makes a copy of the string, if no-memory leaks, fine with the comment)
    // delete [] data;

    Py_RETURN_NONE;
    
  }

  // ========================================================================================
  
  static const char DocString_set_params[] =
    "set_params( self, params )\n"
    "--\n\n"
    "Function for setting the parameters of the model.\n"
    "\nParameters"
    "\n----------"
    "\nparams : 1d-array\n"
    "\tarray containing the list of parameters of the model. "
    "Parameters have to be passed in the correct order:\n\n"
    "\t#. :code:`'insitu': [psi_max, tau_star]`\n"
    "\t#. :code:`'constant': [psi]`\n"
    "\t#. :code:`'delayedexp': [psi_norm, k_shape, tau_star]`\n"
    "\t#. :code:`'lognormal': [psi_norm, sigma_star, tau_star]`\n"
    "\t#. :code:`'burst': [...]`\n"
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
  
  static const char DocString_set_tau_quench[] =
    "set_tau_quench( self, tau_quench )\n"
    "--\n\n"
    "Function for setting the :math:`\\tau_\\text{quench}` parameter.\n"
    "\nParameters"
    "\n----------"
    "\ntau_quench : float (optional)\n"
    "\tnew value of the `tau_quench` parameter, it regulates the age (expressed in years) "
    "at which the star formation is stopped by an abrupt quenching event. "
    "If no value is passed, it is resetted to the default value (:math:`\\tau_\\text{quench}=20\\ \\text{Gyr}`)\n"
    "\nReturns"
    "\n-------"
    "\n: None"; 
  static PyObject * CPySFH_set_tau_quench ( CPySFH * self, PyObject * args ) {

    double tau_quench = 2.e+10;

    /* Parse arguments */
    if ( !PyArg_ParseTuple( args, "|d", &tau_quench ) ) return NULL;

    /* Call C++ member function */
    self->ptrObj->set_tau_quench( tau_quench );

    // equivalent to return None
    Py_RETURN_NONE;

  }

  // ========================================================================================
  
  static const char DocString_Mstar[] =
    "Mstar( self, tau, npoints = 100 )\n"
    "--\n\n"
    "computes the stellar mass at a given age of the galaxy.\n"
    "It approximates the integral:\n"
    "\n.. math::\n"
    "\n\tM_\\ast(\\tau') = \\int_0^{\\tau'}\\text{d}\\tau "
    "\\bigl[1 - \\mathcal{R}_\\text{IMF}(\\tau)\\bigr]\\psi(\\tau)\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tgalaxy age in years.\n"
    "\nnpoints : int\n"
    "\tthinness for approximated integral computation (default is 100)\n"
    "\nReturns"
    "\n-------"
    "\n: float\n"
    "\tthe stellar mass of the galaxy at time :math:`\\tau`\n"; 
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
    "Mdust( self, tau )\n"
    "--\n\n"
    "Returns the dust mass at a given age of the galaxy.\n"
    "For empirical models of star formation (i.e. :code:`const`, "
    ":code:`delayedexp`, :code:`lognorm`) this is a free parameter.\n"
    "For the :code:`insitu` model, the dust mass is given by "
    ":math:`M_\\text{gas}(\\tau)D(\\tau)` where :math:`M_\\text{gas}=\\psi(\\tau)\\tau_\\ast` "
    "and where :math:`D(\\tau)` is the gas mass ratio "
    "(for an analytic expression of this quantity see Pantoni et al. 2019 and Lapi et al. 2020).\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tage of the galaxy in years.\n"
    "\nReturns"
    "\n-------"
    "\n: float \n"
    "\tDust content in solar masses (:math:`M_\\odot`) at give time :math:`\\tau`.\n"; 
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
    "Mgas( self, tau )\n"
    "--\n\n"
    "Returns the gas mass at a given age of the galaxy.\n"
    "For empirical models of star formation (i.e. :code:`const`, "
    ":code:`delayedexp`, :code:`lognorm`) this is given by :math:`M_\\text{gas} = M_\\text{dust}/D` "
    "where :math:`D \\sim 0.01 (Z_\\text{gas}/Z_\\odot)^{-0.85}` is the "
    "dust-to-gas mass ratio, derived from observations."
    "For the :code:`insitu` model, the gas mass is given by"
    " :math:`M_\\text{gas}=\\psi(\\tau)\\tau_\\ast`\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tage of the galaxy in years.\n"
    "\nReturns"
    "\n-------"
    "\n: float \n"
    "\tGas content in solar masses (:math:`M_\\odot`) at give time :math:`\\tau`.\n";
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
    "Zgas( self, tau )\n"
    "--\n\n"
    "Returns the gas metallicity at a given age of the galaxy.\n"
    "For empirical models of star formation (i.e. :code:`const`, "
    ":code:`delayedexp`, :code:`lognorm`) this is a free parameter with "
    ":math:`Z_\\text{gas} = Z_\\ast`.\n"
    "For the :code:`insitu` model, it is instead given by\n"
    "\n.. math::\n"
    "\n\tZ_\\text{gas}=\\dfrac{s y_Z}{s\\gamma-1}"
    "\\biggl[1 - \\dfrac{(s\\gamma-1)x}{e^{(s\\gamma-1)x}-1}\\biggr]\n"
    "\nwhere :math:`x\\equiv\\tau/s\\tau_\\ast` and :math:`y_Z\\approx0.04` "
    "is the metal production yield (including recycling) for a Chabrier IMF.\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tage of the galaxy in years.\n"
    "\nReturns"
    "\n-------"
    "\n: float \n"
    "\tGas absolute metallicity at give time :math:`\\tau`.\n";
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
    "Zstar( self, tau )\n"
    "--\n\n"
    "Returns the stellar metallicity at a given age of the galaxy.\n"
    "For empirical models of star formation (i.e. :code:`const`, "
    ":code:`delayedexp`, :code:`lognorm`) this is a free parameter with "
    ":math:`Z_\\ast = Z_\\text{gas}`.\n"
    "For the :code:`insitu` model, it is instead given by\n"
    "\n.. math::\n"
    "\n\tZ_\\ast=\\dfrac{y_Z}{\\gamma-1}"
    "\\biggl[1 - \\dfrac{s\\gamma}{(s\\gamma-1}"
    "\\dfrac{e^{-x}-e^{-s\\gamma x}[1 + (s\\gamma -1)x]}"
    "{s\\gamma -1 + e^{-s\\gamma x}- s\\gamma e^{-x}}\\biggr]\n"
    "\nwhere :math:`x\\equiv\\tau/s\\tau_\\ast` and :math:`y_Z\\approx0.04` "
    "is the metal production yield (including recycling) for a Chabrier IMF.\n"
    "\nParameters"
    "\n----------"
    "\ntau : float\n"
    "\tage of the galaxy in years.\n"
    "\nReturns"
    "\n-------"
    "\n: float \n"
    "\tStellar absolute metallicity at give time :math:`\\tau`.\n";
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
    PyObject * NPyPsiGrid = NULL, * NPyZGrid = NULL, * NPyZidxGrid = NULL;
    
    if ( !PyArg_ParseTuple( args, "dO!O!",
			    &age,
			    &PyArray_Type, &tgrid_NPyBuf,
			    &PyArray_Type, &Zgrid_NPyBuf ) ) return NULL;

    /* Convert Numpy-array to C-array */
    double * tgrid, * Zgrid;
    std::size_t tgridsize, Zgridsize;
    if ( NPyArrayToCArray1D< double >( (PyArrayObject*)tgrid_NPyBuf,
				       &tgrid, &tgridsize ) == -1 ) return NULL;
    if ( NPyArrayToCArray1D< double >( (PyArrayObject*)Zgrid_NPyBuf,
				       &Zgrid, &Zgridsize ) == -1 ) return NULL;

    /* Allocate memory */
    double * out_psigrid   = new double [ tgridsize ];
    double * out_Zgrid     = new double [ tgridsize ];
    std::size_t * out_Zidx = new std::size_t [ tgridsize ];
    std::size_t last_idx;

    /* Call member functions and convert to NumPy-arrays */
    self->ptrObj->time_grid( age,
			     tgrid, tgridsize,
			     Zgrid, Zgridsize,
			     &out_psigrid,
			     &out_Zgrid,
			     &out_Zidx,
			     &last_idx );
    // 1) psi-grid
    NPyPsiGrid  = PyArray_SimpleNewFromData( 1, (npy_intp*)&tgridsize, NPY_DOUBLE,
    					     reinterpret_cast< void * >( out_psigrid ) );
    PyArray_ENABLEFLAGS((PyArrayObject*) NPyPsiGrid, NPY_ARRAY_OWNDATA);
    // 2) Z-grid
    NPyZGrid    = PyArray_SimpleNewFromData( 1, (npy_intp*)&tgridsize, NPY_DOUBLE,
    					     reinterpret_cast< void * >( out_Zgrid ) );
    PyArray_ENABLEFLAGS((PyArrayObject*) NPyZGrid, NPY_ARRAY_OWNDATA);
    // 3) Z-indices
    NPyZidxGrid = PyArray_SimpleNewFromData( 1, (npy_intp*)&tgridsize, NPY_UINT64,
    					     reinterpret_cast< void * >( out_Zidx ) );
    PyArray_ENABLEFLAGS((PyArrayObject*) NPyZidxGrid, NPY_ARRAY_OWNDATA);
    // 4) Last valid index
    PyObject * BufLast = PyLong_FromSize_t( last_idx );

    /* Clear heap */
    delete [] tgrid;
    delete [] Zgrid;
    
    // Build return tuple
    PyObject * ret = PyTuple_Pack( 4,
				   NPyPsiGrid,
				   NPyZGrid,
				   NPyZidxGrid,
				   BufLast );
    // 'PyTuple_Pack' increases the reference-count to the
    // packed objects, therefore need to decrease the count
    // after having called:
    Py_DECREF( NPyPsiGrid );
    Py_DECREF( NPyZGrid );
    Py_DECREF( NPyZidxGrid );
    Py_DECREF( BufLast );

    // return the tuple
    return ret;
    
  }

  // ========================================================================================
  
  static PyMethodDef CPySFH_Methods[] =
    {
     { "set_params",
       (PyCFunction) CPySFH_set_params,
       METH_VARARGS,
       DocString_set_params },
     { "set_tau_quench",
       (PyCFunction) CPySFH_set_tau_quench,
       METH_VARARGS,
       DocString_set_tau_quench },
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
     { "__getstate__",
       (PyCFunction) CPySFH___getstate__,
       METH_NOARGS,
       "Pickle the Custom object" },
     { "__setstate__",
       (PyCFunction) CPySFH___setstate__,
       METH_O,
       "Un-pickle the Custom object" },
     {NULL, NULL, 0, NULL}
    };

  static PyTypeObject CPySFH_t = { PyVarObject_HEAD_INIT( NULL, 0 )
				   "galapy.SFH_core.CSFH"   /* tp_name */
  };
  
  // ========================================================================================
  // ===================================== BUILD MODULE =====================================
  // ========================================================================================

  static struct PyModuleDef sfh_module = {
					  PyModuleDef_HEAD_INIT,
					  "galapy.SFH_core",
					  "Python wrap of c++ SFH component implementation.\n"
					  "Build an object of type sfh as:\n"
					  "\t>>> import galapy.SFH_core as csfh\n"
					  "\t>>> sfh = galapy.SFH_core.CSFH()",
					  -1,
					  NULL, NULL, NULL, NULL, NULL				  
  }; /* endPyModuleDef sfh_module */

  /* Create the module */
  PyMODINIT_FUNC PyInit_SFH_core( void ) {

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
    CPySFH_t.tp_doc       = DocString_CSFHinit;
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

