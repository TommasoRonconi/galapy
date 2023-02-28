// External includes:
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// STL includes
#include <vector>
// Internal includes
#include <pyb11_serialize.h>
#include <interpolation.h>

// ==========================================================================================

namespace py = pybind11;

template class utl::interpolator< utl::lin_interp >;

// ==========================================================================================

PYBIND11_MODULE( interp, m ) {
  
  py::class_< utl::interpolator< utl::lin_interp > >( m, "lin_interp" )
    .def(py::init< const std::vector< double > &, const std::vector< double > & >())
    .def("get_x", &utl::interpolator< utl::lin_interp >::get_xv,
	 "Return the x-array grid\n")
    .def("get_y", &utl::interpolator< utl::lin_interp >::get_fv,
	 "Return the y-array grid\n")
    .def("__call__", py::vectorize(&utl::interpolator< utl::lin_interp >::operator()),
	 "Evaluate the interpolated function", py::arg("x") )
    .def("integrate", &utl::interpolator< utl::lin_interp >::integrate,
	 "Integrate the interpolated function within a specified interval.\n"
	 "\nParameters"
	 "\n----------"
	 "\naa : float"
	 "\n\tinterval lower limit"
	 "\nbb : float"
	 "\n\tinterval upper limit\n"
	 "\nReturns"
	 "\n------"
	 "\n: float"
	 "\n\tapproximate integral", py::arg("aa"), py::arg("bb") ) //;
    .def(py::pickle(
    		    []( const utl::interpolator< utl::lin_interp > &o ) { //__getstate__
    		      return utl::__getstate__< utl::interpolator< utl::lin_interp > >( o ); },
    		    []( const py::bytes b ) { //__setstate__
    		      return utl::__setstate__< utl::interpolator< utl::lin_interp > >( b ); }
		    ) );
  
} // end PYBIND11_MODULE( interp, m )

// ==========================================================================================
