%module(directors="1", threads="1", allprotected="1") grplinst


%include "attribute.i"
%include "exception.i"
%include "pyabc.i"
%include "stdint.i"
%include "std_array.i"
%include "std_container.i"
%include "std_iostream.i"
%include "std_list.i"
%include "std_map.i"
%include "std_set.i"
%include "std_shared_ptr.i"
%include "std_string.i"
%include "std_vector.i"
%include "stl.i"
%include "typemaps.i"



/* SWIG exceptions */
%inline %{
	class RangeError {
	};
	class StopIterator {
	};
%}


/* C++ containers */
%template() std::pair<double, double>;
%template(PairVector) std::vector<std::pair<double, double>>;




/*************************************************************************************************/
/**	                        			 NumPy interface  		                                **/
/*************************************************************************************************/

%template(VectorInt) std::vector<int>;
%template(VectorFloat) std::vector<float>;
%template(VectorDouble) std::vector<double>;
%template(VectorString) std::vector<std::string>;

%feature("director:except") {
	if ($error != NULL) {
		PyObject* ptype;
		PyObject* pvalue; 
		PyObject* ptraceback;
		PyErr_Fetch(&ptype, &pvalue, &ptraceback);
		PyErr_Restore(ptype, pvalue, ptraceback);
		PyErr_Print();
		Py_Exit(1);
	}
}

/* Exceptions for Python lists and iterators */
%exception __next__ {
	try {
		$action
	} catch (StopIterator) {
		PyErr_SetString(PyExc_StopIteration, "End of iterator");
		return NULL;
	}
}

%exception __getitem__ {
	try {
		$action
	} catch (RangeError) {
		SWIG_exception(SWIG_IndexError, "Index out of bounds");
		return NULL;
	}
};

%{
	#define SWIG_FILE_WITH_INIT
	#include <vector>
	#include <numpy/arrayobject.h>
	#include "numpy/ufuncobject.h"
%}
%include "numpy.i"
%init %{
	import_array();
	import_ufunc();
%}

// typemap for converting std::vector<double> to numpy array
%typemap(out) const std::vector<double>& {
	npy_intp size = $1.size();

	PyObject* obj = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
	if (! obj) {
		SWIG_exception_fail(SWIG_RuntimeError, "Unable to create numpy array");
	}

	double* data = static_cast<double*>(PyArray_DATA((PyArrayObject*) obj));
	for (npy_intp i = 0; i < size; ++i) {
		data[i] = $1[i];
	}

	$result = obj;
}

// typemap for converting numpy array to std::vector<double>
%typemap(in) (std::vector<double>& vec) {
	PyArrayObject* array = (PyArrayObject*) PyArray_FROM_OTF($input, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	if (! array) {
		SWIG_exception_fail(SWIG_TypeError, "Expected a numpy array of type float64");
	}

	npy_intp size = PyArray_SIZE(array);
	double* data = static_cast<double*>(PyArray_DATA(array));
	vec.assign(data, data + size)
	Py_DECREF(array);

	$1 = vec;
}

%apply(double* INPLACE_ARRAY1, int DIM1) { 
	(double* c, int len_c) 
};

%apply(double* ARGOUT_ARRAY1, int DIM1) {
	(double* rangevec, int n)
};

%apply unsigned int &OUTPUT {
	unsigned int&
}; 

/* Python slots */
%feature("python:slot", "sq_length", functype = "lenfunc") __len__;
%feature("python:slot", "mp_subscript", functype = "binaryfunc") __getitem__;
%feature("python:slot", "tp_iter", functype = "unaryfunc") __iter__;
%feature("python:slot", "tp_iternext", functype = "iternextfunc") __next__;

%typemap(directorin, numinputs = 1) (const double* v) {
	npy_intp dim = v.size();
	$input = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void*) $1);
}

%typemap(directorin, numinputs = 1) (const std::vector<double>& v) {
	npy_intp dim = v.size();
	$input = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void*) $1);
}

%fragment("NumPy_Fragments");
%fragment("NumPy_Macros");
%numpy_typemaps(double, NPY_DOUBLE, size_t)
%numpy_typemaps(int, NPY_INT, size_t)




/*************************************************************************************************/
/**                                         CRPropa                                             **/
/*************************************************************************************************/

/* Headers */
%{
	#include "CRPropa.h"
%}

/* Import CRPropa in wrapper */
%import (module = "crpropa") "crpropa.i"


/*************************************************************************************************/
/**                                   grplinst preamble                                         **/
/*************************************************************************************************/

%{
	#include "grplinst.h"
	using namespace grplinst;
%}


/* Include plugin parts to generate wrappers  */
%include "grplinst/Common.h"
%include "grplinst/Flow.h"
%include "grplinst/Geometry.h"
%include "grplinst/Medium.h"
%include "grplinst/PlasmaInstability.h"


/*************************************************************************************************/
/**                          MediumDensity & MediumTemperature                                  **/
/*************************************************************************************************/

%ignore operator grplinst::MediumDensity*;
%ignore operator grplinst::MediumTemperature*;

%implicitconv crpropa::ref_ptr<grplinst::MediumDensity>;
%implicitconv crpropa::ref_ptr<grplinst::MediumTemperature>;

%template(MediumDensityRefPtr) crpropa::ref_ptr<grplinst::MediumDensity>;
%template(MediumTemperatureRefPtr) crpropa::ref_ptr<grplinst::MediumTemperature>;

%feature("director") grplinst::MediumDensity;
%feature("director") grplinst::MediumTemperature;


/*************************************************************************************************/
/**                          				Geometry                           					**/
/*************************************************************************************************/

%ignore operator grplinst::EmissionGeometry*;
%implicitconv crpropa::ref_ptr<grplinst::EmissionGeometry>;
%template(EmissionGeometryRefPtr) crpropa::ref_ptr<grplinst::EmissionGeometry>;


/* provides access to the concrete geometry from the abstract EmissionGeometry object */
%inline %{
	grplinst::Cone* convertToCone(grplinst::EmissionGeometry* geo) {
		return dynamic_cast<grplinst::Cone*>(geo);
	}

	grplinst::Cone* convertToCone(crpropa::ref_ptr<grplinst::EmissionGeometry> geo) {
		return dynamic_cast<grplinst::Cone*>(geo.get());
	}
%}

%feature("director") grplinst::EmissionGeometry;


/*************************************************************************************************/
/**                          				Flow                                 				**/
/*************************************************************************************************/

%ignore operator grplinst::Flow*;
%implicitconv crpropa::ref_ptr<grplinst::Flow>;
%template(FlowRefPtr) crpropa::ref_ptr<grplinst::Flow>;
%feature("director") grplinst::Flow;


/*************************************************************************************************/
/**                          		PlasmaInstability                                			**/
/*************************************************************************************************/

%ignore operator grplinst::PlasmaInstability*;
%implicitconv crpropa::ref_ptr<grplinst::PlasmaInstability>;
%template(PlasmaInstabilityRefPtr) crpropa::ref_ptr<grplinst::PlasmaInstability>;
%feature("director") grplinst::PlasmaInstability;


/*************************************************************************************************/
/*************************************************************************************************/


%clear(double* vector, int length);

/* ignore list */
%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;


/* hide warnings */
#pragma SWIG nowarn=302,312,315,325,361,389,401,508,509






