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


/* SWIG Exceptions */
%exception {
	try {
		$action
	}  catch (Swig::DirectorException &e) {
		SWIG_exception(SWIG_RuntimeError, e.getMessage());
	}  catch (const std::exception& e) {
		SWIG_exception(SWIG_RuntimeError, e.what());
	}  catch (const char *e) {
		SWIG_exception(SWIG_RuntimeError, e);
	}
}

/* Exceptions for Python lists and iterators */
#ifdef SWIG_PYTHON3
%exception __next__ {
#else
%exception next {
#endif
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
	}
	catch (RangeError) {
		SWIG_exception(SWIG_IndexError, "Index out of bounds");
		return NULL;
	}
};

%{
	#include <iostream>
	#include <iomanip>
	#include <vector>
	#include "CRPropa.h"
	using namespace crpropa;
%}


%ignore operator grplinst::MediumDensity*;
%ignore operator grplinst::MediumTemperature*;
%ignore operator grplinst::Flow*;
%ignore operator grplinst::PlasmaInstability*;

%{
	#include "grplinst/Flow.h"
	#include "grplinst/Medium.h"
	#include "grplinst/PlasmaInstability.h"
%}

/* import crpropa in wrapper */
// %import (module="crpropa") "crpropa-builtin.i"  // doesn't run even if builtin option is active in CMake
%import (module="crpropa") "crpropa.i"


/* handler pointers to containers */
%implicitconv crpropa::ref_ptr<grplinst::MediumDensity>;
%implicitconv crpropa::ref_ptr<grplinst::MediumTemperature>;
%implicitconv crpropa::ref_ptr<grplinst::Flow>;
%implicitconv crpropa::ref_ptr<grplinst::PlasmaInstability>;
%template(MediumDensityRefPtr) crpropa::ref_ptr<grplinst::MediumDensity>;
%template(MediumTemperatureRefPtr) crpropa::ref_ptr<grplinst::MediumTemperature>;
%template(FlowRefPtr) crpropa::ref_ptr<grplinst::Flow>;
%template(PlasmaInstabilityRefPtr) crpropa::ref_ptr<grplinst::PlasmaInstability>;
%feature("director") grplinst::MediumDensity;
%feature("director") grplinst::MediumTemperature;
%feature("director") grplinst::Flow;
%feature("director") grplinst::PlasmaInstability;


/* include plugin parts to generate wrappers */
%include "grplinst/Flow.h"
%include "grplinst/Medium.h"
%include "grplinst/PlasmaInstability.h"


/* Hide warnings */
#pragma SWIG nowarn=312,325,361,389,401,508,509




