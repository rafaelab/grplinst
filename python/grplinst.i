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

%inline %{
	class RangeError {};
	// class StopIterator {};
%}

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
	#include <cmath>
	#include <fstream>
	#include <iomanip>
	#include <iostream>
	#include <vector>
	#include "CRPropa.h"
	using namespace crpropa;
%}


%ignore operator grplinst::MediumDensity*;
%ignore operator grplinst::MediumTemperature*;
%ignore operator grplinst::Flow*;
%ignore operator grplinst::PlasmaInstability*;
%ignore operator grplinst::Simulation*;
%ignore operator grplinst::Scenario*;
%ignore operator grplinst::EmissionGeometry*;

%{
	#include "grplinst/Auxiliary.h"
	#include "grplinst/Flow.h"
	#include "grplinst/Geometry.h"
	#include "grplinst/Medium.h"
	#include "grplinst/PlasmaInstability.h"
	#include "grplinst/Simulation.h"
%}

/* import crpropa in wrapper */
// %import (module="crpropa") "crpropa-builtin.i"  // doesn't run even if builtin option is active in CMake
%import (module="crpropa") "crpropa.i"


/* handler pointers to containers */
%implicitconv crpropa::ref_ptr<grplinst::MediumDensity>;
%implicitconv crpropa::ref_ptr<grplinst::MediumTemperature>;
%implicitconv crpropa::ref_ptr<grplinst::Flow>;
%implicitconv crpropa::ref_ptr<grplinst::PlasmaInstability>;
%implicitconv crpropa::ref_ptr<grplinst::Scenario>;
%implicitconv crpropa::ref_ptr<grplinst::EmissionGeometry>;
%implicitconv crpropa::ref_ptr<grplinst::EmissionObservables>;
%implicitconv crpropa::ref_ptr<grplinst::Simulation>;
%template(MediumDensityRefPtr) crpropa::ref_ptr<grplinst::MediumDensity>;
%template(MediumTemperatureRefPtr) crpropa::ref_ptr<grplinst::MediumTemperature>;
%template(FlowRefPtr) crpropa::ref_ptr<grplinst::Flow>;
%template(PlasmaInstabilityRefPtr) crpropa::ref_ptr<grplinst::PlasmaInstability>;
%template(ScenarioRefPtr) crpropa::ref_ptr<grplinst::Scenario>;
%template(EmissionGeometryRefPtr) crpropa::ref_ptr<grplinst::EmissionGeometry>;
%template(EmissionObservablesRefPtr) crpropa::ref_ptr<grplinst::EmissionObservables>;
%template(SimulationVector) std::vector<crpropa::ref_ptr<grplinst::Simulation>>;
%template(SimulationRefPtr) crpropa::ref_ptr<grplinst::Simulation>;
%feature("director") grplinst::MediumDensity;
%feature("director") grplinst::MediumTemperature;
%feature("director") grplinst::Flow;
%feature("director") grplinst::PlasmaInstability;
// %feature("director") grplinst::Scenario;
// %feature("director") grplinst::Simulation;
%feature("director") grplinst::EmissionGeometry;
// %feature("director") grplinst::EmissionObservables;


/* provides access to the concrete geometry from the abstract EmissionGeometry object */
%inline %{
	grplinst::Cone* convertToCone(grplinst::EmissionGeometry *geo) {
		return dynamic_cast<grplinst::Cone*>(geo);
	}

	grplinst::Cone* convertToCone(crpropa::ref_ptr<grplinst::EmissionGeometry> geo) {
		return dynamic_cast<grplinst::Cone*>(geo.get());
	}
%}


/* include plugin parts to generate wrappers */
%include "grplinst/Auxiliary.h"
%include "grplinst/Flow.h"
%include "grplinst/Geometry.h"
%include "grplinst/Medium.h"
%include "grplinst/PlasmaInstability.h"
%include "grplinst/Simulation.h"




/* NumPy support */
%{
	#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
	import_array();
%}
%template(vectorDouble) std::vector<double>;
%template(vectorInt) std::vector<int>;



/* Hide warnings */
#pragma SWIG nowarn=312,325,361,389,401,508,509




