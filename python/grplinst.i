%module(directors="1", threads="1", allprotected="1") grplinst

%include "std_string.i"
%include "std_iostream.i"
%include "std_container.i"
%include "exception.i"

%{
#include "CRPropa.h"
#include "GRPlInst.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/* include plugin parts to generate wrappers */
%include "GRPlInst.h"
%include "grplinst/PlasmaInstability.h"

/* Support for C++ vector */
%include "std_vector.i"
%include "stl.i"

/* Hide warnings */
#pragma SWIG nowarn=312,325,361,389,401,508,509




