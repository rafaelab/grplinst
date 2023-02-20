%module(directors="1", threads="1", allprotected="1") grplinst

%include "exception.i"
%include "std_container.i"
%include "std_iostream.i"
%include "std_string.i"
%include "std_vector.i"
%include "stl.i"

%{
#include "CRPropa.h"
#include "grplinst.h"
%}

%import (module="crpropa") "crpropa.i"


%include "grplinst.h"



/* Ignore list */
%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;


/* Hide warnings */
#pragma SWIG nowarn=312,325,361,389,401,508,509





