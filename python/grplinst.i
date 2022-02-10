%module(directors="1", threads="1", allprotected="1") grplinst

%include "exception.i"

%{
#include "CRPropa.h"
#include "grplinst.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/* include plugin parts to generate wrappers for */
%include "grplinst.h"





