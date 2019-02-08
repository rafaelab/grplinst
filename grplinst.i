%module(directors="1", threads="1", allprotected="1") GRPlInst

%include "exception.i"

%{
#include "CRPropa.h"
#include "GRPlInst.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/* include plugin parts to generate wrappers for */
%include "GRPlInst.h"





