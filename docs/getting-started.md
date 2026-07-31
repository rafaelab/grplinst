---
layout: default
title: Getting Started
---

# Getting Started

This page walks through building `grplinst` against CRPropa and running the test suite. 
`grplinst` is a CRPropa *plugin*: it compiles to a shared library and an optional Python module that are loaded alongside CRPropa itself.

## Prerequisites
- **CRPropa** (≥ 3.3): Either an existing installation discoverable through `find_package(CRPropa)`, or built on the fly with the `USE_OWN_CRPROPA` option (see below).
- **C++17**: The compiler, compatible with C++17 standards. Tested for some versions of GCC and Clang. 
- **Fortran compiler** (`gfortran`): Needed by CRPropa's dependencies when it is built from source. 
- **CMake** ≥ 3.14
- **SWIG 4+**, **Python 3**, **NumPy**: Only required for the Python bindings (`ENABLE_PYTHON`, on by default). 

On Debian/Ubuntu the toolchain is:
```bash
sudo apt-get update
sudo apt-get install --yes cmake g++ gfortran make swig
python3 -m pip install --upgrade numpy
```

On macOS with Homebrew:
```bash
brew install gcc swig cmake
python3 -m pip install --upgrade numpy
```

## Getting the code

```bash
git clone https://github.com/rafaelab/grplinst.git
cd grplinst
```

## Build options

The top-level `CMakeLists.txt` exposes the following switches (defaults in brackets):
- `ENABLE_PYTHON` (default: `On`): Build the SWIG Python module.
- `ENABLE_TESTING` (default: `On`): Build the C++ (GoogleTest) and Python unit tests.
- `USE_OWN_CRPROPA` (default: `Off`): Fetch and build CRPropa from source via CMake `FetchContent` instead of using an installed copy.
- `ENABLE_COVERAGE` (default: `Off`):  Instrument the build for coverage analysis. (Not yet implemented).
- `ENABLE_PYBIND11_STUBGEN` (default: `On`): Generate Python stub files for the bindings.
Pass any option with `-D<NAME>=On|Off`.

## Building against an installed CRPropa (recommended)

If you already have CRPropa installed, point CMake at it and build:
```bash
cmake -S . -B build
cmake --build build --parallel 4
```

If CRPropa is not in a default location, help CMake find it, for example:
```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="$HOME/.local"
```

The bundled `cmake/FindCRPropa.cmake` module locates the CRPropa headers, library, and SWIG interface files.

## Building a self-contained copy of CRPropa

For a clean-room build (this is what continuous integration uses), let the project download and build CRPropa for you:
```bash
cmake -S . -B build -DUSE_OWN_CRPROPA=On
cmake --build build --parallel 4
```

This is convenient but slower, because it compiles CRPropa in full the first time.

## Using the Python module

With `ENABLE_PYTHON=On`, the build produces the `grplinst` Python package inside
the build tree. Add that directory to your `PYTHONPATH`, or extend `sys.path`
from within your script (as the bundled examples do):
```python
import sys
sys.path.append("build")   # path to the build directory

from crpropa import *
from grplinst import *
```
The same rule applies to CRPropa: its Python module must be importable. If you
built CRPropa yourself, make sure its build directory is on the `PYTHONPATH` too.

## Running the tests

The repository ships both C++ and Python tests, wired into CTest:
```bash
ctest --test-dir build --output-on-failure
```

The individual test files are worth reading as compact, executable documentation:

| File | Coverage |
| --- | --- |
| [`test/testPlasmaInstability.cpp`](https://github.com/rafaelab/grplinst/blob/v2/test/testPlasmaInstability.cpp) | Every model's cooling time, the base-class energy-loss machinery, and the helper functions. |
| [`test/testFlow.cpp`](https://github.com/rafaelab/grplinst/blob/v2/test/testFlow.cpp) | `FlowHomogeneous`, `FlowJet1D`, and the Miniati (2013) beam profile. |
| [`test/testMedium.cpp`](https://github.com/rafaelab/grplinst/blob/v2/test/testMedium.cpp) | Homogeneous density and temperature profiles. |
| [`test/testPython.py`](https://github.com/rafaelab/grplinst/blob/v2/test/testPython.py) | Import checks, `std::vector` ↔ list/NumPy conversions, and Python sub-classing through SWIG directors. |

If the C++ tests build but do not run, verify that CRPropa was itself built with its bundled GoogleTest available (CRPropa ≥ 3.3 exposes it), since the test target reuses it.

Once the tests pass you are ready to build a simulation. 
Continue with [Usage](usage.html), or read the [Physics Background](physics.html) first to understand what the module actually computes.
