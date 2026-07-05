---
layout: default
title: Getting Started
---

# Getting Started

This page covers the practical setup needed to build and test `grplinst` against CRPropa.

## Requirements

- CRPropa 3 installation available through `find_package(CRPropa)` or built through the `USE_OWN_CRPROPA` option.
- A C++17 compiler.
- CMake 3.14 or newer.
- Python 3, SWIG 4+, and NumPy if you want the Python bindings.

## Build Options

The top-level `CMakeLists.txt` exposes a few options that are useful when preparing local builds or CI jobs:

- `ENABLE_PYTHON=On|Off`: build the SWIG Python module.
- `ENABLE_TESTING=On|Off`: compile the unit tests.
- `USE_OWN_CRPROPA=On|Off`: fetch and build CRPropa instead of using an installed copy.

## Typical Local Build

```bash
cmake -S . -B build
cmake --build build
```

If CRPropa is not installed in a default location, point CMake to it with the usual variables such as `CMAKE_PREFIX_PATH` or the package-specific paths used by your local installation.

## Python Bindings

With Python support enabled, the build creates the `grplinst` Python package inside the build tree. The examples in this repository typically extend `sys.path` with the build directory before importing the module.

Example:

```python
import sys
sys.path.append("../build")

from crpropa import *
from grplinst import *
```

## Running Tests

The repository contains both C++ and Python tests.

```bash
ctest --test-dir build --output-on-failure
```

Useful test files to inspect while integrating the library:

- `test/testPlasmaInstability.cpp`: model behaviour and helper-function tests.
- `test/testMedium.cpp`: medium-profile coverage.
- `test/testFlow.cpp`: beam-model coverage.
- `test/testPlugin.py`: Python binding and SWIG-director coverage.
