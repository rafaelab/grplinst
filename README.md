# grplinst

[![Documentation](https://github.com/rafaelab/grplinst/actions/workflows/docs.yml/badge.svg)](https://rafaelab.github.io/grplinst/)
[![Tests](https://github.com/rafaelab/grplinst/actions/workflows/tests.yml/badge.svg)](https://github.com/rafaelab/grplinst/actions/workflows/tests.yml)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fmnras%2Fstz2389-blue)](https://doi.org/10.1093/mnras/stz2389)

Module for the CRPropa code to calculate the effect of plasma instabilities on the development of electromagnetic cascades.
It is described in:

R. Alves Batista, A. Saveliev, E. M. de Gouveia Dal Pino, *The impact of plasma instabilities on the spectra of TeV blazars*, MNRAS **489** (2019) 3836. [doi:10.1093/mnras/stz2389](https://doi.org/10.1093/mnras/stz2389) · [arXiv:1904.13345](https://arxiv.org/abs/1904.13345)

If you make use of grplinst, please consider citing the above paper.

## 📖 Documentation

Full, pedagogical documentation is available at **[rafaelab.github.io/grplinst](https://rafaelab.github.io/grplinst/)**:

- [Getting Started](https://rafaelab.github.io/grplinst/getting-started.html) — building and testing.
- [Physics Background](https://rafaelab.github.io/grplinst/physics.html) — blazar cascades, pair beams, and the effective-cooling approximation.
- [Models](https://rafaelab.github.io/grplinst/models.html) — the implemented prescriptions and their cooling-time formulae.
- [Usage](https://rafaelab.github.io/grplinst/usage.html) — Python and C++ integration.
- [API Reference](https://rafaelab.github.io/grplinst/api-reference.html) — classes, methods, and units.

## Science

This code was designed to study the effects of plasma instabilities on the development of blazar-induced gamma-ray cascades in the intergalactic medium.

It provides *approximate* descriptions for some plasma-instability models found in the literature.
They are modelled as a cooling term for electrons/positrons. Note that this is a very rough approximation and a detailed calculation would require particle-in-cell (PIC) simulations.

## Installation

To install *grplinst* you need CRPropa 3 installed. Go to <https://github.com/CRPropa/CRPropa3/> and follow the instructions.

Then install this module:

1. Download the latest version of the code:
   ```
   git clone https://github.com/rafaelab/grplinst.git
   cd grplinst
   ```

2. Configure and build with CMake:
   ```
   cmake -S . -B build
   cmake --build build --parallel 4
   ```
   If you do not have CRPropa installed, you can let CMake fetch and build it for you by adding `-DUSE_OWN_CRPROPA=On` to the first command.

3. If the code compiled, you are ready to go. Make sure the build directory (where the `grplinst` Python module is created) is on your `PYTHONPATH`, or add it to `sys.path` in your script.

See the [Getting Started](https://rafaelab.github.io/grplinst/getting-started.html) page for build options, dependencies, and how to run the tests.

## Disclaimer

This program is provided 'as is', without warranties of any kind.
Please use your discernment to interpret the results obtained with it.

## Acknowledgements

We thank Michael Kachelrieß and Jacob Benestad for providing feedback that led us to finding a few bugs.
