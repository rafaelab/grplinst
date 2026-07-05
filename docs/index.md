---
layout: default
title: grplinst
---

# grplinst

`grplinst` is a CRPropa extension that models energy losses from plasma instabilities in electromagnetic cascades. It adds C++ and Python interfaces for beam models, intergalactic-medium profiles, and several literature-based cooling prescriptions for cascade electrons and positrons.

The code is intended for studies of blazar-induced gamma-ray cascades in the intergalactic medium, where plasma instabilities are treated as an effective cooling term. This is an approximate modelling layer rather than a first-principles particle-in-cell treatment, so model choice and parameter assumptions matter.

## Documentation Map

- [Getting Started](getting-started.html): build requirements, installation, and test commands.
- [Code Structure](code-structure.html): main abstractions, source layout, and public headers.
- [Models](models.html): implemented plasma-instability prescriptions and when each one is used.
- [Usage](usage.html): Python and C++ integration patterns for CRPropa simulations.
- [References](references.html): papers implemented by the code and related project links.

## What The Code Provides

- `PlasmaInstability*` modules that can be inserted into a CRPropa `ModuleList`.
- `Flow` implementations that describe the pair-beam density used by the cooling models.
- `MediumDensity` and `MediumTemperature` interfaces for homogeneous or user-defined environments.
- Helper functions such as `plasmaFrequency` and `maximumLinearGrowthFrequency`.
- SWIG bindings so the same components can be configured from Python.

## Quick Build

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## Entry Points

- Main umbrella header: `include/grplinst.h`
- Plasma-instability models: `include/grplinst/PlasmaInstability.h`
- Beam-density models: `include/grplinst/Flow.h`
- Medium profiles: `include/grplinst/Medium.h`
- Example Python pipeline: `examples/testPlugin.py`

## Citation

If you use `grplinst`, cite the main method paper listed on the [References](references.html) page.
