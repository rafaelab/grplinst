---
layout: default
title: Code Structure
---

# Code Structure

`grplinst` is intentionally small. Most of the public surface is concentrated in three abstractions that are combined by a `PlasmaInstability` module during CRPropa propagation.

## Main Abstractions

### `PlasmaInstability`

The abstract base class in `include/grplinst/PlasmaInstability.h` defines the CRPropa module interface. Its `process` method:

1. selects only electrons and positrons,
2. queries the model-specific energy-loss time,
3. converts that loss time into `dE/dx`,
4. updates the candidate energy, and
5. limits the next propagation step through the configurable `limit` parameter.

Each literature model only needs to implement `energyLossTime(candidate)`.

### `Flow`

The `Flow` hierarchy in `include/grplinst/Flow.h` supplies the beam density seen by the instability model.

- `FlowHomogeneous` is the simplest choice and is the default starting point for homogeneous scans.
- `FlowJet1D` allows a user-supplied one-dimensional profile.
- `createFlowMiniati2013(...)` builds the beam setup used by the `Miniati2013` prescription.

### `MediumDensity` and `MediumTemperature`

The medium-profile interfaces in `include/grplinst/Medium.h` describe the ambient intergalactic medium.

- `MediumDensityHomogeneous` returns a constant density.
- `MediumTemperatureHomogeneous` returns a constant temperature.
- Python subclasses can override the virtual interfaces through SWIG directors, which is exercised in `test/testPlugin.py`.

## Source Layout

- `src/PlasmaInstability.cc`: implemented loss-time formulae and helper functions.
- `src/Flow.cc`: beam-density models.
- `src/Medium.cc`: homogeneous medium models and thermal helper calculations.
- `src/Geometry.cc`: geometry helpers used by the project.
- `include/grplinst.h`: umbrella include for consumers.
- `python/`: SWIG interface files and Python packaging stubs.
- `examples/`: small end-to-end examples.
- `test/`: C++ and Python regression tests.

## Public Headers

Most users only need one of these headers:

- `include/grplinst.h`: include everything.
- `include/grplinst/PlasmaInstability.h`: only the instability models.
- `include/grplinst/Flow.h`: only the beam-density models.
- `include/grplinst/Medium.h`: only medium-density and temperature profiles.

## How Components Fit Together

The runtime composition is:

1. create a beam-density model (`Flow...`),
2. create ambient medium models (`MediumDensity...`, `MediumTemperature...`),
3. instantiate a literature model (`PlasmaInstability...`),
4. add that module to a CRPropa `ModuleList`,
5. run the simulation with standard CRPropa interaction and observer modules.

This composition pattern is shown in more detail on the [Usage](usage.html) page.

## Practical Notes

- The implementation is designed as an effective cooling prescription, not a microscopic kinetic simulation.
- Different models depend on different subsets of the inputs. Some use beam density and ambient density only, while others also depend on temperature.
- The `efficiency` parameter can rescale the effective losses between 0 and 1, and values outside that range are clamped.

For the implemented literature prescriptions and their input dependencies, see [Models](models.html).
