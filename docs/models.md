---
layout: default
title: Models
---

# Models

`grplinst` implements several effective plasma-instability prescriptions from the literature. All of them inherit from `PlasmaInstability`, so they are interchangeable at the CRPropa module level.

## Shared Inputs

All models work on cascade electrons and positrons and use some combination of:

- beam density from a `Flow` object,
- ambient density from a `MediumDensity` object,
- ambient temperature from a `MediumTemperature` object,
- particle energy and redshift from the current CRPropa candidate.

## Implemented Prescriptions

| Class | Reference | Depends on | Notes |
| --- | --- | --- | --- |
| `PlasmaInstabilityBroderick2012` | Broderick, Chang, Pfrommer (2012) | beam density, ambient density, energy | Piecewise behaviour separated by a critical beam density. |
| `PlasmaInstabilitySchlickeiser2012` | Schlickeiser, Ibscher, Supsar (2012) | beam density, ambient density, temperature, energy | Piecewise expression with explicit temperature dependence. |
| `PlasmaInstabilitySironi2014` | Sironi, Giannios (2014) | beam density, ambient density, energy | Piecewise prescription similar in structure to Broderick 2012. |
| `PlasmaInstabilityVafin2018` | Vafin, Pohl, Niemiec, Bret (2018) | beam density, ambient density, temperature, energy | Strong temperature dependence through an inverse scaling. |
| `PlasmaInstabilityBret2010TwoStream` | Bret, Gremillet, Dieckmann (2010) | beam density, ambient density, energy | Two-stream growth-based cooling time. |
| `PlasmaInstabilityBret2010Filamentation` | Bret, Gremillet, Dieckmann (2010) | beam density, energy | Filamentation estimate; current implementation does not use density or temperature in the final formula. |
| `PlasmaInstabilityShalaby2020` | Shalaby et al. (2020) | beam density, ambient density, energy | Effective scaling from the linear-instability analysis. |
| `PlasmaInstabilityMiniati2013` | Miniati, Elyiv (2013) | luminosity-derived flow, ambient density, energy | Builds its own `Flow` from source luminosity via `createFlowMiniati2013`. |

## Choosing A Model

The choice depends on what you want to compare.

- Use `Broderick2012`, `Sironi2014`, or `Schlickeiser2012` for direct comparisons among early effective prescriptions.
- Use `Vafin2018` or `Shalaby2020` when you want later literature fits included in the same CRPropa workflow.
- Use `Bret2010TwoStream` and `Bret2010Filamentation` when comparing specific instability channels.
- Use `Miniati2013` when you want the convenience of constructing the beam model from luminosity rather than passing a separate `Flow` object.

## Model Parameters Exposed By The Base Class

Beyond the literature-specific formula, every model inherits two practical control parameters:

- `efficiency`: multiplicative factor applied to the energy-loss rate. Values are clamped to the interval `[0, 1]`.
- `limit`: fraction used to restrict the next propagation step after a loss update.

These are useful when performing stability tests or when you want to suppress the effective loss rate without changing the model class itself.

## Helper Functions

Two free functions exposed by the library are useful when comparing analytical estimates:

- `plasmaFrequency(density, id=11)`
- `maximumLinearGrowthFrequency(beamDensity, mediumDensity, lorentzFactor, id=11)`

The unit tests in `test/testPlasmaInstability.cpp` cover their expected scaling behaviour.

For papers, code links, and the main project citation, see [References](references.html). For integration examples, continue with [Usage](usage.html).
