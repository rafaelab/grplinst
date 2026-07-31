---
layout: default
title: Models
---

# Models

`grplinst` implements several effective plasma-instability prescriptions from the literature. 
Each one is a subclass of [`PlasmaInstability`](api-reference.html#plasmainstability-base-class) and provides a single function —- the **cooling time** ($$\tau$$) —- which the base class turns into an energy loss.
For details, see [Physics Background §3](physics.html#3-how-grplinst-models-it-effective-cooling).
Because they share the same interface, most of the models are **interchangeable**: swap one class for another and rerun.

## Notation



Every formula below is written in terms of the following *scaled variables*, so the numbers match the source code directly. $$\tau$$ is in *seconds*.

| Symbol | Definition |
| --- | --- | --- |
| $$\tilde{E}$$ | particle energy $$E / \mathrm{TeV}$$ (local frame, $$E = E_\mathrm{obs}(1+z)$$) |
| $$\tilde{n}_b$$ | beam density $$n_b / 10^{-16}\,\mathrm{m^{-3}}$$ (from the [`Flow`](api-reference.html#flow)) |
| $$\tilde{n}$$ | ambient density $$n / 0.1\,\mathrm{m^{-3}}$$ (from the [`MediumDensity`](api-reference.html#mediumdensity)) |
| $$\tilde{T}$$ | ambient temperature $$T / 10^{4}\,\mathrm{K}$$ (from the [`MediumTemperature`](api-reference.html#mediumtemperature)) |

A model uses only the inputs that appear in its formula; the others may be omitted (pass any placeholder object). 
Several models are **piecewise**, switching between a weak/linear branch and a saturated branch at a critical beam density $$n_\mathrm{crit}$$.

## Implemented prescriptions

* `PlasmaInstabilityBret2010TwoStream`
* `PlasmaInstabilityBret2010Filamentation`
* `PlasmaInstabilityBroderick2012`
* `PlasmaInstabilityMiniati2013`
* `PlasmaInstabilitySironi2014`
* `PlasmaInstabilitySchlickeiser2012`
* `PlasmaInstabilityVafin2018`
* `PlasmaInstabilityShalaby2020`

These models generally provide a cooling term ($$\mathrm{d}E / \mathrm{d}x$$).
For details, see the `grplinst` paer.



## Control parameters 

Every model inherits two practical knobs from the base class, both optional constructor arguments:
- *`efficiency`* ($$\eta$$, defaults to 1): multiplies the energy-loss rate; clamped to $$[0, 1]$$ — negative values become `0`, values above `1` become `1` (a warning is logged). Set to `0` to disable the instability, or to a fraction to test sensitivity.
- *`limit`* (defaults to 0.1): the next propagation step is restricted to this fraction of the local energy-loss length, so short cooling times are resolved.
```python
plinst = PlasmaInstabilityVafin2018(beam, density, temperature, 0.3, 0.05) # 30% efficiency, tighter step control
```


---

For the papers behind each model, see [References](references.html). \
For CRPropa integration, continue with [Usage](usage.html).
