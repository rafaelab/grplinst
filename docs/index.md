---
layout: default
title: grplinst
---

# grplinst

**`grplinst`** ("group–plasma instabilities") is an extension for the
[CRPropa 3](https://github.com/CRPropa/CRPropa3) propagation code. It models the
energy losses that **plasma (beam) instabilities** are expected to inflict on the
electron–positron pairs produced in **blazar-induced electromagnetic cascades**
as they travel through the intergalactic medium (IGM).

The module supplies:

- a family of `PlasmaInstability*` modules that plug straight into a CRPropa
  `ModuleList`, each implementing a different literature prescription for the
  instability cooling time;
- `Flow` classes that describe the **pair-beam density** the cascade drives into
  the IGM;
- `MediumDensity` and `MediumTemperature` classes that describe the ambient
  plasma;
- helper functions such as `plasmaFrequency` and
  `maximumLinearGrowthFrequency`;
- full **Python bindings** through SWIG, so every component can be configured,
  sub-classed, and combined from Python.

> **A word of caution.** Plasma instabilities in a dilute relativistic pair beam
> are a genuinely hard kinetic problem. `grplinst` does **not** solve it from
> first principles; it wraps a set of *effective* cooling prescriptions taken
> from the literature and lets you compare them within the same cascade
> simulation. A fully self-consistent treatment requires particle-in-cell (PIC)
> methods. Please read the [Physics Background](physics.html) page before drawing
> quantitative conclusions.

## Where to go next

| Page | What it covers |
| --- | --- |
| [Getting Started](getting-started.html) | Prerequisites, building against CRPropa, running the tests. |
| [Physics Background](physics.html) | Blazar cascades, pair beams, why instabilities matter, the effective-cooling approximation, and its limits. |
| [Models](models.html) | Every implemented prescription, its cooling-time formula, regimes, and inputs. |
| [Usage](usage.html) | Python and C++ integration patterns, parameter scans, and custom profiles. |
| [API Reference](api-reference.html) | Class hierarchy, methods, parameters, units, and redshift conventions. |
| [References](references.html) | The method paper (please cite it) and the source of each model. |

## A thirty-second example

```python
from crpropa import *
from grplinst import *

# ambient intergalactic medium and pair beam
temperature = MediumTemperatureHomogeneous(1e4)     # K
density     = MediumDensityHomogeneous(1e-1)         # m^-3
beam        = FlowHomogeneous(1e38, Vector3d(0, 0, 0))  # source luminosity in W

# one instability prescription, ready to drop into a ModuleList
plinst = PlasmaInstabilityBroderick2012(beam, density, temperature)

sim = ModuleList()
sim.add(SimplePropagation(1e-3 * kpc, 10 * Mpc))
sim.add(Redshift())
sim.add(plinst)
```

A complete, physically meaningful pipeline (with pair production and
inverse-Compton scattering) is given on the [Usage](usage.html) page and in
[`examples/testPlugin.py`](https://github.com/rafaelab/grplinst/blob/v2/examples/testPlugin.py).

## Citing grplinst

If `grplinst` contributes to your work, please cite the method paper:

> R. Alves Batista, A. Saveliev, E. M. de Gouveia Dal Pino,
> *The impact of plasma instabilities on the spectra of TeV blazars*,
> MNRAS **489** (2019) 3836.
> [doi:10.1093/mnras/stz2389](https://doi.org/10.1093/mnras/stz2389) ·
> [arXiv:1904.13345](https://arxiv.org/abs/1904.13345)

See [References](references.html) for the full bibliography, including the origin
of each instability model.
