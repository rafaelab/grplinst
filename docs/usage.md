---
layout: default
title: Usage
---

# Usage

The typical workflow is to build a beam model, define the ambient medium, choose one plasma-instability prescription, and insert it into a normal CRPropa electromagnetic cascade setup.

## Minimal Python Example

```python
import sys
sys.path.append("../build")

from crpropa import *
from grplinst import *

temperature = MediumTemperatureHomogeneous(1e4)
density = MediumDensityHomogeneous(1e-1)
beam = FlowHomogeneous(1e38, Vector3d(0, 0, 0))

plinst = PlasmaInstabilityBroderick2012(beam, density, temperature)

sim = ModuleList()
sim.add(SimplePropagation(1e-3 * kpc, 10 * Mpc))
sim.add(Redshift())
sim.add(plinst)
```

This is the same composition pattern used in `examples/testPlugin.py`, where the instability module is combined with pair production, inverse Compton scattering, and an observer.

## Minimal C++ Example

```cpp
#include <crpropa/ModuleList.h>
#include <crpropa/module/SimplePropagation.h>
#include <crpropa/module/Redshift.h>
#include <grplinst.h>

using namespace crpropa;
using namespace grplinst;

auto density = new MediumDensityHomogeneous(1e-1);
auto temperature = new MediumTemperatureHomogeneous(1e4);
auto flow = new FlowHomogeneous(1e38, Vector3d(0, 0, 0));
auto plasma = new PlasmaInstabilityVafin2018(flow, density, temperature);

ModuleList sim;
sim.add(new SimplePropagation(1e-3 * kpc, 10 * Mpc));
sim.add(new Redshift());
sim.add(plasma);
```

## Common Use Cases

### Homogeneous Parameter Scans

For broad comparisons among literature models, use:

- `FlowHomogeneous`
- `MediumDensityHomogeneous`
- `MediumTemperatureHomogeneous`

This keeps the environmental assumptions fixed while you vary luminosity, density, temperature, redshift, and the chosen instability prescription.

### Luminosity-Driven Miniati Setup

If you want the `Miniati2013` prescription, instantiate `PlasmaInstabilityMiniati2013` directly from source luminosity plus medium properties:

```python
plinst = PlasmaInstabilityMiniati2013(1e38, density, temperature)
```

Internally this constructs the flow model with `createFlowMiniati2013(...)`.

### Custom Profiles From Python

The SWIG bindings support Python subclasses of `Flow`, `MediumDensity`, and `MediumTemperature`. That allows custom medium or beam models without recompiling the C++ core.

```python
class ConstantDensity(MediumDensity):
    def getDensity(self, position, redshift = 0.):
        return 0.1

class ConstantFlow(Flow):
    def getDensity(self, energy, position, redshift = 0.):
        return 1e-16
```

This path is exercised by the Python tests in `test/testPlugin.py`.

## Practical Integration Notes

- Add the plasma-instability module in the same part of the `ModuleList` where you want the effective cooling to compete with the rest of the cascade physics.
- The instability acts only on electrons and positrons; photons pass through unchanged.
- The `limit` parameter constrains the next propagation step after a loss update, which matters when you compare models with very short loss times.
- The `efficiency` parameter is a convenient way to scale the cooling strength for sensitivity tests.

## Example Pipelines In This Repository

- `examples/testPlugin.py`: end-to-end example using `PlasmaInstabilityBroderick2012`.
- `examples/_test.py`: scratch-style example with `PlasmaInstabilityMiniati2013`.
- `test/testPlugin.py`: examples of subclassing the public interfaces from Python.

For the available model classes and their literature background, see [Models](models.html).
