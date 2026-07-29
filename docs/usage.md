---
layout: default
title: Usage
---

# Usage

The recipe is always the same:
1. describe the ambient medium (`MediumDensity`, `MediumTemperature`),
2. describe the pair beam (`Flow`),
3. pick one instability prescription (`PlasmaInstability*`),
4. add it to a normal CRPropa electromagnetic-cascade `ModuleList`.

The instability module is just another CRPropa `Module`, so it composes with pair production, inverse-Compton scattering, redshift, observers, and break conditions
exactly as you would expect.

## Minimal Python example

```python
import sys
sys.path.append("build")   # so that `import grplinst` works

from crpropa import *
from grplinst import *

# 1. ambient medium
temperature = MediumTemperatureHomogeneous(1e4) # K
density = MediumDensityHomogeneous(1e-1) # m^-3

# 2. pair beam (argument is the source luminosity in W)
beam = FlowHomogeneous(1e38, Vector3d(0, 0, 0))

# 3. one prescription
plinst = PlasmaInstabilityBroderick2012(beam, density, temperature)

# 4. drop it into a ModuleList
sim = ModuleList()
sim.add(SimplePropagation(1e-3 * kpc, 10 * Mpc))
sim.add(Redshift())
sim.add(plinst)
```

## A realistic cascade pipeline

The snippet above only demonstrates wiring. In a real study you run the full
cascade so that the instability *competes* with inverse-Compton cooling. This is
the structure of
[`examples/testPlugin.py`](https://github.com/rafaelab/grplinst/blob/v2/examples/testPlugin.py):

```python
from crpropa import *
from grplinst import *

nEvents = 10000

# background photon fields
ebl = IRB_Saldana21()
cmb = CMB()

# source and environment
z = 0.0308  # e.g. Mrk 421
nIGM = 1e-1 # m^-3
L = 1e38 # W
T = 1e4  # K

# plasma-instability configuration
temperature = MediumTemperatureHomogeneous(T)
density = MediumDensityHomogeneous(nIGM)
beam = FlowHomogeneous(L, Vector3d(0, 0, 0))
plinst = PlasmaInstabilityBroderick2012(beam, density, temperature)

# cascade physics
electrons = photons = True
processes = [
    Redshift(),
    EMPairProduction(ebl, electrons),
    EMPairProduction(cmb, electrons),
    plinst,  # <-- the instability cooling
    EMInverseComptonScattering(ebl, photons),
    EMInverseComptonScattering(cmb, photons),
]

# injection: a power-law photon source at the blazar distance
source = Source()
source.add(SourcePowerLawSpectrum(1e9 * eV, 1e14 * eV, -1))
source.add(SourceParticleType(22))
source.add(SourcePosition(Vector3d(redshift2ComovingDistance(z), 0, 0)))
source.add(SourceRedshift1D())
source.add(SourceDirection(Vector3d(-1, 0, 0)))

# assemble
sim = ModuleList()
sim.add(SimplePropagation(1e-3 * kpc, 10 * Mpc))
for p in processes:
    sim.add(p)
# sim.add(MaximumTrajectoryLength(4000 * Mpc)) # only useful in 3D or non-relativistic simulations
sim.add(MinimumEnergy(1e9 * eV))

observer = Observer()
observer.add(Observer1D())
sim.add(observer)

sim.setShowProgress(True)
sim.run(source, nEvents, True)
```

**Where to place the module.** Add the instability in the same part of the `ModuleList` as the other interactions / energy losses. 
Both the instabilities module and inverse Compton scattering both act on the pairs each step. 

## Minimal C++ example

```cpp
#include <crpropa/ModuleList.h>
#include <crpropa/module/SimplePropagation.h>
#include <crpropa/module/Redshift.h>
#include <grplinst.h>

using namespace crpropa;
using namespace grplinst;

auto density = new MediumDensityHomogeneous(1e-1);
auto temperature = new MediumTemperatureHomogeneous(1e4);
auto beam = new FlowHomogeneous(1e38, Vector3d(0, 0, 0));
auto plinst = new PlasmaInstabilityVafin2018(beam, density, temperature);

ModuleList sim;
sim.add(new SimplePropagation(1e-3 * kpc, 10 * Mpc));
sim.add(new Redshift());
sim.add(plinst);
```

## Common patterns

### Homogeneous parameter scans

For broad comparisons, keep the environment homogeneous and vary the physical inputs and the model class:
```python
for L in [1e37, 1e38, 1e39]:
    for Model in [PlasmaInstabilityBroderick2012,
                  PlasmaInstabilitySironi2014,
                  PlasmaInstabilityVafin2018,
                  PlasmaInstabilityShalaby2020]:
        beam   = FlowHomogeneous(L)
        plinst = Model(beam, MediumDensityHomogeneous(1e-1), MediumTemperatureHomogeneous(1e4))
        # ... build and run a simulation, tag the output by (L, Model) ...
```

This is exactly the "run several models over the same grid and report the spread" workflow recommended in [Models](models.html#choosing-a-model).

### Luminosity-driven Miniati setup

`PlasmaInstabilityMiniati2013` constructs its own beam profile, so you pass a **luminosity** instead of a `Flow`:
```python
plinst = PlasmaInstabilityMiniati2013(1e38, density, temperature)
```
Internally this calls `createFlowMiniati2013(...)` to build a tabulated `FlowJet1D` (see [Models](models.html#plasmainstabilityminiati2013)).


### A tabulated beam profile

If you have a beam density as a function of distance from the source, use
`FlowJet1D` directly:
```python
import numpy as np
distances = np.array([1, 10, 100, 1000]) * Mpc # monotonic
densities = np.array([1e-15, 1e-17, 1e-19, 1e-21]) # m^-3, comoving

beam = FlowJet1D(distances, densities, 1e38, Vector3d(0, 0, 0), True) # last arg: log interpolation
```

Distances and densities are interpolated (logarithmically by default); outside the
tabulated range the nearest endpoint value is held, and the density is scaled by
`(1+z)³`.

### Custom profiles from Python (SWIG directors)

Because `Flow`, `MediumDensity`, and `MediumTemperature` are exposed as SWIG **directors**, you can subclass them in pure Python and override the virtual
methods -— no recompilation needed. 
The trivial examples below illustrate how to do that (although they are redundant with respect to existing classes).

```python
class ConstantDensity(MediumDensity):
    def getDensity(self, position, redshift = 0.):
        return 0.1 # m^-3

class ConstantTemperature(MediumTemperature):
    def getTemperature(self, position, redshift = 0.):
        return 1e4 # K

class ConstantBeam(Flow):
    def getDensity(self, energy, position, redshift = 0.):
        return 1e-16  # m^-3   (note the extra `energy` argument)

plinst = PlasmaInstabilityBroderick2012(ConstantBeam(), ConstantDensity(), ConstantTemperature())
```

These subclasses can then be passed anywhere the corresponding C++ type is expected. 
This is what is done in
[`test/testPython.py`](https://github.com/rafaelab/grplinst/blob/v2/test/testPython.py).

> **Signature reminder.** `Flow.getDensity` takes `(energy, position, redshift)`, whereas `MediumDensity.getDensity` takes `(position, redshift)`. 
> Getting these wrong is the most common integration mistake.

## Integration notes

- The instability acts **only on electrons and positrons**; photons are untouched.
- The `limit` parameter constrains the next step to a fraction of the local energy-loss length; this is important when comparing models with very short cooling times.
- The `efficiency` parameter (clamped to `[0, 1]`) scales the cooling strength; use `0` to switch the effect off for a baseline run.
- You can make the beam-density estimate consistent with the rest of the simulation by attaching the same `EMPairProduction` / `EMInverseComptonScattering` modules to the `Flow` with `addPairProduction` / `addInverseCompton` (see [API Reference](api-reference.html#flow)).

For the class-by-class reference and units, see [API Reference](api-reference.html).
