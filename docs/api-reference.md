---
layout: default
title: API Reference
---

# API Reference

This page documents the public classes and functions, their parameters, and the
units and redshift conventions they use. For the physics behind them see
[Physics Background](physics.html); for the cooling-time formulae see
[Models](models.html).

Everything lives in the `grplinst` namespace (C++) or the `grplinst` module
(Python). The umbrella header pulls in all components:

```cpp
#include <grplinst.h>   // Flow, Medium, Geometry, PlasmaInstability
```

## Architecture at a glance

`grplinst` is deliberately small. A run is composed from four ingredients:

```
  MediumDensity ──┐
  MediumTemperature ─┤
  Flow (beam density) ┼──►  PlasmaInstability*  ──►  CRPropa ModuleList
                      │        (a CRPropa Module)
        energy, z  ───┘
```

- **`Flow`** answers "how dense is the pair beam here?"
- **`MediumDensity` / `MediumTemperature`** describe the ambient plasma.
- **`PlasmaInstability`** (abstract) is the CRPropa module. Concrete subclasses
  differ *only* in their `energyLossTime` implementation; the base class handles
  everything else (particle selection, redshift bookkeeping, step limiting).

All classes derive from CRPropa's reference-counted `Referenced`/`Module` types
and are held through `crpropa::ref_ptr`, so memory is managed automatically.

Class hierarchy:

```
crpropa::Module
└── PlasmaInstability                 (abstract)
    ├── PlasmaInstabilityBroderick2012
    ├── PlasmaInstabilitySironi2014
    ├── PlasmaInstabilitySchlickeiser2012
    ├── PlasmaInstabilityVafin2018
    ├── PlasmaInstabilityBret2010TwoStream
    ├── PlasmaInstabilityBret2010Filamentation
    ├── PlasmaInstabilityShalaby2020
    └── PlasmaInstabilityMiniati2013

crpropa::Referenced
├── Flow                              (abstract)
│   ├── FlowHomogeneous
│   └── FlowJet1D
├── MediumDensity                     (abstract)
│   └── MediumDensityHomogeneous
├── MediumTemperature                 (abstract)
│   └── MediumTemperatureHomogeneous
└── EmissionGeometry                  (abstract)
    └── Cone
```

## `PlasmaInstability` (base class) {#plasmainstability-base-class}

Header: `include/grplinst/PlasmaInstability.h`. The abstract CRPropa module that
applies the effective cooling.

**Constructor**

```cpp
PlasmaInstability(ref_ptr<Flow> flow,
                  ref_ptr<MediumDensity> density,
                  ref_ptr<MediumTemperature> temperature,
                  double efficiency = 1.0,
                  double limit = 0.1);
```

| Parameter | Meaning |
| --- | --- |
| `flow` | beam-density model ([`Flow`](#flow)) |
| `density` | ambient density model ([`MediumDensity`](#mediumdensity)) |
| `temperature` | ambient temperature model ([`MediumTemperature`](#mediumtemperature)) |
| `efficiency` | scales `dE/dx`; clamped to `[0, 1]` (default `1`) |
| `limit` | next step ≤ this fraction of the loss length (default `0.1`) |

**Configuration methods**

`setFlowProperties`, `setMediumDensity`, `setMediumTemperature`,
`setEfficiencyFactor`, `setLimit`, and the matching getters
`getFlowProperties`, `getMediumDensity`, `getMediumTemperature`,
`getEfficiencyFactor`, `getLimit`.

`setEfficiencyFactor` enforces the `[0, 1]` clamp and logs a warning if the value
is out of range.

**What `process(candidate)` does**

For each candidate the base class:

1. returns immediately unless the particle is an electron or positron
   (`|id| == 11`);
2. evaluates the local-frame energy `E = E_obs·(1+z)` and de-redshifted step
   `dx = Δx/(1+z)`;
3. asks the subclass for the cooling time `τ = energyLossTime(candidate)`;
4. forms the energy loss per length
   `dE/dx = efficiency · E / (c·τ)` via `computeEnergyLossPerLength`
   (returns `0` if `τ ≤ 0`, which safely disables the loss);
5. updates the energy to `max(0, E − dE/dx·dx)`, stored back as
   `E_new/(1+z)`;
6. limits the next step to `limit · E / (dE/dx)`.

**Virtual interface (implemented by each model)**

```cpp
virtual double energyLossTime(const crpropa::Candidate& candidate) const = 0;
```

Returns the characteristic cooling time `τ` **in seconds**. This is the only
method a new prescription must provide. The concrete subclasses are listed in
[Models](models.html).

## `Flow` {#flow}

Header: `include/grplinst/Flow.h`. Abstract base describing the pair-beam density.

**Key virtual method**

```cpp
virtual double getDensity(double energy,
                          const crpropa::Vector3d& position,
                          double redshift = 0) const = 0;
```

Returns the beam number density (m⁻³) for a pair of the given `energy` at
`position` and `redshift`. *(Note the leading `energy` argument — this
distinguishes it from `MediumDensity::getDensity`.)*

**Shared configuration**

- `setOrigin` / `getOrigin` — the flow origin (source position).
- `setLuminosity` / `getLuminosity` — source luminosity (W). Default `1`.
- `setPairProduction` / `addPairProduction` and
  `setInverseCompton` / `addInverseCompton` — optionally attach CRPropa
  `EMPairProduction` / `EMInverseComptonScattering` modules so that the beam
  density is computed from their **actual rates** instead of the built-in
  analytic approximations.

### `FlowHomogeneous`

Constant-luminosity beam whose density follows the Broderick et al. (2012)
estimate (see [Physics §4.1](physics.html#41-the-pair-beam-density-flow)).

```cpp
FlowHomogeneous();                                       // luminosity = 1
FlowHomogeneous(double luminosity,
                crpropa::Vector3d origin = {0, 0, 0});
```

`getDensity` returns `L / (2π·λγγ³·Γ_IC) / E`, using analytic `λγγ` and `Γ_IC`
unless pair-production / inverse-Compton modules have been attached, in which case
their rates are used. This is an **upper limit** on the true beam density.

### `FlowJet1D`

A one-dimensional beam profile `n_b(r)` tabulated along the line of sight.

```cpp
FlowJet1D(const std::vector<double>& distances,
          const std::vector<double>& beamDensity,
          double luminosity = 1,
          crpropa::Vector3d centre = {0, 0, 0},
          bool interpolateLog = true);
```

- `distances` must be monotonic; `beamDensity` must have the same length
  (a mismatch throws `std::invalid_argument`).
- With `interpolateLog = true`, distances are interpolated in `log10`.
- Outside the tabulated range the nearest endpoint value is held.
- The returned density is scaled by `(1+z)³`.

Accessors: `setDistanceProfile`, `setDensityProfile`, `setInterpolateLog`,
`getDistanceProfile`, `getDensityProfile`.

### `createFlowMiniati2013`

```cpp
ref_ptr<Flow> createFlowMiniati2013(double luminosity,
                                    crpropa::Vector3d centre = {0, 0, 0},
                                    bool logDistance = true);
```

Builds a `FlowJet1D` from the 16-point beam-density profile of Miniati & Elyiv
(2013), scaled linearly with `luminosity / 10³⁸ W`. Used internally by
[`PlasmaInstabilityMiniati2013`](models.html#plasmainstabilityminiati2013).

## `MediumDensity` {#mediumdensity}

Header: `include/grplinst/Medium.h`. Abstract base for the ambient plasma density.

```cpp
virtual double getDensity(const crpropa::Vector3d& position,
                          const double& redshift = 0.) const = 0;
```

### `MediumDensityHomogeneous`

```cpp
MediumDensityHomogeneous(double density);   // comoving density in m^-3
```

`getDensity` returns `n₀·(1+z)³`. Accessors: `setDensityValue`,
`getDensityValue`.

## `MediumTemperature` {#mediumtemperature}

Header: `include/grplinst/Medium.h`. Abstract base for the ambient temperature.

```cpp
virtual double getTemperature(const crpropa::Vector3d& position,
                              const double& redshift = 0.) const = 0;

// convenience: thermal velocity sqrt(k_B T / m) for particle `id`
double getVelocity(int id,
                   const crpropa::Vector3d& position,
                   const double& redshift = 0) const;
```

### `MediumTemperatureHomogeneous`

```cpp
MediumTemperatureHomogeneous(double temperature);   // K
```

`getTemperature` returns `T₀·(1+z)`. Accessors: `setTemperatureValue`,
`getTemperatureValue`.

## Helper functions

```cpp
double plasmaFrequency(double density, int id = 11);

double maximumLinearGrowthFrequency(double beamDensity,
                                    double mediumDensity,
                                    double lorentzFactor,
                                    int id = 11);
```

- **`plasmaFrequency(n, id)`** — `ω_p = sqrt(n·q²/(m·ε₀))` in Hz. For electrons
  (`id = 11`) `q = e` and `m = m_e`; nuclei use their charge number and nuclear
  mass.
- **`maximumLinearGrowthFrequency(n_b, n, γ, id)`** — the maximum linear growth
  rate `ω_p(n)·(n_b/n)·(1/γ)` in Hz. Neglects magnetic fields and assumes an
  angular spread `Δθ = <1/γ>`. This is the rate inverted by
  [`PlasmaInstabilityMiniati2013`](models.html#plasmainstabilityminiati2013).

## Geometry (infrastructure)

Header: `include/grplinst/Geometry.h`. `EmissionGeometry` (abstract) and its
concrete `Cone` describe an emission region and can compute an area and a volume.
These classes are **provided for future emission-geometry calculations** and are
not used by the instability models; you can ignore them for standard cascade runs.

## Constants and units

`grplinst` uses SI units throughout, consistent with CRPropa. Where CRPropa unit
symbols (`TeV`, `Mpc`, `kpc`, `eV`, …) are available they are used.

| Quantity | Unit |
| --- | --- |
| Energy | J internally; CRPropa exposes `eV`, `TeV`, … multipliers |
| Number density (`n`, `n_b`) | m⁻³ |
| Temperature | K |
| Luminosity | W |
| Length / position | m (`Mpc`, `kpc` multipliers) |
| Cooling time `τ`, frequencies | s, Hz |

Two module-internal constants worth noting (`include/grplinst/Common.h`):
`u_CMB = 4.178×10⁻¹⁴ J m⁻³` (CMB energy density) and
`mec2 = m_e c²` (electron rest energy).

## Python bindings

The SWIG interface (`python/grplinst.i`) exposes every public class and function
under the `grplinst` Python module, with a few conveniences:

- `Flow`, `MediumDensity`, `MediumTemperature`, `PlasmaInstability`, and
  `EmissionGeometry` are **directors**, so you can subclass them in Python and
  override their virtual methods (see
  [Usage](usage.html#custom-profiles-from-python-swig-directors)).
- `std::vector<double>` arguments and return values interoperate with Python
  lists and NumPy arrays.
- `ref_ptr` conversions are implicit, so you pass ordinary Python objects wherever
  a `ref_ptr<...>` is expected.
