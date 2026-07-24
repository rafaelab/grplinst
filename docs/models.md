---
layout: default
title: Models
---

# Models

`grplinst` implements several effective plasma-instability prescriptions from the
literature. Each one is a subclass of
[`PlasmaInstability`](api-reference.html#plasmainstability-base-class) and
provides a single function — the **cooling time** `τ` — which the base class turns
into an energy loss (see [Physics Background §3](physics.html#3-how-grplinst-models-it-effective-cooling)).
Because they share the same interface, the models are **interchangeable**: swap
one class for another and rerun.

## Notation

Every formula below is written in terms of the following **scaled variables**, so
the numbers match the source code directly. `τ` is in **seconds**.

| Symbol | Definition | Fiducial value |
| --- | --- | --- |
| `Ẽ` | particle energy `E / TeV` (in the local frame, `E = E_obs·(1+z)`) | 1 |
| `ñ_b` | beam density `n_b / 10⁻¹⁶ m⁻³` (from the [`Flow`](api-reference.html#flow)) | 1 |
| `ñ` | ambient density `n / 0.1 m⁻³` (from the [`MediumDensity`](api-reference.html#mediumdensity)) | 1 |
| `T̃` | ambient temperature `T / 10⁴ K` (from the [`MediumTemperature`](api-reference.html#mediumtemperature)) | 1 |

A model uses only the inputs that appear in its formula; the others may be omitted
(pass any placeholder object). Several models are **piecewise**, switching between
a weak/linear branch and a saturated branch at a critical beam density `n_crit`.

## Implemented prescriptions

### `PlasmaInstabilityBroderick2012`

Broderick, Chang & Pfrommer (2012). Piecewise in the beam density:

```
n_crit = 1.6e-13 · ñ / Ẽ²                             [m⁻³]

n_b < n_crit :   τ = 7e7 · sqrt(ñ) / (Ẽ · ñ_b)
n_b ≥ n_crit :   τ = 5e5 · Ẽ^(1/3) / ( ñ_b^(1/3) · ñ^(1/6) )
```

Uses beam density, ambient density, and energy. The weak branch cools roughly as
`1/(E·n_b)`; above `n_crit` the dependence flattens to cube roots.

### `PlasmaInstabilitySironi2014`

Sironi & Giannios (2014). Same structure as Broderick (2012), different constants:

```
n_crit = 8e-14 · ñ / Ẽ²                               [m⁻³]

n_b < n_crit :   τ = 1.4e7 · sqrt(ñ) / (Ẽ · ñ_b)
n_b ≥ n_crit :   τ = 9.6e5 · Ẽ^(1/3) / ( ñ_b^(1/3) · ñ^(1/6) )
```

Uses beam density, ambient density, and energy.

### `PlasmaInstabilitySchlickeiser2012`

Schlickeiser, Ibscher & Supsar (2012). Piecewise, with an explicit temperature
dependence in both the threshold and the cooling time:

```
n_crit = 2.5e-19 · ñ · T̃² / Ẽ                        [m⁻³]

n_b < n_crit :   τ = 5e14 · Ẽ^(5/3) · ñ_b^(1/3) / ( ñ^(5/6) · T̃² )
n_b ≥ n_crit :   τ = 8e6 · Ẽ^(1/3) / ( ñ_b^(1/3) · ñ^(1/6) )
                       · [ 1 + 1.25·ln(T̃) − 0.25·ln(ñ) ]
```

Uses beam density, ambient density, temperature, and energy.

### `PlasmaInstabilityVafin2018`

Vafin, Pohl, Niemiec & Bret (2018). A single regime with a strong inverse
temperature dependence:

```
τ = 1.9e11 · Ẽ^(4/3) · ñ^(1/3) / ( ñ_b^(1/3) · T̃ )
```

Uses beam density, ambient density, temperature, and energy. Because `τ ∝ 1/T`, a
hotter medium cools the beam faster in this prescription.

### `PlasmaInstabilityBret2010TwoStream`

Bret, Gremillet & Dieckmann (2010), two-stream channel:

```
τ = 1.6e10 · Ẽ / ( ñ_b^(1/3) · ñ^(1/6) )
```

Uses beam density, ambient density, and energy.

### `PlasmaInstabilityBret2010Filamentation`

Bret, Gremillet & Dieckmann (2010), filamentation channel:

```
τ = 2.5e9 · sqrt(Ẽ) / sqrt(ñ_b)
```

Uses beam density and energy only. **Note:** although a `MediumDensity` and
`MediumTemperature` must be supplied to construct the module, the final
filamentation expression does not use them.

### `PlasmaInstabilityShalaby2020`

Shalaby et al. (2020), effective scaling from the linear-instability analysis:

```
τ = 3.2e12 · Ẽ^(6/5) · ñ_b^(-2/5) · ñ^(-1/10)
```

Uses beam density, ambient density, and energy.

### `PlasmaInstabilityMiniati2013`

Miniati & Elyiv (2013). Rather than a fitted cooling time, this model uses the
**inverse of the maximum linear growth rate**:

```
τ = 1 / ω_growth ,   ω_growth = ω_p(n) · (n_b / n) · (1/γ)
```

where `γ = E / (m_e c²)` and `ω_p(n) = sqrt(n e² / (m_e ε₀))` is the plasma
frequency of the ambient medium (see [Helper functions](#helper-functions)).

This model is special in that it **builds its own beam density**: the constructor

```python
PlasmaInstabilityMiniati2013(luminosity, density, temperature)
```

calls `createFlowMiniati2013(luminosity)` internally, which sets up a
16-point tabulated `FlowJet1D` profile (a numerical `n_b(r)` from the paper)
scaled linearly with `luminosity / 10³⁸ W`. You therefore pass a **luminosity**,
not a `Flow`. See [Usage](usage.html#luminosity-driven-miniati-setup).

## Summary table

| Class | Reference | Inputs used | Regimes |
| --- | --- | --- | --- |
| `PlasmaInstabilityBroderick2012` | Broderick+ 2012 | `n_b`, `n`, `E` | piecewise |
| `PlasmaInstabilitySironi2014` | Sironi & Giannios 2014 | `n_b`, `n`, `E` | piecewise |
| `PlasmaInstabilitySchlickeiser2012` | Schlickeiser+ 2012 | `n_b`, `n`, `T`, `E` | piecewise |
| `PlasmaInstabilityVafin2018` | Vafin+ 2018 | `n_b`, `n`, `T`, `E` | single |
| `PlasmaInstabilityBret2010TwoStream` | Bret+ 2010 | `n_b`, `n`, `E` | single |
| `PlasmaInstabilityBret2010Filamentation` | Bret+ 2010 | `n_b`, `E` | single |
| `PlasmaInstabilityShalaby2020` | Shalaby+ 2020 | `n_b`, `n`, `E` | single |
| `PlasmaInstabilityMiniati2013` | Miniati & Elyiv 2013 | luminosity → `n_b`, `n`, `E` | growth rate |

## Choosing a model

- For direct comparisons among the **early effective prescriptions**, use
  `Broderick2012`, `Sironi2014`, and `Schlickeiser2012`.
- For **later literature fits**, add `Vafin2018` and `Shalaby2020`.
- To isolate a **specific instability channel**, compare `Bret2010TwoStream`
  against `Bret2010Filamentation`.
- Use `Miniati2013` when you prefer to specify a **source luminosity** and let the
  module derive the beam profile.

Because the prescriptions can disagree by orders of magnitude, the recommended
workflow is to run **several models over the same grid** of luminosity, density,
temperature, and redshift, and report the spread — not to rely on a single class.

## Control parameters (inherited)

Every model inherits two practical knobs from the base class, both optional
constructor arguments:

- **`efficiency`** (`η`, default `1`): multiplies the energy-loss rate.
  Clamped to `[0, 1]` — negative values become `0`, values above `1` become `1`
  (a warning is logged). Set to `0` to disable the instability, or to a fraction
  to test sensitivity.
- **`limit`** (default `0.1`): the next propagation step is restricted to this
  fraction of the local energy-loss length, so short cooling times are resolved.

```python
# 30% efficiency, tighter step control
plinst = PlasmaInstabilityVafin2018(beam, density, temperature, 0.3, 0.05)
```

## Helper functions

Two free functions are exposed for analytical estimates and are used internally by
`Miniati2013`:

- `plasmaFrequency(density, id=11)` — the plasma frequency
  `ω_p = sqrt(n q² / (m ε₀))` in Hz, for the given particle species (electron by
  default; nuclei use their charge and mass).
- `maximumLinearGrowthFrequency(beamDensity, mediumDensity, lorentzFactor, id=11)`
  — the maximum linear growth rate `ω_p(n) · (n_b/n) · (1/γ)` in Hz. Neglects
  magnetic fields and assumes an angular spread `Δθ = <1/γ>`.

Their scaling behaviour (`ω_p ∝ sqrt(n)`, growth rate linear in `n_b` and inverse
in `γ`) is checked in
[`test/testPlasmaInstability.cpp`](https://github.com/rafaelab/grplinst/blob/v2/test/testPlasmaInstability.cpp).

For the papers behind each model, see [References](references.html). For CRPropa
integration, continue with [Usage](usage.html).
