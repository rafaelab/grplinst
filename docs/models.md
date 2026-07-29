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

Every formula below is written in terms of the following **scaled variables**, so the numbers match the source code directly. $$\tau$$ is in **seconds**.

| Symbol | Definition | Fiducial value |
| --- | --- | --- |
| $$\tilde{E}$$ | particle energy $$E / \mathrm{TeV}$$ (local frame, $$E = E_\mathrm{obs}(1+z)$$) | 1 |
| $$\tilde{n}_b$$ | beam density $$n_b / 10^{-16}\,\mathrm{m^{-3}}$$ (from the [`Flow`](api-reference.html#flow)) | 1 |
| $$\tilde{n}$$ | ambient density $$n / 0.1\,\mathrm{m^{-3}}$$ (from the [`MediumDensity`](api-reference.html#mediumdensity)) | 1 |
| $$\tilde{T}$$ | ambient temperature $$T / 10^{4}\,\mathrm{K}$$ (from the [`MediumTemperature`](api-reference.html#mediumtemperature)) | 1 |

A model uses only the inputs that appear in its formula; the others may be omitted (pass any placeholder object). 
Several models are **piecewise**, switching between a weak/linear branch and a saturated branch at a critical beam density $$n_\mathrm{crit}$$.

## Implemented prescriptions

### `PlasmaInstabilityBroderick2012`

Broderick, Chang & Pfrommer (2012). Piecewise in the beam density:

$$ n_\mathrm{crit} = \frac{1.6\times10^{-13}\,\tilde{n}}{\tilde{E}^{2}}\ \mathrm{m^{-3}} , $$

$$ \tau = \begin{cases} 7\times10^{7}\,\dfrac{\sqrt{\tilde{n}}}{\tilde{E}\,\tilde{n}_b} & n_b < n_\mathrm{crit} , \\[2.5ex] 5\times10^{5}\,\dfrac{\tilde{E}^{1/3}}{\tilde{n}_b^{1/3}\,\tilde{n}^{1/6}} & n_b \ge n_\mathrm{crit} . \end{cases} $$

Uses beam density, ambient density, and energy. The weak branch cools roughly as
$$1/(E\,n_b)$$; above $$n_\mathrm{crit}$$ the dependence flattens to cube roots.

### `PlasmaInstabilitySironi2014`

Sironi & Giannios (2014). Same structure as Broderick (2012), different constants:

$$ n_\mathrm{crit} = \frac{8\times10^{-14}\,\tilde{n}}{\tilde{E}^{2}}\ \mathrm{m^{-3}} , $$

$$ \tau = \begin{cases} 1.4\times10^{7}\,\dfrac{\sqrt{\tilde{n}}}{\tilde{E}\,\tilde{n}_b} & n_b < n_\mathrm{crit} , \\[2.5ex] 9.6\times10^{5}\,\dfrac{\tilde{E}^{1/3}}{\tilde{n}_b^{1/3}\,\tilde{n}^{1/6}} & n_b \ge n_\mathrm{crit} . \end{cases} $$

Uses beam density, ambient density, and energy.

### `PlasmaInstabilitySchlickeiser2012`

Schlickeiser, Ibscher & Supsar (2012). Piecewise, with an explicit temperature
dependence in both the threshold and the cooling time:

$$ n_\mathrm{crit} = \frac{2.5\times10^{-19}\,\tilde{n}\,\tilde{T}^{2}}{\tilde{E}}\ \mathrm{m^{-3}} , $$

$$ \tau = \begin{cases} 5\times10^{14}\,\dfrac{\tilde{E}^{5/3}\,\tilde{n}_b^{1/3}}{\tilde{n}^{5/6}\,\tilde{T}^{2}} & n_b < n_\mathrm{crit} , \\[2.5ex] 8\times10^{6}\,\dfrac{\tilde{E}^{1/3}}{\tilde{n}_b^{1/3}\,\tilde{n}^{1/6}}\,\Big[\,1 + 1.25\ln\tilde{T} - 0.25\ln\tilde{n}\,\Big] & n_b \ge n_\mathrm{crit} . \end{cases} $$

Uses beam density, ambient density, temperature, and energy.

### `PlasmaInstabilityVafin2018`

Vafin, Pohl, Niemiec & Bret (2018). A single regime with a strong inverse
temperature dependence:

$$ \tau = 1.9\times10^{11}\,\frac{\tilde{E}^{4/3}\,\tilde{n}^{1/3}}{\tilde{n}_b^{1/3}\,\tilde{T}} . $$

Uses beam density, ambient density, temperature, and energy. Because
$$\tau \propto 1/T$$, a hotter medium cools the beam faster in this prescription.

### `PlasmaInstabilityBret2010TwoStream`

Bret, Gremillet & Dieckmann (2010), two-stream channel:

$$ \tau = 1.6\times10^{10}\,\frac{\tilde{E}}{\tilde{n}_b^{1/3}\,\tilde{n}^{1/6}} . $$

Uses beam density, ambient density, and energy.

### `PlasmaInstabilityBret2010Filamentation`

Bret, Gremillet & Dieckmann (2010), filamentation channel:

$$ \tau = 2.5\times10^{9}\,\sqrt{\frac{\tilde{E}}{\tilde{n}_b}} . $$

Uses beam density and energy only. **Note:** although a `MediumDensity` and
`MediumTemperature` must be supplied to construct the module, the final
filamentation expression does not use them.

### `PlasmaInstabilityShalaby2020`

Shalaby et al. (2020), effective scaling from the linear-instability analysis:

$$ \tau = 3.2\times10^{12}\,\tilde{E}^{6/5}\,\tilde{n}_b^{-2/5}\,\tilde{n}^{-1/10} . $$

Uses beam density, ambient density, and energy.

### `PlasmaInstabilityMiniati2013`

Miniati & Elyiv (2013). Rather than a fitted cooling time, this model uses the
**inverse of the maximum linear growth rate**,

$$ \tau = \frac{1}{\omega_\mathrm{growth}} , \qquad \omega_\mathrm{growth} = \omega_p(n)\,\frac{n_b}{n}\,\frac{1}{\gamma} , $$

where $$\gamma = E / (m_e c^2)$$ and
$$\omega_p(n) = \sqrt{n e^2 / (m_e \varepsilon_0)}$$ is the plasma frequency of the
ambient medium (see [Helper functions](#helper-functions)).

This model is special in that it **builds its own beam density**: the constructor

```python
PlasmaInstabilityMiniati2013(luminosity, density, temperature)
```

calls `createFlowMiniati2013(luminosity)` internally, which sets up a 16-point
tabulated `FlowJet1D` profile (a numerical $$n_b(r)$$ from the paper) scaled
linearly with $$L / 10^{38}\,\mathrm{W}$$. You therefore pass a **luminosity**,
not a `Flow`. See [Usage](usage.html#luminosity-driven-miniati-setup).

## Summary table

| Class | Reference | Inputs used | Regimes |
| --- | --- | --- | --- |
| `PlasmaInstabilityBroderick2012` | Broderick+ 2012 | $$n_b, n, E$$ | piecewise |
| `PlasmaInstabilitySironi2014` | Sironi & Giannios 2014 | $$n_b, n, E$$ | piecewise |
| `PlasmaInstabilitySchlickeiser2012` | Schlickeiser+ 2012 | $$n_b, n, T, E$$ | piecewise |
| `PlasmaInstabilityVafin2018` | Vafin+ 2018 | $$n_b, n, T, E$$ | single |
| `PlasmaInstabilityBret2010TwoStream` | Bret+ 2010 | $$n_b, n, E$$ | single |
| `PlasmaInstabilityBret2010Filamentation` | Bret+ 2010 | $$n_b, E$$ | single |
| `PlasmaInstabilityShalaby2020` | Shalaby+ 2020 | $$n_b, n, E$$ | single |
| `PlasmaInstabilityMiniati2013` | Miniati & Elyiv 2013 | $$L \rightarrow n_b,\, n,\, E$$ | growth rate |

## Choosing a model

- For direct comparisons among the **early effective prescriptions**, use `Broderick2012`, `Sironi2014`, and `Schlickeiser2012`.
- For **later literature fits**, add `Vafin2018` and `Shalaby2020`.
- To isolate a **specific instability channel**, compare `Bret2010TwoStream` against `Bret2010Filamentation`.
- Use `Miniati2013` when you prefer to specify a **source luminosity** and let the module derive the beam profile.


## Control parameters (inherited)

Every model inherits two practical knobs from the base class, both optional constructor arguments:
- **`efficiency`** ($$\eta$$, default `1`): multiplies the energy-loss rate. Clamped to $$[0, 1]$$ — negative values become `0`, values above `1` become `1` (a warning is logged). Set to `0` to disable the instability, or to a fraction to test sensitivity.
- **`limit`** (default `0.1`): the next propagation step is restricted to this fraction of the local energy-loss length, so short cooling times are resolved.
```python
plinst = PlasmaInstabilityVafin2018(beam, density, temperature, 0.3, 0.05) # # 30% efficiency, tighter step control
```

## Helper functions

Two free functions are exposed for analytical estimates and are used internally by `Miniati2013`:
- `plasmaFrequency(density, id=11)`: the plasma frequency $$\omega_p = \sqrt{n q^2 / (m \varepsilon_0)}$$ in Hz, for the given particle species (electron by default; nuclei use their charge and mass).
- `maximumLinearGrowthFrequency(beamDensity, mediumDensity, lorentzFactor, id=11)`: the maximum linear growth rate $$\omega_p(n)\,(n_b/n)\,(1/\gamma)$$ in Hz.
Neglects magnetic fields and assumes an angular spread $$\Delta\theta = \langle 1/\gamma \rangle$$.

Their scaling behaviour ($$\omega_p \propto \sqrt{n}$$, growth rate linear in $$n_b$$ and inverse in $$\gamma$$) is checked in
[`test/testPlasmaInstability.cpp`](https://github.com/rafaelab/grplinst/blob/v2/test/testPlasmaInstability.cpp).

For the papers behind each model, see [References](references.html). 
For CRPropa integration, continue with [Usage](usage.html).
