---
layout: default
title: Physics Background
---

# Physics Background

This page explains *what* `grplinst` computes and *why*. 
It is meant to be read before you run anything, so that the numbers you get out carry a physical
meaning. 
None of it is needed to make the code compile, but it is needed to interpret the results honestly.

## 1. Blazar-induced electromagnetic cascades

TeV blazars emit very-high-energy (VHE) gamma rays. 
On their way to us, these photons do not travel through empty space: the Universe is filled, among other things, by the **extragalactic background light** (EBL) and the **cosmic microwave background** (CMB). 
A VHE photon can collide with a low-energy background photon and produce an electron–positron pair,
$$ \gamma + \gamma_\mathrm{bg} \;\rightarrow\; e^{+} + e^{-} . $$
The pairs are ultra-relativistic (Lorentz factors $$\gamma \sim 10^{6}$$) and inherit the direction of the parent gamma ray, forming a narrow **pair beam**.
Each pair then up-scatters CMB photons through **inverse-Compton (IC) scattering**,
$$ e^{\pm} + \gamma_\mathrm{bg} \;\rightarrow\; e^{\pm} + \gamma , $$
producing secondary gamma rays at GeV energies. 
Those can pair-produce again, and so on: the result is an **electromagnetic cascade** that reprocesses TeV power down to GeV energies.

The observable consequences —- a GeV "pair echo", a spatially extended "pair halo", and the reprocessed spectrum —- are sensitive to the **intergalactic magnetic field** (IGMF), because even a tiny field deflects the pairs. 
This makes blazar cascades one of the few probes of magnetic fields in cosmic voids. 
CRPropa's `EMPairProduction` and `EMInverseComptonScattering` modules implement exactly this chain.

## 2. The plasma-instability complication

The picture above assumes the pairs lose energy **only** through inverse-Compton scattering. 
But a relativistic pair beam streaming through the ionised IGM is a two-stream configuration, and such configurations are generically **unstable**:
the beam excites electrostatic (and other) plasma waves in the background medium, feeding its kinetic energy into those waves instead of into IC photons.

If this **beam–plasma instability** grows fast enough, it becomes a competitive (or even dominant) energy-loss channel for the pairs. 
The consequences are considerable:
- the pairs cool *before* they can up-scatter many CMB photons, so the **GeV cascade emission is suppressed**;
- constraints on the IGMF derived by attributing the missing GeV flux to magnetic deflection are **weakened or invalidated**, because the flux may be missing for a plasma reason instead;
- the energy deposited in the IGM may contribute to its **thermal history**.

Whether instabilities actually win over IC cooling is genuinely debated. 
The growth rate depends on the beam density, its angular and energy spread, the temperature and density of the background plasma, and non-linear saturation effects that linear theory does not capture. 
The literature spans "instabilities dominate" to "instabilities are irrelevant", which is precisely why a tool that lets you **swap between prescriptions** is useful.

## 3. How `grplinst` models it: effective cooling

`grplinst` does not simulate the instability itself.
Instead, it represents each literature prescription as an **effective cooling term** acting on the pairs, exactly like any other continuous energy-loss process in CRPropa.

Each model provides a single quantity, a **characteristic cooling time** $$\tau(E, n_b, n, T, z)$$, which the timescale on which a pair of energy $$E$$ loses
its energy to the instability. 
The base module turns this into an energy-loss *length* and applies it during propagation,
$$ \frac{\mathrm{d}E}{\mathrm{d}x} = \eta\,\frac{E}{c\,\tau} , $$
where $$c$$ is the speed of light and $$\eta \in [0, 1]$$ is a user-defined efficiency factor that lets you scale the effect (for sensitivity tests) or switch it off. 
For details, see [`PlasmaInstability::process`](api-reference.html#plasmainstability-base-class).

**Units convention.** The module interprets $$\tau$$ as a time in **seconds**. In fact, internally the code use S.I. units throughout.


Only **electrons and positrons** are affected ($$|\mathrm{id}| = 11$$); photons pass through untouched. 
This is the physically correct choice: the instability acts on the charged pair beam, not on the cascade photons.

## 4. The three physical inputs

Every prescription is a function of some subset of three environmental inputs, which you supply as separate, composable objects.

### 4.1 The pair-beam density — `Flow`

The single most important (and most uncertain) input is the **beam number density** ($$n_b$$). 
A `Flow` object answers the question *"how dense is the pair beam here?"* hrough `getDensity(energy, position, redshift)`.

The default estimate, `FlowHomogeneous`, follows Broderick et al. (2012). 
It assumes the source injects all of its luminosity $$L$$ into the beam and that IC cooling sets the pair lifetime, giving an **upper limit** of sorts on the beam density (their eq. 7):
$$ n_b(E, z) = \frac{L}{2\pi\,\lambda_{\gamma\gamma}^{3}\,\Gamma_\mathrm{IC}}\,\frac{1}{E} \,,$$
wherein
- $$\lambda_{\gamma\gamma}$$ is the **pair-production mean free path** of the parent photon on the EBL. With no CRPropa module attached, the code uses the
  analytic approximation
  $$ \lambda_{\gamma\gamma} \approx 35\ \mathrm{Mpc} \times \frac{0.5\,\mathrm{TeV}}{E} \times \left(\frac{1+z}{2}\right)^{-4.5} ; $$
  the factor $$0.5$$ reflects that a pair member carries about half the parent photon energy.
- $$\Gamma_\mathrm{IC}$$ is the **inverse-Compton loss rate** in the Thomson regime,
  $$ \Gamma_\mathrm{IC} = \frac{4}{3}\,\sigma_\mathrm{T}\,c\,u_\mathrm{CMB}\,\frac{E}{m_e c^2}\,\frac{(1+z)^{4}}{m_e c^2} . $$

If you instead attach CRPropa's `EMPairProduction` and `EMInverseComptonScattering` modules to the `Flow` (via `setPairProduction` / `setInverseCompton`), the corresponding **actual rates** are used in place of the analytic approximations, making the beam-density estimate consistent with the cross-sections and photon fields used in the rest of the simulation.

Because the beam density enters the cooling time strongly (see the exponents in [Models](models.html)), the choice of `Flow` model is often the dominant systematic. `FlowJet1D` lets you impose a tabulated $$n_b(r)$$ profile instead, and you can sub-class `Flow` in Python for anything else.

### 4.2 The ambient density — `MediumDensity`

$$n$$ is the number density of the background IGM plasma. 
It sets the local plasma frequency $$\omega_p = \sqrt{n e^2 / (m_e \varepsilon_0)}$$, which controls how the medium responds to the beam. `MediumDensityHomogeneous` returns a constant comoving value scaled to the proper density, $$n(z) = n_0 (1+z)^3$$.

### 4.3 The ambient temperature — `MediumTemperature`

$$T$$ is the temperature of the background plasma. 
It matters because the thermal spread of the background electrons can **quench** the instability: a hot medium is harder to destabilise. 
Models such as Schlickeiser (2012) and Vafin (2018) depend on $$T$$ explicitly. 
`MediumTemperatureHomogeneous` returns $$T(z) = T_0 (1+z)$$, and the base class can convert temperature into a thermal velocity $$v = \sqrt{k_B T / m}$$ through `getVelocity`.

## 5. Redshift conventions

`grplinst` follows CRPropa's convention that a candidate stores its **observed** (present-day) energy and that the environment is described in **comoving** terms.
Internally each model evaluates the physics in the local frame at the candidate's redshift $$z$$:

| Quantity | Scaling used in the code |
| --- | --- |
| Particle energy | $$E_\mathrm{local} = E_\mathrm{obs}\,(1+z)$$ |
| Medium density | $$n(z) = n_0\,(1+z)^3$$ |
| Beam density (`FlowJet1D`) | $$n_b(z) = n_b^\mathrm{com}\,(1+z)^3$$ |
| Medium temperature | $$T(z) = T_0\,(1+z)$$ |

The step length is likewise de-redshifted, $$\mathrm{d}x = \Delta x /(1+z)$$, before the loss is applied. You normally do not apply these factors yourself — you specify present-day/comoving values and the module handles the scaling.

## 6. What this model is *not*

Please keep the following limitations in mind:
- **It is not a kinetic calculation.** 
  The cooling times are analytic fits or order-of-magnitude estimates from the literature, each derived under its own assumptions (monochromatic vs. spread beams, cold vs. warm plasma, linear vs. saturated growth). 
  A first-principles answer needs particle-in-cell (PIC) simulations.
- **The prescriptions disagree.** Different models can give cooling times that differ by orders of magnitude for the same inputs. 
  That spread is a feature to be explored, not a bug; remember that `grplinst`'s goal is to compare models rather than trusting a single one.
- **The beam-density estimate is an upper limit.** `FlowHomogeneous` assumes 100% of the luminosity goes into the beam and IC-limited lifetimes; 
realistic densities are lower, which lengthens the cooling time.
- **Non-linear saturation is not modelled** except insofar as it is baked into a given fit.

Used with these caveats in mind, `grplinst` is a controlled way to ask *"how much could plasma instabilities change my cascade result, and how much does the answer
depend on which prescription I believe?"*

Continue to [Models](models.html) for the explicit cooling-time formulae, or to [Usage](usage.html) to build a simulation.
