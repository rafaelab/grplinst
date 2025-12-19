#ifndef GRPLINST_PLASMAINSTABILITY_H
#define GRPLINST_PLASMAINSTABILITY_H


#include <crpropa/Common.h>
#include <crpropa/Cosmology.h>
#include <crpropa/Module.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>

#include "grplinst/Medium.h"
#include "grplinst/Flow.h"


namespace grplinst {




/*****************************************************************************/
/*                   PlasmaInstability (abstract base class)                 */
/*****************************************************************************/

class PlasmaInstability : public crpropa::Module {
	protected:
		crpropa::ref_ptr<Flow> flowProperties;
		crpropa::ref_ptr<MediumDensity> mediumDensity;
		crpropa::ref_ptr<MediumTemperature> mediumTemperature;
		double limit;

	public:
		PlasmaInstability(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		virtual ~PlasmaInstability() = default;
		void setFlowProperties(crpropa::ref_ptr<Flow> flow);
		void setMediumDensity(crpropa::ref_ptr<MediumDensity> density);
		void setMediumTemperature(crpropa::ref_ptr<MediumTemperature> temperature);
		void setLimit(double limit);
		crpropa::ref_ptr<MediumDensity> getMediumDensity() const;
		crpropa::ref_ptr<Flow> getFlowProperties() const;
		crpropa::ref_ptr<MediumTemperature> getMediumTemperature() const;
		void process(crpropa::Candidate* candidate) const;
		double computeEnergyLossPerLength(crpropa::Candidate* candidate) const;
		virtual double energyLossTime(double energy, double beamDensity, double mediumDensity, double mediumTemperature) const = 0;
};

using PlasmaInstabilityPtr = std::unique_ptr<PlasmaInstability>;
// using PlasmaInstabilityRefPtr = crpropa::ref_ptr<PlasmaInstability>;



/*****************************************************************************/
/*                         Plasma Instability Models                         */
/*****************************************************************************/

/**
 * @brief Broderick et al. 2012 model for plasma instability energy loss time.
 * @see Broderick, Chang, Pfrommer. Astrophys. J. 752 (2012) 22.
 */
class PlasmaInstabilityBroderick2012 : public PlasmaInstability {
	public:
		using PlasmaInstability::PlasmaInstability;
		double energyLossTime(double energy, double beamDensity, double mediumDensity, double mediumTemperature) const override;
};

/**
 * @brief Schlickeiser et al. 2012 model for plasma instability energy loss time.
 * @see Schlickeiser, Ibscher, Supsar. Astrophys. J. 758 (2012) 102.
 */
class PlasmaInstabilitySchlickeiser2012 : public PlasmaInstability {
	public:
		using PlasmaInstability::PlasmaInstability;
		double energyLossTime(double energy, double beamDensity, double mediumDensity, double mediumTemperature) const override;
};

/**
 * @brief Sironi & Giannios 2014 model for plasma instability energy loss time.
 * @see Sironi, Giannios. Astrophys. J. 787 (2014) 49. arXiv:1312.4538
 */
class PlasmaInstabilitySironi2014 : public PlasmaInstability {
	public:
		using PlasmaInstability::PlasmaInstability;
		double energyLossTime(double energy, double beamDensity, double mediumDensity, double mediumTemperature) const override;
};

/**
 * @brief Vafin et al. 2018 model for plasma instability energy loss time.
 * @see Vafin, Pohl, Niemiec, Bret. Astrophys. J. 865 (2018) 23. arXiv:1807.04203
 */
class PlasmaInstabilityVafin2018 : public PlasmaInstability {
	public:
		using PlasmaInstability::PlasmaInstability;
		double energyLossTime(double energy, double beamDensity, double mediumDensity, double mediumTemperature) const override;
};

/**
 * @brief Bret et al. 2010 two-stream model for plasma instability energy loss time.
 * @see Bret, Gremillet, Dieckmann. Phys. Plasmas 17 (2010) 120501. arXiv:1010.5763
 */
class PlasmaInstabilityBret2010TwoStream : public PlasmaInstability {
	public:
		using PlasmaInstability::PlasmaInstability;
		double energyLossTime(double energy, double beamDensity, double mediumDensity, double mediumTemperature) const override;
};

/**
 * @brief Bret et al. 2010 filamentation model for plasma instability energy loss time.
 * @see Bret, Gremillet, Dieckmann. Phys. Plasmas 17 (2010) 120501. arXiv:1010.5763
 */
class PlasmaInstabilityBret2010Filamentation : public PlasmaInstability {
	public:
		using PlasmaInstability::PlasmaInstability;
		double energyLossTime(double energy, double beamDensity, double mediumDensity, double mediumTemperature) const override;
};

/**
 * @brief Shalaby et al. 2020 model for plasma instability energy loss time.
 * @see Shalaby et al. Phys. Rev. Lett. 124 (2020) 105101. arXiv:1907.13350
 */
class PlasmaInstabilityShalaby2020 : public PlasmaInstability {
	public:
		using PlasmaInstability::PlasmaInstability;
		double energyLossTime(double energy, double beamDensity, double mediumDensity, double mediumTemperature) const override;
};


/*****************************************************************************/
/*                             Helper functions                              */
/*****************************************************************************/

/**
 * @brief Computes the plasma frequency for a given density and particle type.
 * @param density Particle number density in m^-3.
 * @param id Particle ID (default: electron).
 * @return Plasma frequency in Hz.
 */
inline double plasmaFrequency(double density, int id = 11);


/**
 * @brief Computes the maximum linear growth frequency of the instability.
 * Neglects magnetic fields and assumes ΔΘ = <1/γ>.
 * @param beamDensity Beam particle number density in m^-3.
 * @param mediumDensity Medium particle number density in m^-3.
 * @param lorentzFactor Mean Lorentz factor of the beam particles.
 * @param id Particle ID of the beam particles (default: electron).
 * @return Maximum linear growth frequency in Hz.
 */
inline double maximumLinearGrowthFrequency(double beamDensity, double mediumDensity, double lorentzFactor, int id = 11);



} // namespace grplinst

#endif // GRPLINST_PLASMAINSTABILITY_H