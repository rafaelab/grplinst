#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Cosmology.h>
#include <crpropa/Referenced.h>

#ifdef _OPENMP
    #include "omp.h"
#endif

#include "grplinst/Medium.h"
#include "grplinst/Flow.h"



namespace grplinst {

class PlasmaInstability : public crpropa::Module {
	private:
		crpropa::ref_ptr<MediumDensity> mediumDensity;
		crpropa::ref_ptr<MediumTemperature> mediumTemperature;
		crpropa::ref_ptr<Flow> flowProperties;
		std::string plasmaInstabilityModel;
		double limit;
	public:
		PlasmaInstability(MediumDensity *density, MediumTemperature *temperature, Flow *flow, const std::string instabilityModel, double limit = 0.1);
		~PlasmaInstability();
		void setMediumDensity(MediumDensity *density);
		void setMediumTemperature(MediumTemperature *temperature);
		void setFlowProperties(Flow *flow);
		void setInstabilityModel(const std::string model);
		void setLimit(double limit);
		crpropa::ref_ptr<MediumDensity> getMediumDensity() const;
		crpropa::ref_ptr<MediumTemperature> getMediumTemperature() const;
		crpropa::ref_ptr<Flow> getFlowProperties() const;
		std::string getInstabilityModel() const;
		void process(crpropa::Candidate *candidate) const;
		virtual double energyLoss(crpropa::Candidate *candidate) const = 0;
};


// class PlasmaInstabilityBroderick2012 : public PlasmaInstability {

// }



} // namespace grplinst