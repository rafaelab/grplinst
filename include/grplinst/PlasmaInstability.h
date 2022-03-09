#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Cosmology.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Module.h>



#ifdef _OPENMP
    #include "omp.h"
#endif

#include "grplinst/Medium.h"
#include "grplinst/Flow.h"



namespace grplinst {

class PlasmaInstability : public crpropa::Module {
	protected:
		crpropa::ref_ptr<Flow> flowProperties;
		crpropa::ref_ptr<MediumDensity> mediumDensity;
		crpropa::ref_ptr<MediumTemperature> mediumTemperature;
		double limit;
	public:
		PlasmaInstability(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		~PlasmaInstability();
		void setFlowProperties(crpropa::ref_ptr<Flow> flow);
		void setMediumDensity(crpropa::ref_ptr<MediumDensity> density);
		void setMediumTemperature(crpropa::ref_ptr<MediumTemperature> temperature);
		void setLimit(double limit);
		crpropa::ref_ptr<MediumDensity> getMediumDensity() const;
		crpropa::ref_ptr<Flow> getFlowProperties() const;
		crpropa::ref_ptr<MediumTemperature> getMediumTemperature() const;
		void process(crpropa::Candidate *candidate) const;
		virtual double energyLoss(crpropa::Candidate *candidate) const = 0;
};

/**********************************************/

// Broderick, Chang, Pfrommer. Astrophys. J. 752 (2012) 22. arXiv:1106.5494
class PlasmaInstabilityBroderick2012B : public PlasmaInstability {
	public:
		PlasmaInstabilityBroderick2012B(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

//  Schlickeiser, Ibscher, Supsar. Astrophys. J. 758 (2012) 102.
class PlasmaInstabilitySchlickeiser2012B : public PlasmaInstability {
	public:
		PlasmaInstabilitySchlickeiser2012B(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

// Miniati & Elyiv. Astrophys. J. 770 (2013) 54. arXiv:1208.1761
class PlasmaInstabilityMiniati2013B : public PlasmaInstability {
	public:
		PlasmaInstabilityMiniati2013B(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

//  Vafin, Rafighi, Pohl, Niemiec. Astrophys. J. 857 (2018) 43. arXiv:1803.02990
class PlasmaInstabilityVafin2018B : public PlasmaInstability {
	public:
		PlasmaInstabilityVafin2018B(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};


/**********************************************/
// Things here follow Mahmoud's implementation

// Broderick, Chang, Pfrommer. Astrophys. J. 752 (2012) 22. arXiv:1106.5494
class PlasmaInstabilityBroderick2012C : public PlasmaInstability {
	public:
		PlasmaInstabilityBroderick2012C(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

//  Schlickeiser, Ibscher, Supsar. Astrophys. J. 758 (2012) 102.
class PlasmaInstabilitySchlickeiser2012C : public PlasmaInstability {
	public:
		PlasmaInstabilitySchlickeiser2012C(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

//  Vafin, Rafighi, Pohl, Niemiec. Astrophys. J. 857 (2018) 43. arXiv:1803.02990
// As provided by Mahmoud (might be wrong).
class PlasmaInstabilityVafin2018C : public PlasmaInstability {
	public:
		PlasmaInstabilityVafin2018C(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

// Miniati & Elyiv. Astrophys. J. 770 (2013) 54. arXiv:1208.1761
class PlasmaInstabilityMiniati2013C : public PlasmaInstability {
	public:
		PlasmaInstabilityMiniati2013C(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

/**********************************************/

// Broderick, Chang, Pfrommer. Astrophys. J. 752 (2012) 22. arXiv:1106.5494
class PlasmaInstabilityBroderick2012 : public PlasmaInstability {
	public:
		PlasmaInstabilityBroderick2012(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

// Miniati & Elyiv. Astrophys. J. 770 (2013) 54. arXiv:1208.1761
class PlasmaInstabilityMiniati2013 : public PlasmaInstability {
	private:
		std::vector<double> _w;
		std::vector<double> _d;
	public:
		PlasmaInstabilityMiniati2013(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
		void initTable();
};

//  Schlickeiser, Ibscher, Supsar. Astrophys. J. 758 (2012) 102.
class PlasmaInstabilitySchlickeiser2012 : public PlasmaInstability {
	public:
		PlasmaInstabilitySchlickeiser2012(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

//  Sironi, Giannios. Astrophys. J. 787 (2014) 49. arXiv:1312.4538
class PlasmaInstabilitySironi2014 : public PlasmaInstability {
	public:
		PlasmaInstabilitySironi2014(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};

//  Vafin, Rafighi, Pohl, Niemiec. Astrophys. J. 857 (2018) 43. arXiv:1803.02990
class PlasmaInstabilityVafin2018 : public PlasmaInstability {
	public:
		PlasmaInstabilityVafin2018(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit = 0.1);
		double energyLoss(crpropa::Candidate *candidate) const;
};


// Helper function to compute the plasma frequency for a given density and particle type.
double plasmaFrequency(double density, int id = 11);

// Computes the maximum linear growth rate (in the absence of magnetic fields).
// Assumes ΔΘ = <1/γ>.
double maximumLinearGrowthFrequency(double beamDensity, double mediumDensity, double lorentzFactor, int id = 11);

} // namespace grplinst