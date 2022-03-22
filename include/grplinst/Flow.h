#ifndef GRPLINST_FLOW_H
#define GRPLINST_FLOW_H


#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <crpropa/Candidate.h>
#include <crpropa/Common.h>
#include <crpropa/Cosmology.h>
#include <crpropa/Grid.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>

namespace grplinst {

/**
 @class Beam
 @brief Abstract base class to define properties related to the emitting object.
 */
class Flow : public crpropa::Referenced {
	protected:
		crpropa::Vector3d origin;
	public:
		Flow(crpropa::Vector3d origin = crpropa::Vector3d(0, 0, 0));
		~Flow();
		void setOrigin(const crpropa::Vector3d &origin);
		crpropa::Vector3d getOrigin() const;
		virtual double getDensity(const crpropa::Vector3d &position, double redshift = 0) const = 0;
		virtual double getMeanLorentzFactor(const crpropa::Vector3d &position, double redshift = 0, double lorentzFactorParticle = 1) const = 0;
		virtual double getMeanInverseLorentzFactor(const crpropa::Vector3d &position, double redshift = 0, double lorentzFactorParticle = 1) const = 0;
};

/**
 @class FlowHomogeneous
 @brief Class to describe the  
 */
class FlowHomogeneous : public Flow {
	private:
		double density;
	public:
		FlowHomogeneous(double density, crpropa::Vector3d origin = crpropa::Vector3d(0, 0, 0));
		FlowHomogeneous();
		~FlowHomogeneous();
		void setTotalDensity(double density);
		double getTotalDensity() const;
		double getDensity(const crpropa::Vector3d &position, double redshift = 0) const;
		double getMeanLorentzFactor(const crpropa::Vector3d &position, double redshift = 0, double lorentzFactorParticle = 1) const;
		double getMeanInverseLorentzFactor(const crpropa::Vector3d &position, double redshift = 0, double lorentzFactorParticle = 1) const;
};


/**
 */
class FlowJet1D : public Flow {
	friend class FlowJet1DMahmoud;
	// friend class FlowJetSimulation;
	protected:
		double densityNormalisation;
		std::vector<double> distance;
		std::vector<double> densityProfile;
		std::vector<double> meanLorentzFactor;
		std::vector<double> meanInverseLorentzFactor;
	public:
		FlowJet1D(const std::vector<double> &distances, const std::vector<double> &beamDensity, const std::vector<double> &lorentzFactor, const std::vector<double> &inverseLorentzFactor, double densityNorm = 1, crpropa::Vector3d centre = crpropa::Vector3d(0, 0, 0));
		FlowJet1D(const std::string &filename, double densityNormalisation = 1, crpropa::Vector3d origin = crpropa::Vector3d(0, 0, 0));
		FlowJet1D();
		~FlowJet1D();
		void setDensityNormalisation(double densityNorm);
		void setDensityProfile(const std::vector<double> &density);
		void setDistanceProfile(const std::vector<double> &distance);
		void setLorentzFactorProfile(const std::vector<double> &lorentzFactor);
		void setInverseLorentzFactorProfile(const std::vector<double> &inverseLorentzFactor);
		double getDensityNormalisation() const;
		std::vector<double> getDistanceProfile() const;
		std::vector<double> getDensityProfile() const;
		std::vector<double> getLorentzFactorProfile() const;
		std::vector<double> getInverseLorentzFactorProfile() const;
		double getDensity(const crpropa::Vector3d &position, double redshift = 0) const;
		double getMeanLorentzFactor(const crpropa::Vector3d &position, double redshift = 0, double lorentzFactorParticle = 1) const;
		double getMeanInverseLorentzFactor(const crpropa::Vector3d &position, double redshift = 0, double lorentzFactorParticle = 1) const;
};


/**
 */
class FlowJet1DMahmoud : public FlowJet1D {
	public:
		FlowJet1DMahmoud(double densityNormalisation = 1, crpropa::Vector3d origin = crpropa::Vector3d(0, 0, 0));
		~FlowJet1DMahmoud();
};



/**
 */
// void analyseSimulation(std::string filename, double cutoffEnergy, double spectralIndex, double luminosity, double spectralIndexSimulation = 1);




/**
Converts a given density to luminosity and vice-versa.
This is essentially eq. 13 from:
  R. Alves Batista, A. Saveliev, E. M. de Gouveia Dal Pino. 
  Mon. Not. R. Astron. Soc. 489 (2019) 3836-3849.
  arXiv:1904.13345
which is based on:
  Broderick et al. 
  Astrophys. J. 752 (2012) 22. 
  arXiv:1106.5494
This is an upper limit based on analytical estimates of cascade development assuming inverse Compton,
and considering the plasma cooling dominated by the kinetic oblique mode.
 */

inline double beamDensityToLuminosity(double n, double Ee, double z, double a = 4.5) {
	return 1e38 * (n / 3.7e-16) * (crpropa::TeV / Ee) * pow((1. + z) / 2., 4. - 3. * a);
}

inline double beamLuminosityToDensity(double L, double Ee, double z, double a = 4.5) {
	return 3.7e-16 * (L / 1e38) * (Ee / crpropa::TeV) * pow((1. + z) / 2., 3. * a - 4.);
}

} // namespace grplinst

#endif // GRPLINST_FLOW_H