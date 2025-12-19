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
#include <crpropa/module/EMPairProduction.h>
#include <crpropa/module/EMInverseComptonScattering.h>

namespace grplinst {



/*****************************************************************************/
/*                         Flow (abstract base class)                        */
/*****************************************************************************/

/**
 * @class Flow
 * @brief Abstract base class to define properties related to the emitting object. 
 *
 * The Flow interface defines the spatial origin and provides pure virtual methods to query the local particle density.
 * It also computes the moments of the Lorentz factor distribution at a given position and redshift, for some specific classes.
 */
class Flow : public crpropa::Referenced {
	protected:
		crpropa::ref_ptr<crpropa::EMPairProduction> pairProduction;
		crpropa::ref_ptr<crpropa::EMInverseComptonScattering> inverseCompton;
		crpropa::Vector3d origin = crpropa::Vector3d(0, 0, 0);
		double luminosity = 1.;

	public:
		/** @brief Virtual default destructor. 
		*/
		virtual ~Flow() = default;

		/**
		 * @brief Sets the origin of the flow.
		 * @param centre The position vector representing the origin of the flow.
		 */
		void setOrigin(const crpropa::Vector3d& origin);

		/**
		 * @brief Sets the luminosity of the flow.
		 * @param lum The luminosity value to set.
		 */
		void setLuminosity(double lum);

		/**
		 * @brief Sets the pair production module to be used by CRPropa.
		 * @param pp Reference to the EMPairProduction module.
		 */
		void setPairProduction(crpropa::ref_ptr<crpropa::EMPairProduction> pp);

		/**
		 * @brief Sets the inverse Compton scattering module to be used by CRPropa.
		 * @param ic Reference to the EMInverseComptonScattering module.
		 */
		void setInverseCompton(crpropa::ref_ptr<crpropa::EMInverseComptonScattering> ic);

		/**
		 * @brief Gets the origin of the flow.
		 * @return The position vector representing the origin of the flow.
		 */
		crpropa::Vector3d getOrigin() const;

		/**
		 * @brief Gets the luminosity of the flow.
		 * @return The luminosity value of the flow.
		 */
		double getLuminosity() const;

		/**
		 * @brief Gets the local particle density at a given position and redshift.
		 * @param position The position vector where the density is queried.
		 * @param redshift The redshift at which the density is evaluated (default is 0).
		 * @return The local particle density.
		 */
		virtual double getDensity(const crpropa::Vector3d& position, double redshift = 0) const = 0;

		/**
		 * @brief Gets the mean Lorentz factor of the flow at a given position and redshift.
		 * @param position The position vector where the mean Lorentz factor is queried.
		 * @param redshift The redshift at which the mean Lorentz factor is evaluated (default is 0).
		 * @param lorentzFactorParticle The Lorentz factor of the particle (default is 1).
		 * @return The mean Lorentz factor of the flow.
		 */
		virtual double getMeanLorentzFactor(const crpropa::Vector3d& position, double redshift = 0, double lorentzFactorParticle = 1) const = 0;

		/**
		 * @brief Gets the mean inverse Lorentz factor of the flow at a given position and redshift.
		 * @param position The position vector where the mean inverse Lorentz factor is queried.
		 * @param redshift The redshift at which the mean inverse Lorentz factor is evaluated (default is 0).
		 * @param lorentzFactorParticle The Lorentz factor of the particle (default is 1).
		 * @return The mean inverse Lorentz factor of the flow.
		 */
		virtual double getMeanInverseLorentzFactor(const crpropa::Vector3d& position, double redshift = 0, double lorentzFactorParticle = 1) const = 0;

	// protected:
		// double getPairProductionMeanFreePath(double energy, double redshift = 0) const;
};



/*****************************************************************************/
/*                               FlowHomogeneous                             */
/*****************************************************************************/

/**
 * @class FlowHomogeneous
 * @brief Implementation of a homogeneous flow with constant density everywhere.
 * 
 * The FlowHomogeneous class provides a simple implementation of the Flow interface, assuming a constant density throughout space.
*/
class FlowHomogeneous : public Flow {
	private:
		double density = 0.;

	public:
		 /** @brief Default constructor (density = 0, origin = (0,0,0)). */
		FlowHomogeneous();

		/**
		 * @brief Construct a homogeneous flow with given density and optional origin.
		 * @param density Constant particle density.
		 * @param origin Flow origin position (default: (0,0,0)).
		 */

		FlowHomogeneous(double density, crpropa::Vector3d origin = crpropa::Vector3d(0, 0, 0));

		/**
		 * @brief Set the (total) density of the homogeneous flow.
		 * @param density New constant density value.
		 */
		void setDensityValue(double density);

		/**
		 * @brief Get the (total) density of the homogeneous flow.
		 * @return The constant density value.
		 */
		double getDensityValue() const;


		/**  */
		double getDensity(const crpropa::Vector3d& position, double redshift = 0) const;
		double getMeanLorentzFactor(const crpropa::Vector3d& position, double redshift = 0, double lorentzFactorParticle = 1) const;
		double getMeanInverseLorentzFactor(const crpropa::Vector3d& position, double redshift = 0, double lorentzFactorParticle = 1) const;


		/**  
		 * @brief Simple estimate of the beam density based on luminosity.
		 * 
		 * This basically follows equation 7 of:
		 *   "The Cosmological Impact of Luminous TeV Blazars. I. Implications of Plasma Instabilities for the Intergalactic Magnetic Field and Extragalactic Gamma-Ray Background"
		 *   A. E. Broderick, P. Chang, C. Pfrommer
		 *   The Astrophysical Journal 752 (2012) 22.
		 *   doi:10.1088/0004-637X/752/1/22
		 *   arXiv:1206.0731
		 * This is an upper limit, assuming all luminosity goes into beam particles, and that Compton cooling dominates electron energy losses.
		 * Ideally one should compute the beam density self-consistently during the simulation, or at least estimate it based on the pair production rate.
		 * 
		 * @param energy Particle energy.
		 * @param redshift Redshift.
		 * @return Estimated beam density.
		*/
		double estimateBeamDensity(double energy, double redshift) const;
};


/*****************************************************************************/
/*                                  FlowJet1D                                */
/*****************************************************************************/

/**
 * @class FlowJet1D
 * @brief Implementation of a 1D jet flow with profiles defined along the jet axis.
 * The FlowJet1D class provides an implementation of the Flow interface for a one-dimensional jet structure.
 * It allows defining profiles for density, Lorentz factor, and inverse Lorentz factor along the jet axis.
 * 
 * NOTE: This class has not been tested.
 */
class FlowJet1D : public Flow {
	protected:
		double densityNormalisation = 1.;
		std::vector<double> distance;
		std::vector<double> densityProfile;
		std::vector<double> meanLorentzFactor;
		std::vector<double> meanInverseLorentzFactor;

	public:
		/**
		 * @brief Construct a FlowJet1D from explicit profiles.
		 * @param distances Vector of distances (monotonic) defining the profile grid.
		 * @param beamDensity Density values corresponding to distances.
		 * @param lorentzFactor Mean Lorentz factor values corresponding to distances.
		 * @param inverseLorentzFactor Mean inverse Lorentz factor values corresponding to distances.
		 * @param densityNorm Global density normalisation factor (default 1).
		 * @param centre Origin/centre of the jet (default: (0,0,0)).
		 */
		FlowJet1D(const std::vector<double>& distances, const std::vector<double>& beamDensity, const std::vector<double>& lorentzFactor, const std::vector<double>& inverseLorentzFactor, double densityNorm = 1, crpropa::Vector3d centre = crpropa::Vector3d(0, 0, 0));

		/**
		 * @brief Construct a FlowJet1D from a data file.
		 * @param filename Path to the input file containing the profiles.
		 * @param densityNormalisation Global density normalisation factor (default 1).
		 * @param origin Origin/centre of the jet (default: (0,0,0)).
		 * 
		 * The input file should contain columns with distance, density, Lorentz factor, and inverse Lorentz factor values.
		 */
		FlowJet1D(const std::string &filename, double densityNormalisation = 1, crpropa::Vector3d origin = crpropa::Vector3d(0, 0, 0));

		/** @brief Default constructor (densityNorm = 1, origin = (0,0,0)). 
		 */
		FlowJet1D();

		/**
		 * @brief Set the global density normalisation.
		 * @param densityNorm Multiplicative normalisation factor.
		 */
		void setDensityNormalisation(double densityNorm);

		/**
		 * @brief Get the global density normalisation.
		 * @return The multiplicative normalisation factor.
		 */
		void setDensityProfile(const std::vector<double>& density);

		/**
		 * @brief Set the distance profile along the jet axis.
		 * @param distance Vector of distances defining the profile grid.
		 */
		void setDistanceProfile(const std::vector<double>& distance);

		/**
		 * @brief Set the Lorentz factor profile along the jet axis.
		 * @param lorentzFactor Vector of mean Lorentz factor values.
		 */
		void setLorentzFactorProfile(const std::vector<double>& lorentzFactor);

		/**
		 * @brief Set the inverse Lorentz factor profile along the jet axis.
		 * @param inverseLorentzFactor Vector of mean inverse Lorentz factor values.
		 */
		void setInverseLorentzFactorProfile(const std::vector<double>& inverseLorentzFactor);

		/** 
		 * @brief Get the global density normalisation.
		 * @return The multiplicative normalisation factor.
		 */
		double getDensityNormalisation() const;

		/** 
		 * @brief Get the distance profile along the jet axis.
		 * @return Vector of distances defining the profile grid.
		 */
		std::vector<double> getDistanceProfile() const;

		/** 
		 * @brief Get the density profile along the jet axis.
		 * @return Vector of density values defining the profile grid.
		 */
		std::vector<double> getDensityProfile() const;

		/** 
		 * @brief Get the Lorentz factor profile along the jet axis.
		 * @return Vector of mean Lorentz factor values.
		 */
		std::vector<double> getLorentzFactorProfile() const;

		/** 
		 * @brief Get the inverse Lorentz factor profile along the jet axis.
		 * @return Vector of mean inverse Lorentz factor values.
		 */
		std::vector<double> getInverseLorentzFactorProfile() const;

		/**  */
		double getDensity(const crpropa::Vector3d& position, double redshift = 0) const;
		double getMeanLorentzFactor(const crpropa::Vector3d& position, double redshift = 0, double lorentzFactorParticle = 1) const;
		double getMeanInverseLorentzFactor(const crpropa::Vector3d& position, double redshift = 0, double lorentzFactorParticle = 1) const;
};



/*****************************************************************************/
/*                             Other definitions                             */
/*****************************************************************************/




} // namespace grplinst


#endif // GRPLINST_FLOW_H