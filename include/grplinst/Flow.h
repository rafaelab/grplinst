#ifndef GRPLINST_FLOW_H
#define GRPLINST_FLOW_H


#include <stdexcept>

#include <crpropa/Common.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>
#include <crpropa/module/EMPairProduction.h>
#include <crpropa/module/EMInverseComptonScattering.h>
#include <kiss/string.h>
#include <kiss/logger.h>

#include "grplinst/Common.h"

namespace grplinst {



/*****************************************************************************/
/*                         Flow (abstract base class)                        */
/*****************************************************************************/

/**
 * @class Flow
 * @brief Abstract base class to define properties related to the emitting object. 
 *
 * The Flow interface defines the spatial origin and provides a pure virtual method to query the local particle density.
 * Future versions of CRPropa will implement a `getRate` (or similar) method for `EMPairProduction` and `EMInverseComptonScattering` modules.
 * This will be used to self-consistently compute the beam density based on the pair production rate. 
 */
class Flow : public crpropa::Referenced {
	protected:
		std::vector<crpropa::ref_ptr<crpropa::EMPairProduction>> pairProduction = {};
		std::vector<crpropa::ref_ptr<crpropa::EMInverseComptonScattering>> inverseCompton = {};
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
		void setPairProduction(std::vector<crpropa::ref_ptr<crpropa::EMPairProduction>> pp);
		void addPairProduction(crpropa::ref_ptr<crpropa::EMPairProduction> pp);

		/**
		 * @brief Sets the inverse Compton scattering module to be used by CRPropa.
		 * @param ic Reference to the EMInverseComptonScattering module.
		 */
		void setInverseCompton(std::vector<crpropa::ref_ptr<crpropa::EMInverseComptonScattering>> ic);
		void addInverseCompton(crpropa::ref_ptr<crpropa::EMInverseComptonScattering> ic);

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
		virtual double getDensity(double energy, const crpropa::Vector3d& position, double redshift = 0) const = 0;

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
	public:
		 /** @brief Default constructor (density = 0, origin = (0,0,0)). */
		FlowHomogeneous();

		/**
		 * @brief Construct a homogeneous flow with given luminosity and optional origin.
		 * @param luminosity Constant luminosity.
		 * @param origin Flow origin position (default: (0,0,0)).
		 */
		FlowHomogeneous(double luminosity, crpropa::Vector3d origin = crpropa::Vector3d(0, 0, 0));


		/**  
		 * @brief Get the local particle density at a given position and redshift.
		 * @param energy Particle energy.
		 * @param position The position vector where the density is queried.
		 * @param redshift The redshift at which the density is evaluated (default is 0).
		 * @return The local particle density, computed based on luminosity and redshift.
		 * 
		 * This method computes the local particle density based on the luminosity of the flow and the redshift, following a specific scaling relation. 
		 * It assumes that all luminosity goes into beam particles and that Compton cooling dominates electron energy losses, providing an upper limit on the beam density.
		 * Ideally, one should compute the beam density self-consistently during the simulation or estimate it based on the pair production rate.
		*/
		double getDensity(double energy, const crpropa::Vector3d& position, double redshift = 0) const;
};


/*****************************************************************************/
/*                                  FlowJet1D                                */
/*****************************************************************************/

/**
 * @class FlowJet1D
 * @brief Implementation of a 1D jet flow with profiles defined along the jet axis.
 * The FlowJet1D class provides an implementation of the Flow interface for a one-dimensional jet structure.
 * It allows defining profiles, along the jet axis, for the density
 * Note all of these properties are used; just the infrastructure is implemented.
 * NOTE: This class has not been tested.
 */
class FlowJet1D : public Flow {
	protected:
		std::vector<double> distance;
		std::vector<double> densityProfile;
		bool interpolateLog = true;

	public:
		/**
		 * @brief Construct a FlowJet1D from explicit profiles.
		 * @param distances Vector of distances (monotonic) defining the profile grid.
		 * @param beamDensity Density values corresponding to distances.
		 * @param centre Origin/centre of the jet (default: (0,0,0))
		 * @param interpolateLog Whether to interpolate profiles logarithmically (default: true).
		 */
		FlowJet1D(const std::vector<double>& distances, const std::vector<double>& beamDensity, double luminosity = 1, crpropa::Vector3d centre = crpropa::Vector3d(0, 0, 0), bool interpolateLog = true);

		/** @brief Default constructor 
		 */
		FlowJet1D();


		/**
		 * @brief Set the density profile along the jet axis.
		 * @param density Vector of density values defining the profile grid.
		 */
		void setDensityProfile(const std::vector<double>& density);

		/**
		 * @brief Set the distance profile along the jet axis.
		 * @param distance Vector of distances defining the profile grid.
		 */
		void setDistanceProfile(const std::vector<double>& distance);

		/**
		 * @brief Set whether to interpolate profiles logarithmically.
		 * @param interpolateLog True to interpolate logarithmically, false for linear interpolation.
		 */
		void setInterpolateLog(bool interpolateLog);

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

		double getDensity(double energy, const crpropa::Vector3d& position, double redshift = 0) const;

};




/*****************************************************************************/
/*                                  Others                                   */
/*****************************************************************************/

crpropa::ref_ptr<Flow> createFlowMiniati2013(double luminosity, crpropa::Vector3d centre = crpropa::Vector3d(0, 0, 0), bool logDistance = true);


} // namespace grplinst


#endif // GRPLINST_FLOW_H
