#ifndef GRPLINST_MEDIUM_H
#define GRPLINST_MEDIUM_H

#include <crpropa/Common.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>
#include <crpropa/Grid.h>
#include <crpropa/Referenced.h>
#include <crpropa/ParticleId.h>
#include <crpropa/ParticleMass.h>


namespace grplinst {

/**
 @class MediumTemperature
 @brief Abstract base class to define temperature distributions in the medium.
 */
class MediumTemperature : public crpropa::Referenced {
	public:
		MediumTemperature();
		~MediumTemperature();
		double getVelocity(int id, const crpropa::Vector3d &position, double redshift = 0) const;
		virtual double getTemperature(const crpropa::Vector3d &position, double redshift = 0.) const = 0;
};

/**
 @class MediumDensity
 @brief Abstract base class to define density distributions in the medium.
 */
class MediumDensity : public crpropa::Referenced {
	public:
		MediumDensity();
		~MediumDensity();
		virtual double getDensity(const crpropa::Vector3d &position, double redshift = 0.) const = 0;
};


/**
 @class MediumTemperatureHomogeneous
 @brief Medium temperature is the same at all positions.
 */
class MediumTemperatureHomogeneous : public MediumTemperature {
	private:
		double temperature;
	public:
		MediumTemperatureHomogeneous(double temperature);
		~MediumTemperatureHomogeneous();
		double getTemperature(const crpropa::Vector3d &position, double redshift = 0.) const;
};

/**
 @class MediumDensityHomogeneous
 @brief Medium density is the same at all positions.
 */
class MediumDensityHomogeneous : public MediumDensity {
	private:
		double density;
	public:
		MediumDensityHomogeneous(double density);
		~MediumDensityHomogeneous();
		double getDensity(const crpropa::Vector3d &position, double redshift = 0.) const;
};


} // namespace grplinst

#endif // GRPLINST_MEDIUM_H