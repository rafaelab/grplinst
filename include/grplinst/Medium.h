#ifndef GRPLINST_MEDIUM_H
#define GRPLINST_MEDIUM_H

#include <crpropa/Common.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>
#include <crpropa/Grid.h>
#include <crpropa/Referenced.h>


namespace grplinst {

/**
 @class MediumTemperature
 @brief Abstract base class to define temperature distributions in the medium.
 */
class MediumTemperature : public crpropa::Referenced {
	public:
		MediumTemperature() {
		}
		~MediumTemperature() {
		}
		virtual double getValue(crpropa::Vector3d position) const = 0;
		virtual double getValue(crpropa::Vector3d position, double redshift) const = 0;
};

/**
 @class MediumDensity
 @brief Abstract base class to define density distributions in the medium.
 */
class MediumDensity : public crpropa::Referenced {
	public:
		MediumDensity(){
		}
		~MediumDensity(){
		}
		virtual double getValue(crpropa::Vector3d position) const = 0;
		virtual double getValue(crpropa::Vector3d position, double redshift) const = 0;
};


/**
 @class HomogeneousMediumTemperature
 @brief Medium temperature is the same at all positions.
 */
class HomogeneousMediumTemperature : public MediumTemperature {
	private:
		double temperature;
	public:
		HomogeneousMediumTemperature(const double temperature);
		void setTemperature(double temperature);
		double getTemperature() const;
		double getValue(crpropa::Vector3d position) const;
		double getValue(crpropa::Vector3d position, double redshift) const;
};

/**
 @class HomogeneousMediumDensity
 @brief Medium density is the same at all positions.
 */
class HomogeneousMediumDensity : public MediumDensity {
	private:
		double density;
	public:
		HomogeneousMediumDensity(double density);
		void setDensity(double density);
		double getDensity() const; 
		double getValue(crpropa::Vector3d position) const;
		double getValue(crpropa::Vector3d position, double redshift) const;
};


} // namespace grplinst

#endif // GRPLINST_MEDIUM_H