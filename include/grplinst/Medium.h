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


/*****************************************************************************/
/*                             MediumTemperature                             */
/*****************************************************************************/

/**
 *  @class MediumTemperature
 *  @brief Abstract base class to define temperature distributions in the medium.
 */
class MediumTemperature : public crpropa::Referenced {
	public:
		virtual ~MediumTemperature() = default;
		virtual double getTemperature(const crpropa::Vector3d& position, const double& redshift = 0.) const = 0;
		double getVelocity(int id, const crpropa::Vector3d& position, const double& redshift = 0) const;
};



/*****************************************************************************/
/*                       MediumTemperatureHomogeneous                        */
/*****************************************************************************/

/**
 *  @class MediumTemperatureHomogeneous
 *  @brief Medium temperature is the same at all positions.
 */
class MediumTemperatureHomogeneous : public MediumTemperature {
	protected:
		double temperature;

	public:
		MediumTemperatureHomogeneous(double temperature);
		void setTemperatureValue(double T);
		double getTemperatureValue() const;
		double getTemperature(const crpropa::Vector3d& position = crpropa::Vector3d(0., 0., 0.), const double& redshift = 0.) const;
};


/*****************************************************************************/
/*                           MediumTemperatureGrid                           */
/*****************************************************************************/

// /**
//  *  @class MediumTemperatureGrid
//  *  @brief Medium temperature is defined on a grid.
//  */
// class MediumTemperatureGrid<T> : public MediumTemperature {
// 	protected:
// 		crpropa::ref_ptr<crpropa::Grid1f> grid;

// 	public:
// 		MediumTemperatureGrid(const crpropa::Grid3d<double>& grid);
// 		void setGrid(ref_ptr<Grid1f> grid);
// 		ref_ptr<Grid1f> getGrid();
// 		double getTemperature(const crpropa::Vector3d& position, const double& redshift = 0.) const;
// };


/*****************************************************************************/
/*                             MediumDensity                                 */
/*****************************************************************************/

/**
 *  @class MediumDensity
 *  @brief Abstract base class to define density distributions in the medium.
 */
class MediumDensity : public crpropa::Referenced {
	public:
		virtual ~MediumDensity() = default;
		virtual double getDensity(const crpropa::Vector3d& position, const double& redshift = 0.) const = 0;
};


/**
 *  @class MediumDensityHomogeneous
 *  @brief Medium density is the same at all positions.
 */
class MediumDensityHomogeneous : public MediumDensity {
	protected:
		double density;

	public:
		MediumDensityHomogeneous(double density);
		void setDensityValue(double n);
		double getDensityValue() const;
		double getDensity(const crpropa::Vector3d& position = crpropa::Vector3d(0., 0., 0.), const double& redshift = 0.) const;
};


} // namespace grplinst

#endif // GRPLINST_MEDIUM_H