#ifndef GRPLINST_GEOMETRY_H
#define GRPLINST_GEOMETRY_H


#include <cmath>
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


/*****************************************************************************/
/*                   EmissionGeometry (abstract base class)                  */
/*****************************************************************************/

/**
 * @class EmissionGeometry
 * @brief Abstract base class holding an arbitrary shape.
*/
class EmissionGeometry : public crpropa::Referenced {
	protected:
		std::string shapeName;

	public:
		virtual ~EmissionGeometry() = default;
		void setShapeName(std::string shape);
		std::string getShapeName() const;
		virtual double computeArea() const = 0;
		virtual double computeVolume() const = 0;
		// static EmissionGeometry* materialise(std::string shape);
		// auto materialise();
};




/*****************************************************************************/
/*                          Cone (<: EmissionGeometry)                       */
/*****************************************************************************/

/**
 * @class Cone
 * @brief Geometrical shape: cone (regular).
 * Here we define a cone by its axis, angle, height, and radius:
 * - axis: direction vector of the cone axis
 * - angle: opening angle of the cone (in radians)
 * - height: height of the cone
 */
class Cone : public EmissionGeometry {
	protected:
		crpropa::Vector3d axis;
		double angle;
		double height;
		double radius;

	public:
		Cone();
		Cone(double angle, double height, crpropa::Vector3d axis = crpropa::Vector3d(-1, 0, 0));
		~Cone() override = default;
		void setAngle(double angle);
		void setHeight(double height);
		void setRadius(double radius);
		void setAxis(crpropa::Vector3d axis);
		double getHeight() const;
		double getAngle() const;
		double getRadius() const;
		crpropa::Vector3d getAxis() const;
		double computeRadius() const;
		double computeArea() const override;
		double computeVolume() const override;
};



} // namespace grplinst

#endif // GRPLINST_GEOMETRY_H