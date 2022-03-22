#ifndef GRPLINST_GEOMETRY_H
#define GRPLINST_GEOMETRY_H


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
 @class EmissionGeometry
 @brief Abstract base class holding an arbitrary shape.
 */
class EmissionGeometry : public crpropa::Referenced {
	protected:
		std::string shapeName;
	public:
		// EmissionGeometry();
		// virtual ~EmissionGeometry() = default;
		void setShapeName(std::string shape);
		std::string getShapeName() const;
		virtual double computeArea() const = 0;
		virtual double computeVolume() const = 0;
		// static EmissionGeometry* materialise(std::string shape);
		// auto materialise();
};


/**
 @class Cone
 @brief Geometrical shape: cone (regular).
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
		~Cone();
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

// Cone* materialiseCone(EmissionGeometry* geometry);


} // namespace grplinst

#endif // GRPLINST_FLOW_H