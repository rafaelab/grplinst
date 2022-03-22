#include "grplinst/Geometry.h"


namespace grplinst {

/***************************************************************************/
/**/
void EmissionGeometry::setShapeName(std::string s) {
	shapeName = s;
}

std::string EmissionGeometry::getShapeName() const {
	return shapeName;
}


/***************************************************************************/
/**/
Cone::Cone() {
	setAngle(0);
	setHeight(0);
	setRadius(computeRadius());
	setAxis(crpropa::Vector3d(-1, 0, 0));
	setShapeName("cone");
}

Cone::Cone(double angle, double height, crpropa::Vector3d axis) {
	setAngle(angle);
	setAxis(axis);
	setHeight(height);
	setRadius(computeRadius());
	setShapeName("cone");
}

Cone::~Cone() {
}

void Cone::setAngle(double a) {
	angle = a;
}

void Cone::setHeight(double h) {
	height = h;
}

void Cone::setRadius(double r) {
	radius = r;
}

void Cone::setAxis(crpropa::Vector3d a) {
	axis = a;
}

double Cone::getRadius() const {
	return radius;
}

double Cone::getAngle() const {
	return angle;
}

double Cone::getHeight() const {
	return height;
}

crpropa::Vector3d Cone::getAxis() const {
	return axis;
}

double Cone::computeRadius() const {
	return tan(angle) * height;
}

double Cone::computeArea() const {
	return M_PI * radius * (radius + sqrt(radius * radius + height * height));
}

double Cone::computeVolume() const {
	return M_PI * radius * radius * height / 3.;
}




} // namespace grplinst