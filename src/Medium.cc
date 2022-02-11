#include "grplinst/Medium.h"


namespace grplinst {

HomogeneousMediumTemperature::HomogeneousMediumTemperature(const double T) {
	setTemperature(T);
}

void HomogeneousMediumTemperature::setTemperature(double T) {
	temperature = T;
}

double HomogeneousMediumTemperature::getTemperature() const {
	return temperature;
}

double HomogeneousMediumTemperature::getValue(crpropa::Vector3d position) const {
	return temperature;
}

double HomogeneousMediumTemperature::getValue(crpropa::Vector3d position, double redshift) const {
	return temperature * (1 + redshift);
}


//////////////////
HomogeneousMediumDensity::HomogeneousMediumDensity(double n) : MediumDensity() {
	setDensity(n);
}

void HomogeneousMediumDensity::setDensity(double n) {
	density = n;
}

double HomogeneousMediumDensity::getDensity() const {
	return density;
}

double HomogeneousMediumDensity::getValue(crpropa::Vector3d position) const {
	return density;
}

double HomogeneousMediumDensity::getValue(crpropa::Vector3d position, double redshift) const {
	return density * crpropa::pow_integer<3>(1 + redshift);
}

} // namespace grplinst