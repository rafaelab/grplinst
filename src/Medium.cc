#include "grplinst/Medium.h"


namespace grplinst {


MediumTemperature::MediumTemperature() {
}

MediumTemperature::~MediumTemperature() {
}

double MediumTemperature::getVelocity(int id, crpropa::Vector3d position, double redshift) const {
	// mass
	double m = 0;
	if (fabs(id) == 11) {
		m = crpropa::mass_electron;
	} else if (crpropa::isNucleus(id)) {
		m = crpropa::nuclearMass(id);
	} else {
		std::cerr << "Mass undefined for particle with id: " << id << std::endl;
	}
	
	double T = getValue(position, redshift);

	return sqrt(crpropa::k_boltzmann * T / m);
}

MediumTemperatureHomogeneous::MediumTemperatureHomogeneous(const double T) {
	setTemperature(T);
}

MediumTemperatureHomogeneous::~MediumTemperatureHomogeneous() {
}

void MediumTemperatureHomogeneous::setTemperature(double T) {
	temperature = T;
}

double MediumTemperatureHomogeneous::getTemperature() const {
	return temperature;
}

double MediumTemperatureHomogeneous::getValue(crpropa::Vector3d position, double redshift) const {
	return temperature * (1 + redshift);
}


//////////////////
MediumDensity::MediumDensity() {
}

MediumDensity::~MediumDensity() {
}

MediumDensityHomogeneous::MediumDensityHomogeneous(double n) : MediumDensity() {
	setDensity(n);
}

MediumDensityHomogeneous::~MediumDensityHomogeneous() {
	// MediumDensity::~MediumDensity();
}

void MediumDensityHomogeneous::setDensity(double n) {
	density = n;
}

double MediumDensityHomogeneous::getDensity() const {
	return density;
}

double MediumDensityHomogeneous::getValue(crpropa::Vector3d position, double redshift) const {
	return density * crpropa::pow_integer<3>(1 + redshift);
}

} // namespace grplinst