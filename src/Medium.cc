#include "grplinst/Medium.h"


namespace grplinst {


MediumTemperature::MediumTemperature() {
}

MediumTemperature::~MediumTemperature() {
}

double MediumTemperature::getVelocity(int id, const crpropa::Vector3d &position, double redshift) const {
	// mass
	double m = 0;
	if (fabs(id) == 11) {
		m = crpropa::mass_electron;
	} else if (crpropa::isNucleus(id)) {
		m = crpropa::nuclearMass(id);
	} else {
		std::cerr << "Mass undefined for particle with id: " << id << std::endl;
	}
	
	double T = getTemperature(position, redshift);

	return sqrt(crpropa::k_boltzmann * T / m);
}

MediumTemperatureHomogeneous::MediumTemperatureHomogeneous(double T) : MediumTemperature() {
	temperature = T;
}

MediumTemperatureHomogeneous::~MediumTemperatureHomogeneous() {
}

double MediumTemperatureHomogeneous::getTemperature(const crpropa::Vector3d &position, double redshift) const {
	return temperature * (1 + redshift);
}


//////////////////
MediumDensity::MediumDensity() {
}

MediumDensity::~MediumDensity() {
}

MediumDensityHomogeneous::MediumDensityHomogeneous(double n) : MediumDensity() {
	density = n;
}

MediumDensityHomogeneous::~MediumDensityHomogeneous() {
}

double MediumDensityHomogeneous::getDensity(const crpropa::Vector3d &position, double redshift) const {
	return density * crpropa::pow_integer<3>(1 + redshift);
}

} // namespace grplinst