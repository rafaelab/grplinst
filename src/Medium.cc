#include "grplinst/Medium.h"


namespace grplinst {


/*****************************************************************************/
/*                             MediumTemperature                             */
/*****************************************************************************/

double MediumTemperature::getVelocity(int id, const crpropa::Vector3d& position, const double& redshift) const {
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


/*****************************************************************************/
/*                       MediumTemperatureHomogeneous                        */
/*****************************************************************************/

MediumTemperatureHomogeneous::MediumTemperatureHomogeneous(double T) : MediumTemperature() {
	setTemperatureValue(T);
}

void MediumTemperatureHomogeneous::setTemperatureValue(double T) {
	temperature = T;
}

double MediumTemperatureHomogeneous::getTemperatureValue() const {
	return temperature;
}

double MediumTemperatureHomogeneous::getTemperature(const crpropa::Vector3d& position, const double& redshift) const {
	return temperature * (1 + redshift);
}


/*****************************************************************************/
/*                          MediumTemperatureGrid                            */
/*****************************************************************************/

// MediumTemperatureGrid::MediumTemperatureGrid(const crpropa::ref_ptr<crpropa::Grid1f>& grid) : MediumTemperature() {
// 	setGrid(grid);
// }

// void MediumTemperatureGrid::setGrid(const crpropa::ref_ptr<crpropa::Grid1f>& g) {
// 	grid = g;
// }

// crpropa::ref_ptr<crpropa::Grid1f> MediumTemperatureGrid::getGrid() const {
// 	return grid;
// }

// double MediumTemperatureGrid::getTemperature(const crpropa::Vector3d& position, const double& redshift) const {
// 	return grid->getValue(position) * (1 + redshift);
// }


/*****************************************************************************/
/*                              MediumDensity                                */
/*****************************************************************************/



/*****************************************************************************/
/*                         MediumDensityHomogeneous                          */
/*****************************************************************************/

MediumDensityHomogeneous::MediumDensityHomogeneous(double n) : MediumDensity() {
	setDensityValue(n);
}

void MediumDensityHomogeneous::setDensityValue(double n) {
	density = n;
}

double MediumDensityHomogeneous::getDensityValue() const {
	return density;
}

double MediumDensityHomogeneous::getDensity(const crpropa::Vector3d& position, const double& redshift) const {
	return density * crpropa::pow_integer<3>(1 + redshift);
}


} // namespace grplinst