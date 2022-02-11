#include "grplinst/Flow.h"


namespace grplinst {

HomogeneousFlow::HomogeneousFlow(double n) {
	setDensity(n);
}

void HomogeneousFlow::setDensity(double n) {
	density = n;
}

double HomogeneousFlow::getDensity() const {
	return density;
}

double HomogeneousFlow::getValue(crpropa::Vector3d position) const {
	return density;
}

double HomogeneousFlow::getValue(crpropa::Vector3d position, double redshift) const {
	return density * pow(1 + redshift, 3);
}

} // namespace grplinst