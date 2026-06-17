#include "grplinst/Flow.h"


namespace grplinst {



/*****************************************************************************/
/*                        Flow (abstract base class)                         */
/*****************************************************************************/

void Flow::setOrigin(const crpropa::Vector3d& centre) {
	origin = centre;
}

void Flow::setLuminosity(double lum) {
	luminosity = lum;
}

void Flow::setPairProduction(crpropa::ref_ptr<crpropa::EMPairProduction> pp) {
	pairProduction = pp;
}

void Flow::setInverseCompton(crpropa::ref_ptr<crpropa::EMInverseComptonScattering> ic) {
	inverseCompton = ic;
}

crpropa::Vector3d Flow::getOrigin() const {
	return origin;
}

double Flow::getLuminosity() const {
	return luminosity;
}

// double Flow::getPairProductionMeanFreePath(double energy, double redshift) const {
// 	if (pairProduction != nullptr) {
// 		return pairProduction->getRate(energy, redshift);
// 	}
// }

// double FlowHomogeneous::estimateBeamDensity(double E, double z) const {
// 	return 3.7e-16 * pow((1. + z) / 2., 9.5) * E * (luminosity / 1e38) * (E / crpropa::TeV);
// 	// double d = crpropa::redshift2LightTravelDistance(z);
// 	// double rIC = 1.2 * crpropa::pow_integer<3>(1 + z) * crpropa::kpc;
// 	// double L = luminosity;
// 	// return L / (2 * M_PI * crpropa::pow_integer<3>(d) * rIC) / E;
// }



/*****************************************************************************/
/*                              FlowHomogeneous                              */
/*****************************************************************************/

FlowHomogeneous::FlowHomogeneous() {
}

FlowHomogeneous::FlowHomogeneous(double L, crpropa::Vector3d centre) { 
	setOrigin(centre);
	setLuminosity(L);
}

double FlowHomogeneous::getDensity(double energy, const crpropa::Vector3d& position, double redshift) const {

	// pair-production mean free path on EBL (Broderick+2012 eq. 1), photon energy = 2 * electron energy
	// factor 0.5 in energy comes from the average energy of the parent
	double lambdaPP = 35. * crpropa::Mpc * (0.5 * crpropa::TeV / energy) * pow((1. + redshift) / 2., -4.5);

	// IC energy-loss rate in Thomson regime
	double GammaIC = (4. / 3.) * crpropa::sigma_thomson * crpropa::c_light * u_CMB * (energy / mec2) * pow(1. + redshift, 4.) / mec2;

	// eq. 7 of Broderick et al. 2012
	return luminosity / (2. * M_PI * crpropa::pow_integer<3>(lambdaPP) * GammaIC) / energy;
}

// double FlowHomogeneous::getDensity(double energy, const crpropa::Vector3d& position, double redshift) const {
// 	// pair-production mean free path on EBL (Broderick+2012 eq. 1), photon energy = 2 * electron energy
// 	// factor 0.5 in energy comes from the average energy of the parent
// 	double lambdaPP = 35. * crpropa::Mpc * (0.5 * crpropa::TeV / energy) * pow((1. + redshift) / 2., -4.5);

// 	// IC energy-loss rate in Thomson regime
// 	static const double u_CMB = 4.178e-14; 
// 	static const double mec2 = crpropa::mass_electron * crpropa::c_squared;
// 	double GammaIC = (4. / 3.) * crpropa::sigma_thomson * crpropa::c_light * u_CMB * (energy / mec2) * pow(1. + redshift, 4.) / mec2;

// 	// eq. 7 of Broderick et al. 2012
// 	return luminosity / (2. * M_PI * crpropa::pow_integer<3>(lambdaPP) * GammaIC) / energy;
// }



/*****************************************************************************/
/*                                  FlowJet1D                                */
/*****************************************************************************/

FlowJet1D::FlowJet1D() {
}

FlowJet1D::FlowJet1D(const std::vector<double>& distances, const std::vector<double>& beamDensity, double luminosity, crpropa::Vector3d centre, bool interpolateLog) {
	if (beamDensity.size() != distances.size()) {
		KISS_LOG_ERROR << "Vectors containing beam profile information should have the same size.";
		throw std::invalid_argument("Vectors containing beam profile information should have the same size.");
	}

	setOrigin(centre);
	setLuminosity(luminosity);
	setDistanceProfile(distances);
	setDensityProfile(beamDensity);
	setInterpolateLog(interpolateLog);
}

void FlowJet1D::setDistanceProfile(const std::vector<double>& distances) {
	distance = distances;
}

void FlowJet1D::setDensityProfile(const std::vector<double>& density) {
	densityProfile = density;
}

void FlowJet1D::setInterpolateLog(bool b) {
	interpolateLog = b;
}

std::vector<double> FlowJet1D::getDistanceProfile() const {
	return distance;
}

std::vector<double> FlowJet1D::getDensityProfile() const {
	return densityProfile;
}


double FlowJet1D::getDensity(double energy, const crpropa::Vector3d& position, double redshift) const {
	double x = (position - origin).getR();
	if (interpolateLog)
		x = log10(x);
	double n = crpropa::interpolate(x, distance, densityProfile);
	return n * crpropa::pow_integer<3>(1 + redshift);
}


/*****************************************************************************/
/*                           FlowMiniati2013                                 */
/*****************************************************************************/


crpropa::ref_ptr<Flow> createFlowMiniati2013(double luminosity, crpropa::Vector3d centre) {
	std::vector<double> distance = {0.87, 1.39, 2.22, 3.55, 5.68, 9.09, 14.55, 23.28, 37.25, 59.60, 95.37, 152.59, 244.14, 390.63, 625., 1000.}; // Mpc
	std::vector<double> density = {2.81e-18, 1.17e-18, 4.73e-19, 1.79e-19, 7.48e-20, 2.93e-20, 1.14e-20, 4.48e-21, 1.65e-21, 5.25e-22, 1.86e-22, 6.31e-23, 2.03e-23, 6.13e-24, 1.75e-24, 4.71e-25}; // cm^-3

	for (size_t i = 0; i < distance.size(); i++) {
		distance[i] *= crpropa::Mpc; // to m
		density[i] *= 1e7; // to m^-3
		density[i] *= (luminosity / 1e38); // scale with luminosity
	}
	
	crpropa::ref_ptr<FlowJet1D> flow = new FlowJet1D(distance, density, luminosity, centre, true);

	return flow;
}




} // namespace grplinst
