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


// double Flow::getPairProductionMeanFreePath(double E, double z) const {
// 	double rate = interpolate(E, tabEnergy, tabRate);
// 	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);
// }


/*****************************************************************************/
/*                              FlowHomogeneous                              */
/*****************************************************************************/

FlowHomogeneous::FlowHomogeneous(double L, crpropa::Vector3d centre) { 
	setOrigin(centre);
	setLuminosity(L);
}

FlowHomogeneous::FlowHomogeneous() {
}

double FlowHomogeneous::getDensity(double energy, const crpropa::Vector3d& position, double redshift) const {
	// pair-production mean free path on EBL (Broderick+2012 eq. 1), photon energy = 2 * electron energy
	// factor 0.5 in energy comes from the average energy of the parent
	double lambdaPP = 35. * crpropa::Mpc * (0.5 * crpropa::TeV / energy) * pow((1. + redshift) / 2., -4.5);

	// IC energy-loss rate in Thomson regime
	static const double u_CMB = 4.178e-14; 
	static const double mec2 = crpropa::mass_electron * crpropa::c_squared;
	double GammaIC = (4. / 3.) * crpropa::sigma_thomson * crpropa::c_light * u_CMB * (energy / mec2) * pow(1. + redshift, 4.) / mec2;

	// eq. 7 of Broderick et al. 2012
	return luminosity / (2. * M_PI * crpropa::pow_integer<3>(lambdaPP) * GammaIC) / energy;
}


double FlowHomogeneous::getMeanLorentzFactor(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	return lorentzFactorParticle * (1 + redshift);
}

double FlowHomogeneous::getMeanInverseLorentzFactor(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	return 1. / lorentzFactorParticle / (1 + redshift);
}

double FlowHomogeneous::getMeanLorentzFactorSquared(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	return pow(lorentzFactorParticle * (1 + redshift), 2);
}

double FlowHomogeneous::getAngularSpread(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	return 1. / (lorentzFactorParticle * (1 + redshift));
}

// double FlowHomogeneous::estimateBeamDensity(double E, double z) const {
// 	return 3.7e-16 * pow((1. + z) / 2., 9.5) * E * (luminosity / 1e38) * (E / crpropa::TeV);
// 	// double d = crpropa::redshift2LightTravelDistance(z);
// 	// double rIC = 1.2 * crpropa::pow_integer<3>(1 + z) * crpropa::kpc;
// 	// double L = luminosity;
// 	// return L / (2 * M_PI * crpropa::pow_integer<3>(d) * rIC) / E;
// }

/*****************************************************************************/
/*                                  FlowJet1D                                */
/*****************************************************************************/

FlowJet1D::FlowJet1D(const std::vector<double>& distances, const std::vector<double>& beamDensity, const std::vector<double>& lorentzFactor, const std::vector<double>& inverseLorentzFactor, const std::vector<double>& meanLorentzFactorSquared, const std::vector<double>& angularSpread, double densityNorm, crpropa::Vector3d centre, bool interpolateLog) {
	if ((beamDensity.size() != distances.size()) or (lorentzFactor.size() != distances.size()) or (inverseLorentzFactor.size() != distances.size()) or (meanLorentzFactorSquared.size() != distances.size()) or (angularSpread.size() != distances.size())) {
		throw std::length_error("Vectors containing beam profile information should have the same size.");
	}
	setOrigin(centre);
	setDensityNormalisation(densityNorm);
	setDistanceProfile(distances);
	setDensityProfile(beamDensity);
	setLorentzFactorProfile(lorentzFactor);
	setInverseLorentzFactorProfile(inverseLorentzFactor);
	setMeanLorentzFactorSquaredProfile(meanLorentzFactorSquared);
	setAngularSpreadProfile(angularSpread);
	setInterpolateLog(interpolateLog);
}

FlowJet1D::FlowJet1D() {
}

void FlowJet1D::setDensityNormalisation(double n) {
	densityNormalisation = n;
}

double FlowJet1D::getDensityNormalisation() const {
	return densityNormalisation;
}

void FlowJet1D::setDistanceProfile(const std::vector<double>& distances) {
	for (size_t i = 0; i < distances.size(); i++) {
		distance.push_back(distances[i]);
	}
}

void FlowJet1D::setLorentzFactorProfile(const std::vector<double>& lorentzFactor) {
	for (size_t i = 0; i < lorentzFactor.size(); i++) {
		meanLorentzFactor.push_back(lorentzFactor[i]);
	}
}

void FlowJet1D::setInverseLorentzFactorProfile(const std::vector<double>& inverseLorentzFactor) {
	for (size_t i = 0; i < inverseLorentzFactor.size(); i++) {
		meanInverseLorentzFactor.push_back(inverseLorentzFactor[i]);
	}
}

void FlowJet1D::setDensityProfile(const std::vector<double>& density) {
	for (size_t i = 0; i < density.size(); i++) {
		densityProfile.push_back(density[i] * densityNormalisation);
	}
}

void FlowJet1D::setMeanLorentzFactorSquaredProfile(const std::vector<double>& lfSquared) {
	for (size_t i = 0; i < lfSquared.size(); i++) {
		meanLorentzFactorSquared.push_back(lfSquared[i]);
	}
}

void FlowJet1D::setAngularSpreadProfile(const std::vector<double>& as) {
	for (size_t i = 0; i < as.size(); i++) {
		angularSpread.push_back(as[i]);
	}
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

std::vector<double> FlowJet1D::getLorentzFactorProfile() const {
	return meanLorentzFactor;
}

std::vector<double> FlowJet1D::getInverseLorentzFactorProfile() const {
	return meanInverseLorentzFactor;
}

std::vector<double> FlowJet1D::getMeanLorentzFactorSquaredProfile() const {
	return meanLorentzFactorSquared;
}

std::vector<double> FlowJet1D::getAngularSpreadProfile() const {
	return angularSpread;
}

double FlowJet1D::getDensity(double energy, const crpropa::Vector3d& position, double redshift) const {
	double x = (position - origin).getR();
	if (interpolateLog)
		x = log10(x);
	double n = crpropa::interpolate(x, distance, densityProfile);
	return n * crpropa::pow_integer<3>(1 + redshift);
}

double FlowJet1D::getMeanLorentzFactor(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	double x = (position - origin).getR();
	if (interpolateLog)
		x = log10(x);
	double lf = crpropa::interpolate(x, distance, meanLorentzFactor);
	return lf * (1 + redshift);
}

double FlowJet1D::getMeanInverseLorentzFactor(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	double x = (position - origin).getR();
	if (interpolateLog)
		x = log10(x);
	double ilf = crpropa::interpolate(x, distance, meanInverseLorentzFactor);
	return ilf / (1 + redshift);
}

double FlowJet1D::getMeanLorentzFactorSquared(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	double x = (position - origin).getR();
	if (interpolateLog)
		x = log10(x);
	double lf2 = crpropa::interpolate(x, distance, meanLorentzFactorSquared);
	return lf2 * pow(1 + redshift, 2);
}

double FlowJet1D::getAngularSpread(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	double x = (position - origin).getR();
	if (interpolateLog)
		x = log10(x);
	double as = crpropa::interpolate(x, distance, angularSpread);
	return as;
}




} // namespace grplinst
