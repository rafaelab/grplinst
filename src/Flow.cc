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

FlowHomogeneous::FlowHomogeneous(double n, crpropa::Vector3d centre) { 
	setOrigin(centre);
	setDensityValue(n);
}

FlowHomogeneous::FlowHomogeneous() {
}

void FlowHomogeneous::setDensityValue(double n) {
	density = n;
}

double FlowHomogeneous::getDensityValue() const {
	return density;
}

double FlowHomogeneous::getDensity(const crpropa::Vector3d& position, double redshift) const {
	return density * crpropa::pow_integer<3>(1 + redshift);
}

double FlowHomogeneous::getMeanLorentzFactor(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	return lorentzFactorParticle * (1 + redshift);
}

double FlowHomogeneous::getMeanInverseLorentzFactor(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	return 1. / lorentzFactorParticle / (1 + redshift);
}

double FlowHomogeneous::estimateBeamDensity(double E, double z) const {
	return 3.7e-16 * pow((1. + z) / 2., 9.5) * E * (luminosity / 1e38) * (E / crpropa::TeV);
}

/*****************************************************************************/
/*                                  FlowJet1D                                */
/*****************************************************************************/

FlowJet1D::FlowJet1D(const std::vector<double>& distances, const std::vector<double>& beamDensity, const std::vector<double>& lorentzFactor, const std::vector<double>& inverseLorentzFactor, double densityNorm, crpropa::Vector3d centre) {
	if ((beamDensity.size() != distances.size()) or (lorentzFactor.size() != distances.size()) or (inverseLorentzFactor.size() != distances.size())) {
		std::length_error("Vectors containing beam profile information should have the same size.");
	}
	setOrigin(centre);
	setDensityNormalisation(densityNorm);
	setDistanceProfile(distances);
	setDensityProfile(beamDensity);
	setLorentzFactorProfile(lorentzFactor);
	setInverseLorentzFactorProfile(inverseLorentzFactor);
	
}

FlowJet1D::FlowJet1D(const std::string &filename, double densityNorm, crpropa::Vector3d centre) {
	setOrigin(centre);
	setDensityNormalisation(densityNorm);

	// read file and store jet profile
	std::ifstream infile(filename.c_str());
	if (! infile.good()) {
		throw std::runtime_error("FlowJet1D could not open file " + filename + ".");
	}
	std::string line;
	while (std::getline(infile, line)) {
		auto firstNon = line.find_first_not_of(" \t");
		if (firstNon == std::string::npos or line[firstNon] == '#') 
			continue;
		std::stringstream stream(line);

		double d;
		double n;
		double gamma;
		double gamma_1;
		stream >> d >> n >> gamma >> gamma_1;

		n *= densityNorm;
	
		distance.push_back(d);
		densityProfile.push_back(n);
		meanLorentzFactor.push_back(gamma);
		meanInverseLorentzFactor.push_back(gamma_1);
	}
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

double FlowJet1D::getDensity(const crpropa::Vector3d& position, double redshift) const {
	double n = crpropa::interpolate((position - origin).getR(), distance, densityProfile);
	return n * crpropa::pow_integer<3>(1 + redshift);
}

double FlowJet1D::getMeanLorentzFactor(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	double lf = crpropa::interpolate((position - origin).getR(), distance, meanLorentzFactor);
	return lf * (1 + redshift);
}

double FlowJet1D::getMeanInverseLorentzFactor(const crpropa::Vector3d& position, double redshift, double lorentzFactorParticle) const {
	double ilf = crpropa::interpolate((position - origin).getR(), distance, meanInverseLorentzFactor);
	return ilf / (1 + redshift);
}





} // namespace grplinst