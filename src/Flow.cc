#include "grplinst/Flow.h"


namespace grplinst {



/***************************************************************************/
/**/
Flow::Flow(crpropa::Vector3d origin) {
}

Flow::~Flow() {
}

void Flow::setOrigin(const crpropa::Vector3d &centre) {
	origin = centre;
}

crpropa::Vector3d Flow::getOrigin() const {
	return origin;
}


/***************************************************************************/
/**/
FlowHomogeneous::FlowHomogeneous(double n, crpropa::Vector3d centre) {
	setOrigin(centre);
	setTotalDensity(n);
}

FlowHomogeneous::FlowHomogeneous() {
}

FlowHomogeneous::~FlowHomogeneous() {
}

void FlowHomogeneous::setTotalDensity(double n) {
	density = n;
}

double FlowHomogeneous::getTotalDensity() const {
	return density;
}

double FlowHomogeneous::getDensity(const crpropa::Vector3d &position, double redshift) const {
	return density * crpropa::pow_integer<3>(1 + redshift);
}

double FlowHomogeneous::getMeanLorentzFactor(const crpropa::Vector3d &position, double redshift, double lorentzFactorParticle) const {
	return lorentzFactorParticle * (1 + redshift);
}

double FlowHomogeneous::getMeanInverseLorentzFactor(const crpropa::Vector3d &position, double redshift, double lorentzFactorParticle) const {
	return 1. / lorentzFactorParticle / (1 + redshift);
}

/***************************************************************************/
/**/
FlowJet1D::FlowJet1D(const std::vector<double> &distances, const std::vector<double> &beamDensity, const std::vector<double> &lorentzFactor, const std::vector<double> &inverseLorentzFactor, double densityNorm, crpropa::Vector3d centre) {
	if ((beamDensity.size() != distances.size()) || (lorentzFactor.size() != distances.size()) || (inverseLorentzFactor.size() != distances.size())) {
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
	if (!infile.good()) {
		throw std::runtime_error("FlowJet1D could not open file " + filename + ".");
	}
	std::string line;
	while (std::getline(infile,line)) {
		std::stringstream stream(line);
		if (stream.peek() == '#')
			continue;

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

FlowJet1D::~FlowJet1D() {
}

void FlowJet1D::setDensityNormalisation(double n) {
	densityNormalisation = n;
}

double FlowJet1D::getDensityNormalisation() const {
	return densityNormalisation;
}

void FlowJet1D::setDistanceProfile(const std::vector<double> &distances) {
	for (size_t i = 0; i < distances.size(); i++) {
		distance.push_back(distances[i]);
	}
}

void FlowJet1D::setLorentzFactorProfile(const std::vector<double> &lorentzFactor) {
	for (size_t i = 0; i < lorentzFactor.size(); i++) {
		meanLorentzFactor.push_back(lorentzFactor[i]);
	}
}

void FlowJet1D::setInverseLorentzFactorProfile(const std::vector<double> &inverseLorentzFactor) {
	for (size_t i = 0; i < inverseLorentzFactor.size(); i++) {
		meanInverseLorentzFactor.push_back(inverseLorentzFactor[i]);
	}
}

void FlowJet1D::setDensityProfile(const std::vector<double> &density) {
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

double FlowJet1D::getDensity(const crpropa::Vector3d &position, double redshift) const {
	double n = crpropa::interpolate((position - origin).getR(), distance, densityProfile);
	return n * crpropa::pow_integer<3>(1 + redshift);
}

double FlowJet1D::getMeanLorentzFactor(const crpropa::Vector3d &position, double redshift, double lorentzFactorParticle) const {
	double lf = crpropa::interpolate((position - origin).getR(), distance, meanLorentzFactor);
	return lf * (1 + redshift);
}

double FlowJet1D::getMeanInverseLorentzFactor(const crpropa::Vector3d &position, double redshift, double lorentzFactorParticle) const {
	double ilf = crpropa::interpolate((position - origin).getR(), distance, meanInverseLorentzFactor);
	return ilf / (1 + redshift);
}


/***************************************************************************/
/**/
FlowJet1DMahmoud::FlowJet1DMahmoud(double densityNorm, crpropa::Vector3d centre) {

	// the first element (0) is going to be inferred from extrapolation
	std::vector<double> d = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
	std::vector<double> n = {0, 2.65e-21, 4.17e-22, 1.19e-22, 4.36e-23, 1.87e-23, 8.94e-24, 4.67e-24, 2.63e-24, 1.58e-24};
	// std::vector<double> lf = {2.1e5, 2.013800e+05, 1.188416e+05, 7.559221e+04, 5.522591e+04, 4.629579e+04, 4.120994e+04, 3.600084e+04, 3.187101e+04, 2.766477e+04, 2.553351e+04, 2.382528e+04, 2.195751e+04, 2.093016e+04, 2.107416e+04, 1.796684e+04, 1.675042e+04, 1.736836e+04, 1.634214e+04, 1.485561e+04};
	// std::vector<double> ilf = {1.066650e-05, 2.045386e-05, 2.976448e-05, 3.951229e-05, 4.828246e-05, 5.684162e-05, 6.584768e-05, 7.508342e-05, 8.467142e-05, 9.500127e-05, 1.039197e-04, 1.133797e-04, 1.237442e-04, 1.328346e-04, 1.426141e-04, 1.527329e-04, 1.606893e-04, 1.696674e-04, 1.782544e-04};
	std::vector<double> lf = {2.013800e+05, 1.188416e+05, 5.522591e+04, 4.120994e+04, 3.187101e+04, 2.553351e+04, 2.195751e+04, 2.107416e+04,  1.675042e+04, 1.634214e+04, 1.485561e+04};
	std::vector<double> ilf = {1.066650e-05, 2.976448e-05, 4.828246e-05, 6.584768e-05, 8.467142e-05, 1.039197e-04, 1.237442e-04, 1.426141e-04, 1.606893e-04, 1.782544e-04};

	// // simple extrapolation
	// w_[0] = exp((log(w_[2] / w_[1]) / (d_[2] - d_[1])) * (d_[1] - d_[0]) - log(w_[1]));
	n[0] = n[1] - ((n[2] - n[1]) / (d[2] - d[1])) * (d[1] - d[0]);

	for (size_t i = 0; i < distance.size(); i++) {
		d[i] *= crpropa::Mpc;
		n[i] /= n[0];
	}

	setDensityNormalisation(densityNorm);
	setDensityProfile(n);
	setDistanceProfile(d);
	setLorentzFactorProfile(lf);
	setInverseLorentzFactorProfile(ilf);
	setOrigin(centre);
}

FlowJet1DMahmoud::~FlowJet1DMahmoud() {



}

} // namespace grplinst