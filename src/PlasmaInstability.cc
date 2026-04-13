#include "grplinst/PlasmaInstability.h"


namespace grplinst {





/*****************************************************************************/
/*                             PlasmaInstability                             */
/*****************************************************************************/


PlasmaInstability::PlasmaInstability(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) {
	setFlowProperties(flow);
	setMediumDensity(density);
	setMediumTemperature(temperature);
	setLimit(limit);
	setDescription("PlasmaInstability::PlasmaInstability");
}


void PlasmaInstability::setFlowProperties(crpropa::ref_ptr<Flow> flow) {
	flowProperties = flow;
}

void PlasmaInstability::setMediumDensity(crpropa::ref_ptr<MediumDensity> density) {
	mediumDensity = density;
}

void PlasmaInstability::setMediumTemperature(crpropa::ref_ptr<MediumTemperature> temperature) {
	mediumTemperature = temperature;
}

void PlasmaInstability::setLimit(double l) {
	limit = l;
}

crpropa::ref_ptr<Flow> PlasmaInstability::getFlowProperties() const {
	return flowProperties;
}

crpropa::ref_ptr<MediumDensity> PlasmaInstability::getMediumDensity() const {
	return mediumDensity;
}

crpropa::ref_ptr<MediumTemperature> PlasmaInstability::getMediumTemperature() const {
	return mediumTemperature;
}

void PlasmaInstability::process(crpropa::Candidate* candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double dx = candidate->getCurrentStep() / (1 + z);

	double dEdx = computeEnergyLossPerLength(candidate);
	if (dEdx <= 0) // prevent overshooting
		return;
	
	double Enew = std::max(0.0, E - dEdx * dx);
	candidate->current.setEnergy(Enew / (1 + z));
	candidate->limitNextStep(limit * E / dEdx);
}


double PlasmaInstability::computeEnergyLossPerLength(crpropa::Candidate* candidate) const {
	double tau = energyLossTime(candidate);
	if (tau <= 0.)
		return 0.;

	return candidate->current.getEnergy() / (crpropa::c_light * tau);
}


/*****************************************************************************/
/*                    Implementation for each model                          */
/*****************************************************************************/

double PlasmaInstabilityBroderick2012::energyLossTime(const crpropa::Candidate& candidate) const {
	crpropa::Vector3d position = candidate.current.getPosition();
	double z = candidate.getRedshift();
	double E = candidate.current.getEnergy() * (1 + z);
	
	double n = mediumDensity->getDensity(position, z);
	double nb = flowProperties->getDensity(E, position, z);

	double nCrit = 1.6e-13 / crpropa::pow_integer<2>(E / crpropa::TeV) * (n / 0.1);
	if (nb < nCrit) {
		return 7e7 / (E / crpropa::TeV) / (nb / 1e-16) * sqrt(n / 0.1);
	} else {
		return 5e5 * cbrt(E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.);
	}
}

double PlasmaInstabilitySchlickeiser2012::energyLossTime(const crpropa::Candidate& candidate) const {
	crpropa::Vector3d position = candidate.current.getPosition();
	double z = candidate.getRedshift();
	double E = candidate.current.getEnergy() * (1 + z);

	double n = mediumDensity->getDensity(position, z);
	double nb = flowProperties->getDensity(E, position, z);
	double T = mediumTemperature->getTemperature(position, z);

	double nCrit = 2.5e-19 / (E / crpropa::TeV) * (n / 0.1) * crpropa::pow_integer<2>(T / 1e4);
	if (nb < nCrit) {
		return 5e14 * pow(E / crpropa::TeV, 5. / 3.) * cbrt(nb / 1e-16) * pow(n / 0.1, -5. / 6.) / crpropa::pow_integer<2>(T / 1e4);
	} else {
		return 8e6 * cbrt(E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.) * (1 + 1.25 * log(T / 1e4) - 0.25 * log(n / 0.1));
	}
}

double PlasmaInstabilitySironi2014::energyLossTime(const crpropa::Candidate& candidate) const {
	crpropa::Vector3d position = candidate.current.getPosition();
	double z = candidate.getRedshift();
	double E = candidate.current.getEnergy() * (1 + z);

	double n = mediumDensity->getDensity(position, z);
	double nb = flowProperties->getDensity(E, position, z);

	double nCrit = 8e-14 / crpropa::pow_integer<2>(E / crpropa::TeV) * (n / 0.1);
	if (nb < nCrit) {
		return 1.4e7 / (E / crpropa::TeV) / (nb / 1e-16) * sqrt(n / 0.1);
	} else {
		return 9.6e5 * cbrt(E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.);
	}
}

double PlasmaInstabilityVafin2018::energyLossTime(const crpropa::Candidate& candidate) const {
	crpropa::Vector3d position = candidate.current.getPosition();
	double z = candidate.getRedshift();
	double E = candidate.current.getEnergy() * (1 + z);

	double n = mediumDensity->getDensity(position, z);
	double nb = flowProperties->getDensity(E, position, z);
	double T = mediumTemperature->getTemperature(position, z);

	 return 1.9e11 * pow(E / crpropa::TeV, 4. / 3.) * cbrt(n / 0.1) / cbrt(nb / 1e-16) / (T / 1e4);
}

double PlasmaInstabilityBret2010Filamentation::energyLossTime(const crpropa::Candidate& candidate) const {
	crpropa::Vector3d position = candidate.current.getPosition();
	double z = candidate.getRedshift();
	double E = candidate.current.getEnergy() * (1 + z);

	double n = mediumDensity->getDensity(position, z);
	double nb = flowProperties->getDensity(E, position, z);
	double T = mediumTemperature->getTemperature(position, z);

	return 2.5e9 * sqrt(E / crpropa::TeV) / sqrt(nb / 1e-16);
}

double PlasmaInstabilityBret2010TwoStream::energyLossTime(const crpropa::Candidate& candidate) const {	 		
	crpropa::Vector3d position = candidate.current.getPosition();
	double z = candidate.getRedshift();
	double E = candidate.current.getEnergy() * (1 + z);

	double n = mediumDensity->getDensity(position, z);
	double nb = flowProperties->getDensity(E, position, z);
	double T = mediumTemperature->getTemperature(position, z);

	return 1.6e10 * (E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.);
}

double PlasmaInstabilityShalaby2020::energyLossTime(const crpropa::Candidate& candidate) const {
	crpropa::Vector3d position = candidate.current.getPosition();
	double z = candidate.getRedshift();
	double E = candidate.current.getEnergy() * (1 + z);

	double n = mediumDensity->getDensity(position, z);
	double nb = flowProperties->getDensity(E, position, z);
	double T = mediumTemperature->getTemperature(position, z);

	return 3.2e12 * pow(E / crpropa::TeV, 6. / 5.) * pow(nb / 1e-16, -2. / 5.) * pow(n / 0.1, -1. / 10.);
}

PlasmaInstabilityMiniati2013::PlasmaInstabilityMiniati2013() : PlasmaInstability() {
	initFlow();
}

void PlasmaInstabilityMiniati2013::initFlow() {
	std::ifstream inputFile(filename.c_str());

	if (! inputFile.good())
		throw std::runtime_error("PlasmaInstabilityMiniati2013: could not open file " + filename + ".");

	std::vector<double> d, n, g, g_1, g2, dQ;

	while (inputFile.good()) {
		if (inputFile.peek() != '#') {
			double _d, _n, _g, _g_1, _g2, _dQ;
			inputFile >> _d >> _n >> _g >> _g_1 >> _g2 >> _dQ;
			
			if (inputFile.good()) {
				_d = log10(_d * crpropa::Mpc);
				_n = _n / crpropa::ccm;
				_g = _g * 1e5;
				_g_1 = _g_1 * 1e4;
				_g2 = _g2 * 1e6 * _g;
				_dQ = _dQ * 1e-5;
				d.push_back(_d);
				n.push_back(_n);
				g.push_back(_g);
				g_1.push_back(_g_1);
				g2.push_back(_g2);
				dQ.push_back(_dQ);

			}
		}
		inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	inputFile.close();

	flow = new FlowJet1D(d, n, g, g_1, g2, dQ);
}


double PlasmaInstabilityMiniati2013::energyLossTime(const crpropa::Candidate& candidate) const {
	crpropa::Vector3d position = candidate.current.getPosition();
	double z = candidate.getRedshift();
	double E = candidate.current.getEnergy() * (1 + z);

	double n = mediumDensity->getDensity(position, z);
	double nb = flowProperties->getDensity(E, position, z);
	double T = mediumTemperature->getTemperature(position, z);
	double lf = flow->getMeanLorentzFactor(position, z);
	double ilf = flow->getMeanInverseLorentzFactor(position, z);
	double dTh = flow->getAngularSpread(position, z);
	return 1.5e9 * crpropa::year / (n / 2e-8 / crpropa::ccm) * (ilf / 1e-4) * (lf / 1e5) * crpropa::pow_integer<2>(dTh / 1e-4) * (T / 3e3);
}


/*****************************************************************************/
/*                         Additional Functions                              */
/*****************************************************************************/

double maximumLinearGrowthFrequency(double beamDensity, double mediumDensity, double inverseLorentzFactor, int id) {
	double wp =  plasmaFrequency(mediumDensity, id);
	return wp * beamDensity / mediumDensity * inverseLorentzFactor;
}

double plasmaFrequency(double density, int id) {
	// charge
	double q = crpropa::eplus;
	if (crpropa::isNucleus(id)) {
		q *= crpropa::chargeNumber(id);
	}

	// mass
	double m = 0;
	if (fabs(id) == 11) {
		m = crpropa::mass_electron;
	} else if (crpropa::isNucleus(id)) {
		m = crpropa::nuclearMass(id);
	} else {
		std::cerr << "Mass undefined for particle with id: " << id << std::endl;
	}

	return sqrt(density * q * q / m / crpropa::epsilon0);
}


} // namespace grplinst