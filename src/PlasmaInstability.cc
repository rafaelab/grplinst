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
	if (dEdx < 0) // prevent overshooting
		dEdx = 0;

	double Enew = E - dEdx * dx;

	candidate->current.setEnergy(Enew / (1 + z));
	candidate->limitNextStep(limit * E / dEdx);
}

double PlasmaInstability::computeEnergyLossPerLength(crpropa::Candidate* candidate) const {
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double dx = candidate->getCurrentStep() / (1 + z);

	double n = mediumDensity->getDensity(candidate->current.getPosition(), z);
	double T = mediumTemperature->getTemperature(candidate->current.getPosition(), z);
	double nb = flowProperties->getDensity(candidate->current.getPosition(), z);
	
	double tau = energyLossTime(E, nb, n, T);


	return E / (crpropa::c_light * tau);
}



/*****************************************************************************/
/*                    Implementation for each model                          */
/*****************************************************************************/

double PlasmaInstabilityBroderick2012::energyLossTime(double E, double nb, double n, double T) const {
	double nCrit = 1.6e-13 / crpropa::pow_integer<2>(E / crpropa::TeV) * (n / 0.1);
	if (nb < nCrit) {
		return 7e7 / (E / crpropa::TeV) / (nb / 1e-16) * sqrt(n / 0.1);
	} else {
		return 5e5 * cbrt(E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.);
	}
}

double PlasmaInstabilitySchlickeiser2012::energyLossTime(double E, double nb, double n, double T) const {
	double nCrit = 2.5e-19 / (E / crpropa::TeV) * (n / 0.1) * crpropa::pow_integer<2>(T / 1e4);
	if (nb < nCrit) {
		return 5e14 * pow(E / crpropa::TeV, 5. / 3.) * cbrt(nb / 1e-16) * pow(n / 0.1, -5. / 6.) * crpropa::pow_integer<2>(T / 1e4);
	} else {
		return 8e6 * cbrt(E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.) * (1 + 1.25 * log(T / 1e4) - 0.25 * log(n / 0.1));
	}
}

double PlasmaInstabilitySironi2014::energyLossTime(double E, double nb, double n, double T) const {
	double nCrit = 8e-14 / crpropa::pow_integer<2>(E / crpropa::TeV) * (n / 0.1);
	if (nb < nCrit) {
		return 1.4e7 / (E / crpropa::TeV) / (nb / 1e-16) * sqrt(n / 0.1);
	} else {
		return 9.6e5 * cbrt(E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.);
	}
}

double PlasmaInstabilityVafin2018::energyLossTime(double E, double nb, double n, double T) const {
	return 1.9e11 * pow(E / crpropa::TeV, 4. / 3.) / cbrt(nb / 1e-16) * cbrt(n / 0.1) / (T / 1e4);
}

double PlasmaInstabilityBret2010Filamentation::energyLossTime(double E, double nb, double n, double T) const {
	return 2.5e9 * sqrt(E / crpropa::TeV) / sqrt(nb / 1e-16);
}

double PlasmaInstabilityBret2010TwoStream::energyLossTime(double E, double nb, double n, double T) const {
	return 1.6e10 * (E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.);
}

double PlasmaInstabilityShalaby2020::energyLossTime(double E, double nb, double n, double T) const {
	return 3.2e12 * pow(E / crpropa::TeV, 6. / 5.) * pow(nb / 1e-16, -2. / 5.) * pow(n / 0.1, -1. / 10.);
}



/*****************************************************************************/
/*                         Additional Functions                              */
/*****************************************************************************/

double maximumLinearGrowthFrequency(double beamDensity, double mediumDensity, double inverseLorentzFactor, int id) {
	double wp =  plasmaFrequency(mediumDensity, id);
	return wp * beamDensity / mediumDensity / inverseLorentzFactor;
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