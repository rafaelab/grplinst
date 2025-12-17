#include "grplinst/PlasmaInstability.h"


namespace grplinst {


/*****************************************************************************/
/*                             PlasmaInstability                             */
/*****************************************************************************/


PlasmaInstability::PlasmaInstability(PlasmaInstabilityModel model, crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) {
	setModel(model);
	setFlowProperties(flow);
	setMediumDensity(density);
	setMediumTemperature(temperature);
	setLimit(limit);
	setDescription("PlasmaInstability::PlasmaInstability");
}

void PlasmaInstability::setModel(PlasmaInstabilityModel m) {
	model = m;
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

PlasmaInstabilityModel PlasmaInstability::getModel() const {
	return model;
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

double PlasmaInstability::computeEnergyLossTime(double E, double nb, double n, double T) const {
	switch(model) {
		case PlasmaInstabilityModel::Broderick2012:
			return computeEnergyLossTime_Broderick2012(E, nb, n, T);
		// case PlasmaInstabilityModel::Miniati2013:
		// 	break;
		case PlasmaInstabilityModel::Schlickeiser2012:
			return computeEnergyLossTime_Schlickeiser2012(E, nb, n, T);
		case PlasmaInstabilityModel::Sironi2014:
			return computeEnergyLossTime_Sironi2014(E, nb, n, T);
		case PlasmaInstabilityModel::Vafin2018:
			return computeEnergyLossTime_Vafin2018(E, nb, n, T);
		case PlasmaInstabilityModel::Bret2010_2s:
			return computeEnergyLossTime_Bret2010_2s(E, nb, n, T);
		case PlasmaInstabilityModel::Bret2010_f:
			return computeEnergyLossTime_Bret2010_f(E, nb, n, T);
		case PlasmaInstabilityModel::Shalaby2020:
			return computeEnergyLossTime_Shalaby2020(E, nb, n, T);
		default:
			std::cerr << "Error: PlasmaInstability model not recognized." << std::endl;
			return std::numeric_limits<double>::infinity();
	}
}

double PlasmaInstability::computeEnergyLossPerLength(crpropa::Candidate* candidate) const {
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double dx = candidate->getCurrentStep() / (1 + z);

	double n = mediumDensity->getDensity(candidate->current.getPosition(), z);
	double T = mediumTemperature->getTemperature(candidate->current.getPosition(), z);
	double nb = flowProperties->getDensity(candidate->current.getPosition(), z);
	
	double tau = computeEnergyLossTime(E, nb, n, T);


	return E / (crpropa::c_light * tau);
}



/*****************************************************************************/
/*                    Implementation for each model                          */
/*****************************************************************************/


double PlasmaInstability::computeEnergyLossTime_Broderick2012(double E, double nb, double n, double T) {
	double nCrit = 1.6e-13 / crpropa::pow_integer<2>(E / crpropa::TeV) * (n / 0.1);
	if (nb < nCrit) {
		return 7e7 / (E / crpropa::TeV) / (nb / 1e-16) * sqrt(n / 0.1);
	} else {
		return 5e5 * cbrt(E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.);
	}
}

double PlasmaInstability::computeEnergyLossTime_Schlickeiser2012(double E, double nb, double n, double T) {
	double nCrit = 2.5e-19 / (E / crpropa::TeV) * (n / 0.1) * crpropa::pow_integer<2>(T / 1e4);
	if (nb < nCrit) {
		return 5e14 * pow(E / crpropa::TeV, 5. / 3.) * cbrt(nb / 1e-16) * pow(n / 0.1, -5. / 6.) * crpropa::pow_integer<2>(T / 1e4);
	} else {
		return 8e6 * cbrt(E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.) * (1 + 1.25 * log(T / 1e4) - 0.25 * log(n / 0.1));
	}
}

double PlasmaInstability::computeEnergyLossTime_Sironi2014(double E, double nb, double n, double T) {
	double nCrit = 8e-14 / crpropa::pow_integer<2>(E / crpropa::TeV) * (n / 0.1);
	if (nb < nCrit) {
		return 1.4e7 / (E / crpropa::TeV) / (nb / 1e-16) * sqrt(n / 0.1);
	} else {
		return 9.6e5 * cbrt(E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.);
	}
}

double PlasmaInstability::computeEnergyLossTime_Vafin2018(double E, double nb, double n, double T) {
	return 1.9e11 * pow(E / crpropa::TeV, 4. / 3.) / cbrt(nb / 1e-16) * cbrt(n / 0.1) / (T / 1e4);
}

double PlasmaInstability::computeEnergyLossTime_Bret2010_f(double E, double nb, double n, double T) {
	return 2.5e9 * sqrt(E / crpropa::TeV) / sqrt(nb / 1e-16);
}

double PlasmaInstability::computeEnergyLossTime_Bret2010_2s(double E, double nb, double n, double T) {
	return 1.6e9 * (E / crpropa::TeV) / cbrt(nb / 1e-16) / pow(n / 0.1, 1. / 6.);
}

double PlasmaInstability::computeEnergyLossTime_Shalaby2020(double E, double nb, double n, double T) {
	return 3.2e12 * pow(E / crpropa::TeV, 6. / 5.) * pow(nb / 1e-16, -2. / 5.) * pow(n / 0.1, -1. / 10.);
}

// double PlasmaInstability::computeEnergyLossTime_Miniati2013(double E, double nb, double n, double T) {
// 	double d0 = log10(crpropa::redshift2LightTravelDistance(z) / crpropa::Mpc); // co-moving?
// 	double w = pow(10, crpropa::interpolate(d0, _d, _w));
// 	return 1.4e-29 * pow(1 + z, 2) * pow(E / crpropa::TeV, 2) / w;
// }

// void initTable_Miniati2013() {
// 	double d_[] = {-0.05, 0.14, 0.35, 0.59, 0.79, 0.96, 1.17, 1.40, 1.57, 1.77, 1.99, 2.20, 2.41, 2.60, 2.80, 3.00};
// 	double w_[] = { 3.28, 3.46, 3.14, 3.08, 3.05, 3.00, 2.84, 2.56, 2.40, 2.55, 2.36, 2.08, 1.73, 1.55, 1.05, 0.75};
// 	std::vector<double> D (d_, d_ + sizeof(d_) / sizeof(double));
// 	std::vector<double> W (w_, w_ + sizeof(w_) / sizeof(double));
// }

// 	_d = D;
// 	_w = W;

// /***************************************************************************/
// /**/
// PlasmaInstabilityMiniati2013::PlasmaInstabilityMiniati2013(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
// 	initTable();
// 	setDescription("PlasmaInstability::PlasmaInstabilityMiniati2013");
// }

// PlasmaInstabilityMiniati2013B::PlasmaInstabilityMiniati2013B(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
// 	setDescription("PlasmaInstability::PlasmaInstabilityMiniati2013B");
// }

// PlasmaInstabilityMiniati2013C::PlasmaInstabilityMiniati2013C(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
// 	setDescription("PlasmaInstability::PlasmaInstabilityMiniati2013C");
// }

// void PlasmaInstabilityMiniati2013::initTable() {
// 	double d_[] = {-0.05, 0.14, 0.35, 0.59, 0.79, 0.96, 1.17, 1.40, 1.57, 1.77, 1.99, 2.20, 2.41, 2.60, 2.80, 3.00};
// 	double w_[] = { 3.28, 3.46, 3.14, 3.08, 3.05, 3.00, 2.84, 2.56, 2.40, 2.55, 2.36, 2.08, 1.73, 1.55, 1.05, 0.75};
// 	std::vector<double> D (d_, d_ + sizeof(d_) / sizeof(double));
// 	std::vector<double> W (w_, w_ + sizeof(w_) / sizeof(double));
// 	_d = D;
// 	_w = W;
// }

// double PlasmaInstabilityMiniati2013::energyLoss(crpropa::Candidate* candidate) const {
	
// 	double z = candidate->getRedshift();
// 	double E = candidate->current.getEnergy(); // multiply by (1 + z) for E(z)
	
// 	crpropa::Vector3d pos = candidate->current.getPosition();
// 	crpropa::Vector3d pos0 = candidate->source.getPosition();

// 	// double nMedium = mediumDensity->getValue(pos);
// 	// double nBeam = flowProperties->getValue(pos);
// 	// double L = beamDensityToLuminosity(nBeam, E, z);

// 	// medium and flow properties useless? 
// 	// recheck paper

// 	E /= (1 + z);
// 	double d0 = log10(crpropa::redshift2LightTravelDistance(z) / crpropa::Mpc); // co-moving?
// 	double w = pow(10, crpropa::interpolate(d0, _d, _w));

// 	return 1.4e-29 * crpropa::pow_integer<2>(1 + z) * crpropa::pow_integer<2>(E / crpropa::TeV) / w;
// }

// double PlasmaInstabilityMiniati2013B::energyLoss(crpropa::Candidate* candidate) const {
// 	int id = candidate->current.getId();
// 	if (fabs(id) != 11)
// 		return 0;

// 	double z = candidate->getRedshift();
// 	double E = candidate->current.getEnergy() * (1 + z);
// 	double lf = candidate->current.getLorentzFactor() * (1 + z);
// 	crpropa::Vector3d pos = candidate->current.getPosition();
// 	crpropa::Vector3d pos0 = candidate->source.getPosition();
// 	crpropa::Vector3d dir = candidate->current.getDirection();
// 	crpropa::Vector3d dir0 = candidate->source.getDirection();
	
// 	double nMedium = mediumDensity->getDensity(pos, z);
// 	double nBeam = flowProperties->getDensity(pos, z);
// 	double L = beamDensityToLuminosity(nBeam, E, z);
// 	double T = mediumTemperature->getTemperature(pos, z);
// 	double ve = mediumTemperature->getVelocity(id, pos, z);
// 	double vi = mediumTemperature->getVelocity(crpropa::nucleusId(1, 1), pos, z); 	// consider only protons

// 	double dTheta = (dir - dir0).getR();
// 	if (dTheta == 0)
// 		dTheta = 1e-3;

// 	double wp = plasmaFrequency(nMedium, id);
// 	double gmax = wp * (nBeam / nMedium) * 1 / (dTheta * dTheta) * lf; // Lorentz factor should be in denominator; this results in the right slope, though
// 	// double gmax = wp * (nBeam / nMedium) * 1 / (dTheta * dTheta * lf);
// 	double gnl = wp / (nMedium * T * crpropa::k_boltzmann) * (ve * ve / vi / crpropa::c_light);
// 	double nuc = 1e-11 * (nMedium / 0.02) * pow(T / 3e3, -1.5);
// 	double qr = nuc / 2. / gnl;
// 	double qnr = gmax / gnl;
// 	double tau = nBeam / 2 / gmax / qr;

// 	return E / tau / crpropa::c_light;
// }

// double PlasmaInstabilityMiniati2013C::energyLoss(crpropa::Candidate* candidate) const {
// 	int id = candidate->current.getId();
// 	if (fabs(id) != 11)
// 		return 0;

// 	double z = candidate->getRedshift();
// 	double E = candidate->current.getEnergy() * (1 + z); // multiply by (1 + z) for E(z)
// 	double lf = candidate->current.getLorentzFactor() * (1 + z);
// 	crpropa::Vector3d pos = candidate->current.getPosition();
// 	crpropa::Vector3d pos0 = candidate->source.getPosition();
	
// 	double nMedium = mediumDensity->getDensity(pos, z);
// 	double nBeam = flowProperties->getDensity(pos, z);
// 	double L = beamDensityToLuminosity(nBeam, E, z);
// 	double T = mediumTemperature->getTemperature(pos, z);
// 	double ve = mediumTemperature->getVelocity(id, pos - pos0, z);
// 	double vi = mediumTemperature->getVelocity(crpropa::nucleusId(1, 1), pos - pos0, z); 	// consider only protons
 
// 	double wp = plasmaFrequency(nMedium, id);
// 	double wnl = wp / (nMedium * T) * (ve * ve / vi / crpropa::c_light);
// 	double nuc = 8.21e-12 * (nMedium / 1e-1) * pow(T / 1e4, -1.5);
// 	double delta = nuc / 2. / wnl / (nBeam * E);
// 	double wmax = maximumLinearGrowthFrequency(nBeam, nMedium, lf, id);
// 	double tau_1 = 4 * delta * wmax;

// 	return E * tau_1 / crpropa::c_light;
// }

/***************************************************************************/
/**/
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