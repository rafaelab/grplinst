#include "grplinst/PlasmaInstability.h"


namespace grplinst {

// PlasmaInstability::PlasmaInstability(Flow *flow, MediumDensity *density, MediumTemperature *temperature, double limit) {
PlasmaInstability::PlasmaInstability(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) {
	setFlowProperties(flow);
	setMediumDensity(density);
	setMediumTemperature(temperature);
	setLimit(limit);
	setDescription("PlasmaInstability::PlasmaInstability");
}

PlasmaInstability::~PlasmaInstability() {
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

void PlasmaInstability::process(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double dx = candidate->getCurrentStep() / (1 + z);

	double dEdx = energyLoss(candidate);
	if (dEdx < 0) // prevent overshooting
		dEdx = 0;

	double Enew = E - dEdx * dx;

	candidate->current.setEnergy(Enew / (1 + z));
	candidate->limitNextStep(limit * E / dEdx);
}

/**************************************************/

PlasmaInstabilityBroderick2012B::PlasmaInstabilityBroderick2012B(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilityBroderick2012C");
}

double PlasmaInstabilityBroderick2012B::energyLoss(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return 0;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z); 
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition(); 
	double lf = candidate->current.getLorentzFactor() * (1 + z);

	double nMedium = mediumDensity->getValue(pos, z);
	double nBeam = flowProperties->getValue(pos, z);

	double wp = plasmaFrequency(nMedium, id);
	double tau_1 = 0.4 * lf * nBeam / nMedium * wp;

	return E * tau_1 / crpropa::c_light;
}

PlasmaInstabilityBroderick2012C::PlasmaInstabilityBroderick2012C(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilityBroderick2012C");
}

double PlasmaInstabilityBroderick2012C::energyLoss(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return 0;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy(); // multiply by (1 + z) for E(z)
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition(); 
	double lf = candidate->current.getLorentzFactor() * (1 + z);

	double nMedium = mediumDensity->getValue(pos, z);
	double nBeam = flowProperties->getValue(pos, z);

	double wp = plasmaFrequency(nMedium, id);
	double tau_1 = 0.4 * lf * nBeam / nMedium * wp;

	return E * tau_1 / crpropa::c_light;
}

PlasmaInstabilityVafin2018B::PlasmaInstabilityVafin2018B(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilityVafin2018B");
}

double PlasmaInstabilityVafin2018B::energyLoss(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return 0;
	
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z); // multiply by (1 + z) for E(z)
	double lf = candidate->current.getLorentzFactor() * (1 + z);
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();
	
	double nMedium = mediumDensity->getValue(pos, z);
	double nBeam = flowProperties->getValue(pos, z);
	double T = mediumTemperature->getValue(pos, z);
	double nRatio = (nBeam / 1e-20) / (nMedium / 1e-7);
	
	double wp = plasmaFrequency(nMedium, id);
	double wmax = maximumLinearGrowthFrequency(nBeam, nMedium, lf, id);
	double deltaMin = 6e-6 * (T / 1e4) * pow(lf / 1e6, -4. / 3.) * pow(nRatio, -2. / 3.);
	double delta = deltaMin * (1 + 4.7 * pow(lf / 1e6 * nRatio * nRatio, 2. / 3.));
	// double wmax = 1e-4 * wp * delta * nRatio * (lf / 1e6) / sqrt(T / 1e4);
	double tau_1 = 2 * wmax * delta;

	// tau_1 = 6e-11 * (T / 1e4) * pow(lf / 1e6, -4. / 3.) * pow(nRatio, 1. / 3.);

	return E * tau_1 / crpropa::c_light;
}

PlasmaInstabilityVafin2018C::PlasmaInstabilityVafin2018C(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilityVafin2018C");
}

double PlasmaInstabilityVafin2018C::energyLoss(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return 0;
	
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z); // multiply by (1 + z) for E(z)
	double lf = candidate->current.getLorentzFactor() * (1 + z);
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();
	
	double nMedium = mediumDensity->getValue(pos, z);
	double nBeam = flowProperties->getValue(pos, z);
	double T = mediumTemperature->getValue(pos, z);

	double wp = plasmaFrequency(nMedium, id);
	double wmax = maximumLinearGrowthFrequency(nBeam, nMedium, lf, id);
	// double wmax = wp * (nBeam / nMedium) * lf;

	// assumes <γ_g> in eq. 27 is γ_5
	double delta = 6e-6 * (T / 1e4) * pow(lf / 1e5, -4. / 3.) * pow((nMedium / 1e-7) / (nBeam / 1e-22), 2. / 3.) * (1. + 2.2e-3 * pow(lf / 1e5, 2. / 3.) * pow((nBeam / 1e-22) / (nMedium / 1e-7), 4. / 3.));

	// // uses γ_g as angular
	// double gamma_g = 1 + 
	// double delta = 6e-6 * (T / 1e4) * pow(lf / 1e5, -4. / 3.) * pow((nMedium / 1e-7) / (nBeam / 1e-22), 2. / 3.) * (1. + 2.2e-3 * pow(lf / 1e5, 2. / 3.) * pow((nBeam / 1e-22) / (nMedium / 1e-7), 4. / 3.));

	double tau_1 = 4 * wmax * delta;

	return E * tau_1 / crpropa::c_light;
}

PlasmaInstabilitySchlickeiser2012B::PlasmaInstabilitySchlickeiser2012B(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilitySchlickeiser2012B");
}

double PlasmaInstabilitySchlickeiser2012B::energyLoss(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return 0;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z); // multiply by (1 + z) for E(z)
	double lf = candidate->current.getLorentzFactor() * (1 + z);
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();

	double nMedium = mediumDensity->getValue(pos, z);
	double nBeam = flowProperties->getValue(pos, z);
	double T = mediumTemperature->getValue(pos, z);

	double wp = plasmaFrequency(nMedium, id);
	double wmax = maximumLinearGrowthFrequency(nBeam, nMedium, lf, id);

	double tau_e_1 = 1.1e-6 * pow(lf / 1e6, -1. / 3) * pow(nMedium / 0.1, 1. / 6.) * pow(nBeam / 1e-16, 1. / 3.);
	double tau = 7.4 / tau_e_1 * (1. + 1.25 * log(T / 1e4) - 0.25 * log(nMedium / 0.1));
	double tau_1 = 1 / tau;

	// double tau_1 = 1. / 7.4 * wmax * (1 + 1.25 * log(T / 1e4) - 0.25 * log(nMedium / 1e-7 / cm_3));

	// double taue_1 = 1.1e-6 * pow(lf / 1e6, -1. / 3.) * pow(nBeam / 1e-22 / cm_3, 1. / 3.) * pow(nMedium / 1e-7 / cm_3, 1. / 6.); 
	// double tau_1 = 1. / 7.4 * wmax * (1 + 1.25 * log(T / 1e4) - 0.25 * log(nMedium / 1e-7 / cm_3));

	return E * tau_1 / crpropa::c_light;
}

PlasmaInstabilitySchlickeiser2012C::PlasmaInstabilitySchlickeiser2012C(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilitySchlickeiser2012C");
}

double PlasmaInstabilitySchlickeiser2012C::energyLoss(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return 0;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z); // multiply by (1 + z) for E(z)
	double lf = candidate->current.getLorentzFactor() * (1 + z);
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();
	
	double nMedium = mediumDensity->getValue(pos, z);
	double nBeam = flowProperties->getValue(pos, z);
	double T = mediumTemperature->getValue(pos, z);

	double wp = plasmaFrequency(nMedium, id);
	double wmax = maximumLinearGrowthFrequency(nBeam, nMedium, lf, id);
	double tau_1 = 1. / 7.4 * wmax * (1 + 1.25 * log(T / 1e4) - 0.25 * log(nMedium / 1e-1));

	return E * tau_1 / crpropa::c_light;
}

PlasmaInstabilityMiniati2013B::PlasmaInstabilityMiniati2013B(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilityMiniati2013B");
}

double PlasmaInstabilityMiniati2013B::energyLoss(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return 0;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double lf = candidate->current.getLorentzFactor() * (1 + z);
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();
	crpropa::Vector3d dir = candidate->current.getDirection();
	crpropa::Vector3d dir0 = candidate->source.getDirection();
	
	double nMedium = mediumDensity->getValue(pos, z);
	double nBeam = flowProperties->getValue(pos, z);
	double L = beamDensityToLuminosity(nBeam, E, z);
	double T = mediumTemperature->getValue(pos, z);
	double ve = mediumTemperature->getVelocity(id, pos, z);
	double vi = mediumTemperature->getVelocity(crpropa::nucleusId(1, 1), pos, z); 	// consider only protons

	double dTheta = (dir - dir0).getR();
	if (dTheta == 0)
		dTheta = 1e-3;

	double wp = plasmaFrequency(nMedium, id);
	double gmax = wp * (nBeam / nMedium) * 1 / (dTheta * dTheta) * lf; // Lorentz factor should be in denominator; this results in the right slope, though
	// double gmax = wp * (nBeam / nMedium) * 1 / (dTheta * dTheta * lf);
	double gnl = wp / (nMedium * T * crpropa::k_boltzmann) * (ve * ve / vi / crpropa::c_light);
	double nuc = 1e-11 * (nMedium / 0.02) * pow(T / 3e3, -1.5);
	double qr = nuc / 2. / gnl;
	double qnr = gmax / gnl;
	double tau = nBeam / 2 / gmax / qr;

	return E / tau / crpropa::c_light;
}

PlasmaInstabilityMiniati2013C::PlasmaInstabilityMiniati2013C(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilityMiniati2013C");
}

double PlasmaInstabilityMiniati2013C::energyLoss(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return 0;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z); // multiply by (1 + z) for E(z)
	double lf = candidate->current.getLorentzFactor() * (1 + z);
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();
	
	double nMedium = mediumDensity->getValue(pos, z);
	double nBeam = flowProperties->getValue(pos, z);
	double L = beamDensityToLuminosity(nBeam, E, z);
	double T = mediumTemperature->getValue(pos, z);
	double ve = mediumTemperature->getVelocity(id, pos - pos0, z);
	double vi = mediumTemperature->getVelocity(crpropa::nucleusId(1, 1), pos - pos0, z); 	// consider only protons
 
	double wp = plasmaFrequency(nMedium, id);
	double wnl = wp / (nMedium * T) * (ve * ve / vi / crpropa::c_light);
	double nuc = 8.21e-12 * (nMedium / 1e-1) * pow(T / 1e4, -1.5);
	double delta = nuc / 2. / wnl / (nBeam * E);
	double wmax = maximumLinearGrowthFrequency(nBeam, nMedium, lf, id);
	double tau_1 = 4 * delta * wmax;

	return E * tau_1 / crpropa::c_light;
}

/**************************************************/


PlasmaInstabilityBroderick2012::PlasmaInstabilityBroderick2012(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilityBroderick2012");
}

double PlasmaInstabilityBroderick2012::energyLoss(crpropa::Candidate *candidate) const {
	
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy(); // multiply by (1 + z) for E(z)
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();
	
	// Redshift not included here: it is accounted for in the equation
	double nMedium = mediumDensity->getValue(pos - pos0);
	double nBeam = flowProperties->getValue(pos - pos0);
	double L = beamDensityToLuminosity(nBeam, E, z);

	double eta = 1.;
	double Ethr = 8.7e-6 * pow(1 + z, -13. / 6.) * pow(L / 1e38, -1. / 3.) * pow(nMedium / 0.1, 1. / 3.);
	double a0, a1, a2, a3, a4;
	if (E < Ethr) {
		a0 = 7.7e-26;
		a1 = 8;
		a2 = 3;
		a3 = 1;
		a4 = -1. / 2;
	} else {
		a0 = 2.3e-22;
		a1 = 11. / 3;
		a2 = 1;
		a3 = 1. / 3;
		a4 = 1. / 6;
	}

	return a0 * eta * pow(E / crpropa::TeV, a2) * pow(L / 1e38, a3) * pow(nMedium / 0.1, a4);
}

PlasmaInstabilityMiniati2013::PlasmaInstabilityMiniati2013(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	initTable();
	setDescription("PlasmaInstability::PlasmaInstabilityMiniati2013");
}

void PlasmaInstabilityMiniati2013::initTable() {
	double d_[] = {-0.05, 0.14, 0.35, 0.59, 0.79, 0.96, 1.17, 1.40, 1.57, 1.77, 1.99, 2.20, 2.41, 2.60, 2.80, 3.00};
	double w_[] = { 3.28, 3.46, 3.14, 3.08, 3.05, 3.00, 2.84, 2.56, 2.40, 2.55, 2.36, 2.08, 1.73, 1.55, 1.05, 0.75};
	std::vector<double> D (d_, d_ + sizeof(d_) / sizeof(double));
	std::vector<double> W (w_, w_ + sizeof(w_) / sizeof(double));
	_d = D;
	_w = W;
}

double PlasmaInstabilityMiniati2013::energyLoss(crpropa::Candidate *candidate) const {
	
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy(); // multiply by (1 + z) for E(z)
	
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();

	// double nMedium = mediumDensity->getValue(pos);
	// double nBeam = flowProperties->getValue(pos);
	// double L = beamDensityToLuminosity(nBeam, E, z);

	// medium and flow properties useless? 
	// recheck paper

	E /= (1 + z);
	double d0 = log10(crpropa::redshift2LightTravelDistance(z) / crpropa::Mpc); // co-moving?
	double w = pow(10, crpropa::interpolate(d0, _d, _w));

	return 1.4e-29 * crpropa::pow_integer<2>(1 + z) * crpropa::pow_integer<2>(E / crpropa::TeV) / w;
}

PlasmaInstabilitySchlickeiser2012::PlasmaInstabilitySchlickeiser2012(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilitySchlickeiser2012");
}

double PlasmaInstabilitySchlickeiser2012::energyLoss(crpropa::Candidate *candidate) const {
	
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy(); // multiply by (1 + z) for E(z)
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();
	
	double nMedium = mediumDensity->getValue(pos - pos0);
	double nBeam = flowProperties->getValue(pos - pos0);
	double L = beamDensityToLuminosity(nBeam, E, z);
	double T = mediumTemperature->getValue(pos);

	double eta = 1.;
	double Ethr = 7.9e-8 * pow(1 + z, -9. / 4) / sqrt(L / 1e38) * sqrt(nMedium / 0.1) * (T / 1e4);
	double F = 1. + 1.25 * log(T / 1e4) - 0.25 * log(nMedium / 0.1) + 0.5 * log(1 + z);

	double a0, a1, a2, a3, a4, b;
	if (E < Ethr) {
		a0 = 4.7e-30;
		a1 = -5;
		a2 = -1;
		a3 = 1. / 3.;
		a4 = 5. / 6.;
		b = crpropa::pow_integer<2>(T / 1e4);
	} else {
		a0 = 1.4e-23;
		a1 = 11. / 3;
		a2 = 1;
		a3 = 1. / 3.;
		a4 = 1. / 6.;
		b = 1. / F;
	}

	return a0 * eta * pow(1 + z, a1) * pow(E / crpropa::TeV, a2) * pow(L / 1e38, a3) * pow(nMedium / 0.1, a4) * b;
}

PlasmaInstabilitySironi2014::PlasmaInstabilitySironi2014(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilitySironi2014");
}

double PlasmaInstabilitySironi2014::energyLoss(crpropa::Candidate *candidate) const {
	
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy(); // multiply by (1 + z) for E(z)
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();
	
	double nMedium = mediumDensity->getValue(pos - pos0);
	double nBeam = flowProperties->getValue(pos - pos0);
	double L = beamDensityToLuminosity(nBeam, E, z);
	double T = mediumTemperature->getValue(pos);

	double eta = 1.; 
	double Ethr = 6.9e-6 * pow(1 + z, -13. / 16) * pow(L / 1e38, -1. / 3) / sqrt(nMedium / 0.1);

	double a0, a1, a2, a3, a4, b;
	if (E < Ethr) {
		a0 = 3.9e-25;
		a1 = 8.;
		a2 = 3.;
		a3 = 1.;
		a4 = -0.5;
	} else {
		a0 = 1.2e-22;
		a1 = 11. / 3;
		a2 = 1.;
		a3 = 1. / 3;
		a4 = 1. / 6;
	}

	return a0 * eta * pow(1 + z, a1) * pow(E / crpropa::TeV, a2) * pow(L / 1e38, a3) * pow(nMedium / 0.1, a4);
}

PlasmaInstabilityVafin2018::PlasmaInstabilityVafin2018(crpropa::ref_ptr<Flow> flow, crpropa::ref_ptr<MediumDensity> density, crpropa::ref_ptr<MediumTemperature> temperature, double limit) : PlasmaInstability(flow, density, temperature, limit) {
	setDescription("PlasmaInstability::PlasmaInstabilityVafin2018");
}

double PlasmaInstabilityVafin2018::energyLoss(crpropa::Candidate *candidate) const {
	
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy(); // multiply by (1 + z) for E(z)
	crpropa::Vector3d pos = candidate->current.getPosition();
	crpropa::Vector3d pos0 = candidate->source.getPosition();
	
	double nMedium = mediumDensity->getValue(pos - pos0);
	double nBeam = flowProperties->getValue(pos - pos0);
	double L = beamDensityToLuminosity(nBeam, E, z);
	double T = mediumTemperature->getValue(pos);

	return 2.7e-20 * pow(0.5 + z / 2., 19. / 6) * E * pow(E / crpropa::TeV, -1.) * pow(L / 1e38, 1. / 3) * pow(nMedium / 0.1,  -1. / 3) * (T / 1e4);
}





double maximumLinearGrowthFrequency(double beamDensity, double mediumDensity, double lorentzFactor, int id) {
	double wp =  plasmaFrequency(mediumDensity, id);
	return wp * beamDensity / mediumDensity * lorentzFactor;
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