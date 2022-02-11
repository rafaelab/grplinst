#include "grplinst/PlasmaInstability.h"


namespace grplinst {

PlasmaInstability::PlasmaInstability(MediumDensity *density, MediumTemperature *temperature, Flow *flow, const std::string model, double limit) {
	setMediumDensity(density);
	setMediumTemperature(temperature);
	setFlowProperties(flow);
	setInstabilityModel(model);
	setDescription("PlasmaInstability::PlasmaInstability");
}

PlasmaInstability::~PlasmaInstability()  {
	delete mediumDensity;
	delete mediumTemperature;
	delete flowProperties;
}

void PlasmaInstability::setMediumDensity(MediumDensity *density) {
	mediumDensity = density;
}

void PlasmaInstability::setMediumTemperature(MediumTemperature *temperature) {
	mediumTemperature = temperature;
}

void PlasmaInstability::setFlowProperties(Flow *flow) {
	flowProperties = flow;
}

void PlasmaInstability::setInstabilityModel(const std::string model) {
	plasmaInstabilityModel = model;
}

void PlasmaInstability::setLimit(double l) {
	limit = l;
}

crpropa::ref_ptr<MediumDensity> PlasmaInstability::getMediumDensity() const {
	return mediumDensity;
}

crpropa::ref_ptr<MediumTemperature> PlasmaInstability::getMediumTemperature() const {
	return mediumTemperature;
}

crpropa::ref_ptr<Flow> PlasmaInstability::getFlowProperties() const {
	return flowProperties;
}

std::string PlasmaInstability::getInstabilityModel() const {
	return plasmaInstabilityModel;
}

void PlasmaInstability::process(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	if (fabs(id) != 11)
		return;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double dx = candidate->getCurrentStep() / (1 + z);

	double dEdx = energyLoss(candidate);
	if (dEdx < 0)
		dEdx = 0;

	candidate->current.setEnergy(E - dEdx * dx);
	candidate->limitNextStep(limit * E / dEdx);
}

// void PlasmaInstability::process(Candidate *candidate) const {

// 	int id = candidate->current.getId();

// 	// only works for electrons and positrons
// 	if (fabs(id) != 11)
// 		return;

// 	double dx = candidate->getCurrentStep();
// 	double z = candidate->getRedshift();
// 	double E = candidate->current.getEnergy();
// 	double dEdx = coolingPower(E, z);

// 	if (dEdx < 0)
// 		dEdx = 0;
// 	double Enew = E - dEdx * dx;
// 	candidate->current.setEnergy(Enew);
// 	candidate->limitNextStep(0.1 * E / dEdx);
// }

// double PlasmaInstability::coolingPower(double E, double z) const {

// 	HomogeneousBeam *beamType = dynamic_cast<HomogeneousBeam*>(&beam);
// 	if!(beamType) {
// 		std::cout << "Only homogeneous beams are supported for now" << std::endl;
// 		return;
// 	}

// 	double dEdx = 0;
// 	if(plasmaInstabilityModel == "A" or plasmaInstabilityModel == "a") {
// 		dEdx = coolingPowerA(E, z);
// 	} else if(plasmaInstabilityModel == "B" or plasmaInstabilityModel == "b") {
// 		dEdx = coolingPowerB(E, z);
// 	} else if(plasmaInstabilityModel == "C" or plasmaInstabilityModel == "c") {
// 		dEdx = coolingPowerC(E, z);
// 	} else if(plasmaInstabilityModel == "D" or plasmaInstabilityModel == "d") {
// 		dEdx = coolingPowerD(E, z);
// 	} else if(plasmaInstabilityModel == "E" or plasmaInstabilityModel == "e") {
// 		dEdx = coolingPowerE(E, z);
// 	}
// 	return dEdx;
// }

// double PlasmaInstability::coolingPowerA(double E, double z) const {

// 	E /= (1 + z);
// 	double eta = 1.;
// 	double Ethr = 8.7e-6 * pow(1 + z, -13. / 6) * pow(luminosityBeam / 1e38, -1. / 3) * pow(densityIGM / 0.1, 1. / 3);

// 	double a0, a1, a2, a3, a4;
// 	if (E < Ethr) {
// 		a0 = 7.7e-26;
// 		a1 = 8;
// 		a2 = 3;
// 		a3 = 1;
// 		a4 = -1. / 2;
// 	} else {
// 		a0 = 2.3e-22;
// 		a1 = 11. / 3;
// 		a2 = 1;
// 		a3 = 1. / 3;
// 		a4 = 1. / 6;
// 	}
// 	return a0 * eta * pow(E / TeV, a2) * pow(luminosityBeam / 1e38, a3) * pow(densityIGM / 0.1, a4);
// }

// double PlasmaInstability::coolingPowerB(double E, double z) const {
	
// 	E /= (1 + z);
// 	double d0 = log10(redshift2LightTravelDistance(z) / Mpc); // co-moving?
// 	double w = pow(10, interpolate(d0, _d, _w));
// 	return 1.4e-29 * pow(1 + z, 2) * pow(E / TeV, 2) / w;
// }

// double PlasmaInstability::coolingPowerC(double E, double z) const {

// 	E /= (1 + z);
// 	double eta = 1.;
// 	double Ethr = 7.9e-8 * pow(1 + z, -9. / 4) / sqrt(luminosityBeam / 1e38) * pow(densityIGM / 0.1, 0.5) * pow(temperatureIGM / 1e4, 1.);
// 	double F = 1. + 1.25 * log(temperatureIGM / 1e4) - 0.25 * log(densityIGM / 0.1) + 0.5 * log(1 + z);

// 	double a0, a1, a2, a3, a4, b;
// 	if (E < Ethr) {
// 		a0 = 4.7e-30;
// 		a1 = -5;
// 		a2 = -1;
// 		a3 = 1. / 3.;
// 		a4 = 5. / 6.;
// 		b = pow(temperatureIGM / 1e4, 2);
// 	} else {
// 		a0 = 1.4e-23;
// 		a1 = 11. / 3;
// 		a2 = 1;
// 		a3 = 1. / 3.;
// 		a4 = 1. / 6.;
// 		b = 1. / F;
// 	}
// 	return a0 * eta * pow(1 + z, a1) * pow(E / TeV, a2) * pow(luminosityBeam / 1e38, a3) * pow(densityIGM / 0.1, a4) * b;
// }

// double PlasmaInstability::coolingPowerD(double E, double z) const {

// 	E /= (1 + z);
// 	double eta = 1.; // for now fixed; should add a function to play with it.
// 	double Ethr = 6.9e-6 * pow(1 + z, -13. / 16) * pow(luminosityBeam / 1e38, -1. / 3) / sqrt(densityIGM / 0.1);

// 	double a0, a1, a2, a3, a4, b;
// 	if (E < Ethr) {
// 		a0 = 3.9e-25;
// 		a1 = 8.;
// 		a2 = 3.;
// 		a3 = 1.;
// 		a4 = -0.5;
// 	} else {
// 		a0 = 1.2e-22;
// 		a1 = 11. / 3;
// 		a2 = 1.;
// 		a3 = 1. / 3;
// 		a4 = 1. / 6;
// 	}
// 	return a0 * eta * pow(1 + z, a1) * pow(E / TeV, a2) * pow(luminosityBeam / 1e38, a3) * pow(densityIGM / 0.1, a4);
// }

// double PlasmaInstability::coolingPowerE(double E, double z) const {

// 	E /= (1 + z);
// 	return 2.7e-20 * pow(0.5 + z / 2., 19. / 6) * E * pow(E / TeV, -1.) * pow(luminosityBeam / 1e38, 1. / 3) * pow(densityIGM / 0.1,  -1. / 3) * (temperatureIGM / 1e4);

// }

} // namespace grplinst