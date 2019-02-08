#include "grplinst.h"


using namespace crpropa;

PlasmaInstability::PlasmaInstability(double nIGM, double temperature, double lumBeam, char model) : Module() {

	setIGMDensity(nIGM);
	setTemperature(temperature);
	setBeamLuminosity(lumBeam);
	setModel(model);
	initTables();
	setDescription("PlasmaInstability::PlasmaInstability");
}

void PlasmaInstability::initTables() {
	
	// Initiate tables for Miniati and Elyv model
	double d_[] = {-0.05, 0.14, 0.35, 0.59, 0.79, 0.96, 1.17, 1.40, 1.57, 1.77, 1.99, 2.20, 2.41, 2.60, 2.80, 3.00};
	double w_[] = { 3.28, 3.46, 3.14, 3.08, 3.05, 3.00, 2.84, 2.56, 2.40, 2.55, 2.36, 2.08, 1.73, 1.55, 1.05, 0.75};
	std::vector<double> D (d_, d_ + sizeof(d_) / sizeof(double));
	std::vector<double> W (w_, w_ + sizeof(w_) / sizeof(double));
	_d = D;
	_w = W;
}

void PlasmaInstability::setIGMDensity(double rho) {
	// Defines the density of the IGM (at z=z?, at z=0?)
	// in units of m^-3	
	densityIGM = rho;
}

void PlasmaInstability::setBeamLuminosity(double rho) {
	// Defines the density of the beam (at z=z?, at z=0?)
	// in units of m^-3
	luminosityBeam = rho;
}

void PlasmaInstability::setTemperature(double t) {
	// Defines the temperature of the IGM (at z=z?, at z=0?)
	temperatureIGM = t;
}

void PlasmaInstability::setModel(char model) {
	// Select the model of plasma instability.
	// Models available:
	//  A: Broderick, Chang, Pfrommer. Astrophys. J. 752 (2012) 22. arXiv:1106.5494
	//  B: Miniati & Elyiv. Astrophys. J. 770 (2013) 54. arXiv:1208.1761
	//  C: Schlickeiser, Ibscher, Supsar. Astrophys. J. 758 (2012) 102.
	//  D: Sironi, Giannios. Astrophys. J. 787 (2014) 49. arXiv:1312.4538
	//  E: Vafin, Rafighi, Pohl, Niemiec. Astrophys. J. 857 (2018) 43. arXiv:1803.02990
	plasmaInstabilityModel = model;
}

void PlasmaInstability::process(Candidate *candidate) const {

	int id = candidate->current.getId();
	// only works for electrons and positrons
	if (fabs(id) != 11)
		return;

	double dx = candidate->getCurrentStep();
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy();
	double dEdx = coolingPower(E, z);

	if (dEdx < 0)
		dEdx = 0;
	double Enew = E - dEdx * dx;
	candidate->current.setEnergy(Enew);
	candidate->limitNextStep(0.1 * E / dEdx);
}

double PlasmaInstability::coolingPower(double E, double z) const {

	double dEdx = 0;
	switch(plasmaInstabilityModel) {
		case 'a':
		case 'A':
			dEdx = coolingPowerA(E, z);
			break;
		case 'b':
		case 'B':
			dEdx = coolingPowerB(E, z);
			break;
		case 'c':
		case 'C':
			dEdx = coolingPowerC(E, z);
			break;
		case 'd':
		case 'D':
			dEdx = coolingPowerD(E, z);
			break;
		case 'e':
		case 'E':
			dEdx = coolingPowerE(E, z);
			break;
		default:
			dEdx = 0;
	}
	return dEdx;
}

double PlasmaInstability::coolingPowerA(double E, double z) const {

	double T = temperatureIGM / (1 + z);
	double nIGM = densityIGM / pow(1 + z, 3);
	double L = luminosityBeam;
	E /= (1 + z);
	double eta = 1.;
	double Ethr = 8.7e-6 * pow(1 + z, -13 / 6) * pow(L / 1e38, -1 / 3) * pow(nIGM / 0.1, 1 / 3);
	if (E < Ethr)
		return 7.7e-26 * eta * pow(1 + z, 8) * pow(E / TeV, 3) * (L / 1e38) * pow(nIGM / 0.1, -0.5);
	else
		return 2.3e-22 * eta * pow(1 + z, 11 / 3) * (E / TeV) * pow(L / 1e38, 1 / 3) * pow(nIGM / 0.1, 1 / 6);
}

double PlasmaInstability::coolingPowerB(double E, double z) const {

	double d0 = log10(redshift2LightTravelDistance(z) / Mpc); // co-moving?
	double w = pow(10, interpolate(d0, _d, _w));
	return 1.4e-29 * pow(1 + z, 2) * pow(E / TeV, 2) / w;
}

double PlasmaInstability::coolingPowerC(double E, double z) const {

	double T = temperatureIGM / (1 + z);
	double nIGM = densityIGM / pow(1 + z, 3);
	double L = luminosityBeam;
	E /= (1 + z);

	double eta = 1.;
	double Ethr = 7.9e-8 * pow(1 + z, -9 / 4) * pow(L / 1e38, -0.5) * pow(nIGM / 0.1, 0.5) * (T / 1e4);
	double F = 1 + 1.25 * log(T / 1e4) - 0.25 * log(nIGM / 0.1) + 0.5 * log(1 + z);
	if (E < Ethr)
		return 4.7e-30 * eta * pow(1 + z, -5) * pow(E / TeV, -1) * pow(L / 1e38, 1 / 3) * pow(nIGM / 0.1, 5 / 6) * pow(T / 1e4, 2);
	else
		return 1.4e-23 * eta * pow(1 + z, 11 / 3) * (E / TeV) * pow(L / 1e38, 1 / 3) * pow(nIGM / 0.1, 1 / 6) / F;
}

double PlasmaInstability::coolingPowerD(double E, double z) const {

	double T = temperatureIGM / (1 + z);
	double nIGM = densityIGM / pow(1 + z, 3);
	double L = luminosityBeam;

	E /= (1 + z);
	double eta = 1.; // for now fixed; should add a function to play with it.

	double Ethr = 6.9e-6 * pow(1 + z, -13 / 16) * pow(L / 1e38, -1/3) * pow(nIGM / 0.1, -0.5);
	if (E < Ethr)
		return 3.9e-25 * eta * pow(1 + z, 8) * pow(E / TeV, 3) * (L / 1e38) * pow(nIGM / 0.1, -1 / 2);
	else
		return 1.2e-22 * eta * pow(1 + z, 11 / 3) * (E / TeV) * pow(L / 1e38, 1 / 3) * pow(nIGM / 0.1, 1 / 6);
}

double PlasmaInstability::coolingPowerE(double E, double z) const {

	double T = temperatureIGM / (1 + z);
	double nIGM = densityIGM / pow(1 + z, 3);
	double L = luminosityBeam;
	E /= (1 + z);
	return 1.3e-26 * pow(1 + z, -13 / 6) * pow(E / TeV, -2 / 3) * pow(L / 1e38, -1 / 3) * pow(nIGM / 0.1,  1 / 3) * (T / 1e4);
}
