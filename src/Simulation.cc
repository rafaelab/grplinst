#include "grplinst/Simulation.h"


namespace grplinst {


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
EmissionGeometry::EmissionGeometry() {
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
Cone::Cone(double angle, double height, crpropa::Vector3d axis) {
	setAngle(angle);
	setAxis(axis);
	setHeight(height);
	setRadius(computeRadius());
}

void Cone::setAngle(double a) {
	angle = a;
}

void Cone::setHeight(double h) {
	height = h;
}

void Cone::setRadius(double r) {
	radius = r;
}

void Cone::setAxis(crpropa::Vector3d a) {
	axis = a;
}

double Cone::computeRadius() const {
	return tan(angle) * height;
}

double Cone::computeArea() const {
	return M_PI * radius * (radius + sqrt(radius * radius + height * height));
}

double Cone::computeVolume() const {
	return M_PI * radius * radius * height / 3.;
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
Scenario::Scenario(double alpha, double Emax, double L, double dist, crpropa::ref_ptr<EmissionGeometry> geo, double alpha0, double z) {
	setSpectralIndex(alpha);
	setSpectralIndexSimulation(alpha0);
	setEnergyCutoff(Emax);
	setLuminosity(L);
	setSourceRedshift(z);
	setGeometry(geo);
	setDistance(dist);
}

void Scenario::setSpectralIndex(double s) {
	spectralIndex = s;
}

void Scenario::setSpectralIndexSimulation(double s) {
	spectralIndexSimulation = s;
}

void Scenario::setEnergyCutoff(double Emax) {
	energyCutoff = Emax;
}

void Scenario::setLuminosity(double L) {
	luminosity = L;
}

void Scenario::setSourceRedshift(double z) {
	redshift = z;
}

void Scenario::setDistance(double d) {
	distance = d;
}

void Scenario::setGeometry(crpropa::ref_ptr<EmissionGeometry> geo) {
	geometry = geo;
}

double Scenario::getSpectralIndex() const {
	return spectralIndex;
}

double Scenario::getSpectralIndexSimulation() const {
	return	spectralIndexSimulation;
}

double Scenario::getEnergyCutoff() const {
	return energyCutoff;
}

double Scenario::getLuminosity() const {
	return luminosity;
}

double Scenario::getDistance() const {
	return distance;
}

double Scenario::getSourceRedshift() const {
	return redshift;
}

double Scenario::getSourceLightTravelDistance() const {
	if (redshift == 0) {
		std::cout << "You are trying to compute a distance for redshift 0." << std::endl;
		return 0.;
	}
	return crpropa::redshift2LightTravelDistance(redshift);
}

double Scenario::getSourceComovingDistance() const {
	if (redshift == 0) {
		std::cout << "You are trying to compute a distance for redshift 0." << std::endl;
		return 0.;
	}
	return crpropa::redshift2ComovingDistance(redshift);
}

double Scenario::getSourceLuminosityDistance() const {
	if (redshift == 0) {
		std::cout << "You are trying to compute a distance for redshift 0." << std::endl;
		return 0.;
	}
	return crpropa::redshift2LuminosityDistance(redshift);
}

crpropa::ref_ptr<EmissionGeometry> Scenario::getGeometry() const {
	return geometry;
}

double Scenario::computeVolume() const {
	return geometry->computeVolume();
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
EmissionObservables::EmissionObservables() {
	setMeanLorentzFactor(0);
	setMeanInverseLorentzFactor(0);
	setWeightSum(0);
	setEnergyTotal(0);
	setEnergyTotal0(0);
	setDensity(0);
}

EmissionObservables::EmissionObservables(double lf, double ilf, double Etot, double Etot0, double wtot, double n) {
	setMeanLorentzFactor(lf);
	setMeanInverseLorentzFactor(ilf);
	setWeightSum(wtot);
	setEnergyTotal(Etot);
	setEnergyTotal0(Etot0);
	setDensity(n);
}

void EmissionObservables::setMeanLorentzFactor(double lf) {
	meanLorentzFactor = lf;
}

void EmissionObservables::setMeanInverseLorentzFactor(double ilf) {
	meanInverseLorentzFactor = ilf;
}

void EmissionObservables::setEnergyTotal(double E) {
	energyTotal = E;
}

void EmissionObservables::setEnergyTotal0(double E0) {
	energyTotal0 = E0;
}

void EmissionObservables::setWeightSum(double w) {
	weightSum = w;
}

void EmissionObservables::setDensity(double n) {
	density = n;
}

void EmissionObservables::incrementMeanLorentzFactor(double lf) {
	meanLorentzFactor += lf;
}

void EmissionObservables::incrementMeanInverseLorentzFactor(double ilf) {
	meanInverseLorentzFactor += ilf;
}

void EmissionObservables::incrementEnergyTotal(double E) {
	energyTotal += E;
}

void EmissionObservables::incrementEnergyTotal0(double E0) {
	energyTotal0 += E0;
}

void EmissionObservables::incrementWeightSum(double w) {
	weightSum += w;
}

double EmissionObservables::getWeightSum() const {
	return weightSum;
}

double EmissionObservables::getMeanLorentzFactor() const {
	return meanLorentzFactor;
}

double EmissionObservables::getMeanInverseLorentzFactor() const {
	return meanInverseLorentzFactor;
}

double EmissionObservables::getEnergyTotal() const {
	return energyTotal;
}

double EmissionObservables::getEnergyTotal0() const {
	return energyTotal0;
}

double EmissionObservables::getDensity() const {
	return density;
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

Simulation::Simulation() {
	setFilename("");
	scenarios.clear();
	observables.clear();
}

Simulation::Simulation(std::string fn) {
	setFilename(fn);
}

Simulation::~Simulation() {
}

void Simulation::setFilename(std::string fn) {
	filename = fn;
}

void Simulation::addScenario(crpropa::ref_ptr<Scenario> scenario) {
	scenarios.push_back(scenario);
}

void Simulation::process() {
	if (scenarios.size() == 0) {
		std::runtime_error("No scenarios available. Use `addScenario(scenario)` first.");
	}
	if (filename == "") {
		std::invalid_argument("File not provided.");
	} else {
		std::ifstream infile(filename);
		if (!infile.good())
			std::runtime_error("Simulation file not good: " + filename + ".");
	}

	std::vector<double> weight;
	std::vector<double> a;
	std::vector<double> a0;
	std::vector<double> Emax;
	std::vector<double> L;
	std::vector<crpropa::ref_ptr<EmissionGeometry>> geo;
	for (size_t k = 0; k < scenarios.size(); k++) {
		observables.push_back(new EmissionObservables());
		weight.push_back(0);
		a.push_back(scenarios[k]->getSpectralIndex());
		a0.push_back(scenarios[k]->getSpectralIndexSimulation());
		Emax.push_back(scenarios[k]->getEnergyCutoff());
		L.push_back(scenarios[k]->getLuminosity());
		geo.push_back(scenarios[k]->getGeometry());
	}

	// read file only once
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("GRPlInst could not open file " + filename);
	std::string line;
	while (std::getline(infile, line)) {
		std::istringstream lineStream(line);
		if (lineStream.peek() == '#')
			continue;
		double _;
		double e, e0;
		double px, py, pz;
		double p0x, p0y, p0z;
		double w0;
		int id;
		lineStream >> _ >> _ >> _ >> id >> e >> _ >> _ >> _ >> px >> py >> pz >> _ >> _ >> e0 >> _ >> _ >> _ >> p0x >> p0y >> p0z >> _ >> _ >> _ >> _ >> _ >> _ >> _ >> _ >> _ >> w0;
		e *= crpropa::eV;
		e0 *= crpropa::eV;
		
		// consider only electrons
		if (fabs(id) != 11)
			continue;

		for (size_t k = 0; k < scenarios.size(); k++) {
			// reweight spectrum
			weight[k] = w0 * pow(e0, a0[k] - a[k]);
			if (e0 > Emax[k])
				weight[k] *= exp(- e0 / Emax[k]);

			// jet geometry
			// w *=

			// compute gamma from angles?
			double gamma = e / (crpropa::mass_electron * crpropa::c_light * crpropa::c_light);
			// double theta = p.getAngleTo(p0);
			// double gamma = 

			observables[k]->incrementWeightSum(weight[k]);
			observables[k]->incrementEnergyTotal(weight[k] * e);
			observables[k]->incrementEnergyTotal0(weight[k] * e0);
			observables[k]->incrementMeanLorentzFactor(weight[k] * gamma);
			observables[k]->incrementMeanInverseLorentzFactor(weight[k] / gamma);
		}
	}
	infile.close();
	
	for (size_t k = 0; k < observables.size(); k++) {
		double lf = observables[k]->getMeanLorentzFactor();
		double ilf = observables[k]->getMeanInverseLorentzFactor();
		double wsum = observables[k]->getWeightSum();
		double Etot0 = observables[k]->getEnergyTotal0();
		double Etot = observables[k]->getEnergyTotal();
		double L = scenarios[k]->getLuminosity();
		double V = scenarios[k]->computeVolume();

		double norm = L / Etot0;
		double Npair = wsum * norm * Etot / (lf * crpropa::mass_electron * crpropa::pow_integer<2>(crpropa::c_light)) / 2.;
		double npair = Npair / V;
		observables[k]->setMeanLorentzFactor(lf / wsum);
		observables[k]->setMeanInverseLorentzFactor(ilf / wsum);
		observables[k]->setDensity(npair);
	}
}

std::vector<crpropa::ref_ptr<Scenario>> Simulation::getScenarios() const {
	return scenarios;
}

std::vector<crpropa::ref_ptr<EmissionObservables>> Simulation::getEmissionObservables() const {
	return observables;
}

std::string Simulation::getFilename() const {
	return filename;
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void saveEmissionProfile(std::vector<crpropa::ref_ptr<Simulation>> simulations, std::string filename) {
	// check if simulations are of same size
	size_t Nscenarios = simulations[0]->getScenarios().size();
	for (size_t k = 1; k < simulations.size(); k++) {
		if (simulations[k]->getScenarios().size() != Nscenarios)
			std::length_error("To create the emission profiles, the scenarios should be  the same.");
	}
	size_t Nsim = simulations.size();

	// work on filenames
	size_t idx = filename.find_last_of(".");
	std::string fn = filename.substr(0, idx);
	std::string extension = filename.substr(idx + 1);

	for (size_t j = 0; j < Nscenarios; j++) { // loop over scenarios
		crpropa::ref_ptr<Scenario> scenario = simulations[0]->getScenarios()[j];
		// output file
		std::stringstream suffixStr;
		suffixStr << std::fixed << std::setprecision(simulations.size() % 10 + 1) << j;
		std::string outputFile = fn + "-" + suffixStr.str() + extension;
		std::ofstream file;
		file.open(outputFile, std::ios::out);
		file << "# D [m]	n_e [1/m^3] 	<gamma> 	<1/gamma>\n";
		file.setf(std::ios::showpoint);
		file.setf(std::ios::scientific);
		file.width(8);
		for (size_t k = 0; k < simulations.size(); k++) { // loop over simulations
			crpropa::ref_ptr<Simulation> sim = simulations[k];
			crpropa::ref_ptr<EmissionObservables> observables = sim->getEmissionObservables()[j];
			double D = scenario->getDistance();
			double n = observables->getDensity();
			double lf = observables->getMeanLorentzFactor();
			double ilf = observables->getMeanInverseLorentzFactor();
			file << D << "\t" << n << "\t" << lf << "\t" << ilf << "\n";
		}
		file.close();
	}
}


} // namespace grplinst