#include "grplinst/Simulation.h"


namespace grplinst {


/***************************************************************************/
/**/
Scenario::Scenario(double alpha, double Emax, double L, crpropa::ref_ptr<EmissionGeometry> geo, double z, double alpha0, std::string cutoff) {
	setSpectralIndex(alpha);
	setSpectralIndexSimulation(alpha0);
	setEnergyCutoff(Emax);
	setLuminosity(L);
	setSourceRedshift(z);
	setCutoff(cutoff);
	setGeometry(geo);
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

void Scenario::setGeometry(crpropa::ref_ptr<EmissionGeometry> geo) {	
	// std::cout << geometry.get()->getHeight() << std::endl;
	// std::cout << c->getHeight();
	// std::cout << typeid(c).name() << std::endl
	geometry = geo;
}

void Scenario::setCutoff(std::string c) {
	if (c != "exp" || c != "sharp" || c != "expBroken" || c != "exp2" || c != "exp1") {
		std::invalid_argument("Cutoff shape unknown. Use: exp (exp1), exp2, sharp, or expBroken.");
	}
	if (c == "exp")
		c = "exp1";
	cutoff = c;
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

std::string Scenario::getCutoff() const {
	return cutoff;
}

crpropa::ref_ptr<EmissionGeometry> Scenario::getGeometry() const {
	return geometry;
}

double Scenario::computeVolume() const {
	return geometry->computeVolume();
}

double Scenario::computeWeight(double energy0) const {
	double w = pow(energy0, spectralIndexSimulation - spectralIndex);
	if (cutoff == "exp") {
		w *= exp(- energy0 / energyCutoff);
	} else if (cutoff == "expBroken") {
		if (energy0 > energyCutoff)
			w *= exp(- energy0 / energyCutoff);
	} else if (cutoff == "exp2") {
		w *= exp(- crpropa::pow_integer<2>(energy0 / energyCutoff));
	}
	return w;
}

// double Scenario::computeEffectiveLuminosity() const {
// 	geo.
// }


/***************************************************************************/
/**/
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

/***************************************************************************/
/**/
Simulation::Simulation() {
	scenarios.clear();
	observables.clear();
	setFilename("");
	setDistance(0);
	setMagneticField(0);
	setCoherenceLength(0);
}

Simulation::Simulation(std::string fn, double dist, double B, double LB) {
	scenarios.clear();
	observables.clear();
	setFilename(fn);
	setDistance(dist);
	setMagneticField(B);
	setCoherenceLength(LB);
}

Simulation::~Simulation() {
}

void Simulation::setDistance(double d) {
	distance = d;
}

void Simulation::setMagneticField(double B) {
	magneticField = B;
}

void Simulation::setCoherenceLength(double L) {
	coherenceLength = L;
}

void Simulation::setFilename(std::string fn) {
	filename = fn;
}

void Simulation::addScenario(crpropa::ref_ptr<Scenario> scenario) {
	crpropa::ref_ptr<EmissionGeometry> geom = nullptr;
	
	if (scenario->getGeometry()->getShapeName() == "cone") {
		Cone *c = static_cast<Cone*>(crpropa::get_pointer(scenario->getGeometry()));
		double angle = c->getAngle();
		crpropa::Vector3d axis = c->getAxis();
		c->setHeight(distance);
		// geom = c;
		geom = new Cone(angle, distance, axis);
	} else {
		std::runtime_error("Distance correction for geometry only implemented for cone.");
	}

	double a = scenario->getSpectralIndex();
	double a0 = scenario->getSpectralIndexSimulation();
	double Emax = scenario->getEnergyCutoff();
	double L = scenario->getLuminosity();
	double z = scenario->getSourceRedshift();
	std::string cutoff = scenario->getCutoff();

	scenarios.push_back(new Scenario(a, Emax, L, geom, z, a0, cutoff));
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

	for (size_t k = 0; k < scenarios.size(); k++) {
		observables.push_back(new EmissionObservables());
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
		

		for (size_t k = 0; k < scenarios.size(); k++) {
			// reweight spectrum
			double w = w0 * scenarios[k]->computeWeight(e0);
			if (fabs(id) == 11) { // only electrons
				// jet geometry
				// w *=

				// compute gamma from angles?
				double gamma = e / (crpropa::mass_electron * crpropa::c_light * crpropa::c_light);
				// double theta = p.getAngleTo(p0);
				// double gamma = 

				observables[k]->incrementWeightSum(w);
				observables[k]->incrementEnergyTotal0(w * e0);
				observables[k]->incrementEnergyTotal(w * e);
				observables[k]->incrementMeanLorentzFactor(w * gamma);
				observables[k]->incrementMeanInverseLorentzFactor(w / gamma);
			}
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
		double Npair = norm * (Etot * wsum) / (lf * crpropa::mass_electron * crpropa::pow_integer<2>(crpropa::c_light) / wsum) / 2.;
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

double Simulation::getDistance() const {
	return distance;
}

double Simulation::getMagneticField() const {
	return magneticField;
}

double Simulation::getCoherenceLength() const {
	return coherenceLength;
}

/***************************************************************************/
/**/
void saveEmissionProfile(std::vector<crpropa::ref_ptr<Simulation>> simulations, std::string filename) {
	// check if simulations are of same size
	size_t Nscenarios = simulations[0]->getScenarios().size();
	for (size_t k = 1; k < simulations.size(); k++) {
		if (simulations[k]->getScenarios().size() != Nscenarios) {
			std::length_error("To create the emission profiles, the scenarios should be the same.");
		}
		if (simulations[k]->getMagneticField() != simulations[0]->getMagneticField() || simulations[k]->getCoherenceLength() != simulations[0]->getCoherenceLength()) {
			std::invalid_argument("The vector of simulations provided should be for the same magnetic field.");
		}
	}
	size_t Nsim = simulations.size();

	// work on filenames
	size_t idx = filename.find_last_of(".");
	std::string fn = filename.substr(0, idx);
	std::string extension = filename.substr(idx + 1);

	for (size_t j = 0; j < Nscenarios; j++) { // loop over scenarios
		crpropa::ref_ptr<Scenario> scenario = simulations[0]->getScenarios()[j];
		
		// output file
		std::string suffixStr = formatString("z_%4.3f-B_%1.0eG-LB_%1.0ekpc-L_%1.0eW-alpha_%2.1f-Emax_%2.1eeV", scenario->getSourceRedshift(), simulations[0]->getMagneticField() / crpropa::gauss, simulations[0]->getCoherenceLength() / crpropa::kpc, scenario->getLuminosity(), scenario->getSpectralIndex(), scenario->getEnergyCutoff() / crpropa::eV);
		std::string outputFile = fn + "-" + suffixStr + "." + extension;
		std::ofstream file;
		file.open(outputFile, std::ios::out);
		file << "# D [m]	n_e [1/m^3] 	<gamma> 	<1/gamma>\n";
		file.setf(std::ios::showpoint);
		file.setf(std::ios::scientific);
		file.width(8);

		for (size_t k = 0; k < simulations.size(); k++) { // loop over simulations
			crpropa::ref_ptr<Simulation> sim = simulations[k];
			crpropa::ref_ptr<EmissionObservables> observables = sim->getEmissionObservables()[j];
			double D = sim->getDistance();
			double n = observables->getDensity();
			double lf = observables->getMeanLorentzFactor();
			double ilf = observables->getMeanInverseLorentzFactor();
			file << D << "\t" << n << "\t" << lf << "\t" << ilf << '\n';
		}
		file.close();
	}
}




} // namespace grplinst