#ifndef GRPLINST_SIMULATION_H
#define GRPLINST_SIMULATION_H


#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <crpropa/Common.h>
#include <crpropa/Cosmology.h>
#include <crpropa/Grid.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>


namespace grplinst {

/**
Given a CRPropa output file of the beam, this function pre-processes it to obtain flow profiles.
Note that this pre-simulation is without instabilities and should be done either for a spectrum E^-1. The spectrum is reweighted to E^-alpha {1, exp(-E / Emax)} (if E < Emax or E > Emax).
The default output style from CRPropa supported is "all", with columns:
 D, z, SN, ID, E, X, Y, Z, Px, Py, Pz, SN0, ID0, E0, X0, Y0, Z0, P0x, P0y, P0z,  SN1, ID1, E1, X1, Y1, Z1, P1x, P1y, P1z, W
For now it is only possible to re-weight the simulation to include jets with opening angles smaller than the one used in the simulation. The user must ensure consistency themself.
If the jetAngle is smaller than 0, then no angular cuts are done.
 */


class EmissionGeometry : public crpropa::Referenced {
	public:
		EmissionGeometry();
		virtual double computeArea() const = 0;
		virtual double computeVolume() const = 0;
};

class Cone : public EmissionGeometry {
	private:
		crpropa::Vector3d axis;
		double angle;
		double height;
		double radius;
	public:
		Cone(double angle, double height, crpropa::Vector3d axis = crpropa::Vector3d(-1, 0, 0));
		void setAngle(double angle);
		void setHeight(double height);
		void setRadius(double radius);
		void setAxis(crpropa::Vector3d axis);
		double computeRadius() const;
		double computeArea() const;
		double computeVolume() const;
};


class Scenario : public crpropa::Referenced {
	private:
		double spectralIndex;
		double spectralIndexSimulation;
		double energyCutoff;
		double luminosity;
		double redshift;
		double distance;
		crpropa::ref_ptr<EmissionGeometry> geometry;
	public:
		Scenario(double spectralIndex, double energyCutoff, double luminosity, double distance, crpropa::ref_ptr<EmissionGeometry> geometry, double spectralIndexSimulation = 1, double redshift = 0);
		void setSpectralIndex(double s);
		void setSpectralIndexSimulation(double s);
		void setEnergyCutoff(double Emax);
		void setLuminosity(double L);
		void setDistance(double d);
		void setSourceRedshift(double z);
		void setGeometry(crpropa::ref_ptr<EmissionGeometry> geometry);
		double getSpectralIndex() const;
		double getSpectralIndexSimulation() const;
		double getEnergyCutoff() const;
		double getLuminosity() const;
		double getDistance() const;
		double getSourceRedshift() const;
		double getSourceLightTravelDistance() const;
		double getSourceComovingDistance() const;
		double getSourceLuminosityDistance() const;
		crpropa::ref_ptr<EmissionGeometry> getGeometry() const;
		double computeVolume() const;
};


class EmissionObservables : public crpropa::Referenced {
	private:
		double meanLorentzFactor;
		double meanInverseLorentzFactor;
		double weightSum;
		double energyTotal;
		double energyTotal0;
		double density;
	public:
		EmissionObservables();
		EmissionObservables(double lorentzFactor, double inverseLorentzFactor, double energyTotal, double energyTotal0, double weightSum, double density = 0.);
		void setMeanLorentzFactor(double lorentzFactor);
		void setMeanInverseLorentzFactor(double inverseLorentzFactor);
		void setWeightSum(double weightSum);
		void setEnergyTotal(double energyTotal);
		void setEnergyTotal0(double energyTotal0);
		void setDensity(double density);
		double getMeanLorentzFactor() const;
		double getMeanInverseLorentzFactor() const;
		double getWeightSum() const;
		double getEnergyTotal() const;
		double getEnergyTotal0() const;
		double getDensity() const;
		void incrementMeanLorentzFactor(double dLorentzFactor);
		void incrementMeanInverseLorentzFactor(double dInverseLorentzFactor);
		void incrementWeightSum(double dWeightSum);
		void incrementEnergyTotal(double dEnergyTotal);
		void incrementEnergyTotal0(double dEnergyTotal0);
};


class Simulation : public crpropa::Referenced {
	protected:
		std::string filename;
		std::vector<crpropa::ref_ptr<Scenario>> scenarios;
		std::vector<crpropa::ref_ptr<EmissionObservables>> observables;
	public:
		Simulation();
		Simulation(std::string filename);
		~Simulation();
		void setFilename(std::string filename);
		void addScenario(crpropa::ref_ptr<Scenario> scenario);
		void process();
		std::string getFilename() const;
		std::vector<crpropa::ref_ptr<Scenario>> getScenarios() const;
		std::vector<crpropa::ref_ptr<EmissionObservables>> getEmissionObservables() const;
};

void saveEmissionProfile(std::vector<crpropa::ref_ptr<Simulation>>, std::string filename);


}

#endif // GRPLINST_SIMULATION_H