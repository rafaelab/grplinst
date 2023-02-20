#ifndef GRPLINST_H
#define GRPLINST_H


#include <crpropa/Cosmology.h>
#include <crpropa/Module.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Units.h>




class PlasmaInstability : public crpropa::Module {
	protected:
		double temperatureIGM;
		double densityIGM;
		double luminosityBeam;
		char plasmaInstabilityModel;
		std::vector<double> _w;
		std::vector<double> _d;

	public:
		PlasmaInstability(double densityIGM, double temperature, double beamLuminosity, char plasmaInstabilityModel);
		void initTables();
		void setModel(char model);
		void setIGMDensity(double density);
		void setBeamLuminosity(double density);
		void setTemperature(double temperature);
		double getIGMDensity() const;
		double getBeamLuminosity() const;
		double getTemperature() const;
		char getModel() const;
		double coolingPower(double energy, double redshift) const;
		double coolingPowerA(double energy, double redshift) const;
		double coolingPowerB(double energy, double redshift) const;
		double coolingPowerC(double energy, double redshift) const;
		double coolingPowerD(double energy, double redshift) const;
		double coolingPowerE(double energy, double redshift) const;
		std::string getModelReference() const;
		std::string getDescription() const;
		void process(crpropa::Candidate *candidate) const;
};


#endif