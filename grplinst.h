#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Cosmology.h>
#include <crpropa/Random.h>

#ifdef _OPENMP
    #include "omp.h"
#endif

class PlasmaInstability : public crpropa::Module
{
	double temperatureIGM;
	double densityIGM;
	double luminosityBeam;
	char plasmaInstabilityModel;
	std::vector<double> _w;
	std::vector<double> _d;

public:
	/// The parent's constructor need to be called on initialization!
	PlasmaInstability(double densityIGM, double temperature, double beamLuminosity, char plasmaInstabilityModel);
	void process(crpropa::Candidate *candidate) const;
	void setModel(char model);
	void setIGMDensity(double density);
	void setBeamLuminosity(double density);
	void setTemperature(double temperature);
	void initTables();
	double coolingPower(double energy, double redshift) const;
	double coolingPowerA(double energy, double redshift) const;
	double coolingPowerB(double energy, double redshift) const;
	double coolingPowerC(double energy, double redshift) const;
	double coolingPowerD(double energy, double redshift) const;
	double coolingPowerE(double energy, double redshift) const;
};
