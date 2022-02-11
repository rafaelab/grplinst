#ifndef GRPLINST_FLOW_H
#define GRPLINST_FLOW_H

#include <crpropa/Units.h>
#include <crpropa/Vector3.h>
#include <crpropa/Grid.h>
#include <crpropa/Referenced.h>


namespace grplinst {

/**
 @class Beam
 @brief Abstract base class to define properties related to the emitting object.
 */
class Flow : public crpropa::Referenced {
	public:
		Flow() {
		}
		virtual double getValue(crpropa::Vector3d position) const = 0;
		virtual double getValue(crpropa::Vector3d position, double redshift) const = 0;
};

/**
 @class HomogeneousFlow
 @brief Class to describe the  
 */
class HomogeneousFlow : public Flow {
	private:
		double density;
	public:
		HomogeneousFlow(double density);
		void setDensity(double density);
		double getDensity() const;
		double getValue(crpropa::Vector3d position) const;
		double getValue(crpropa::Vector3d position, double redshift) const;
};

/**
Converts a given density to luminosity and vice-versa.
This is essentially eq. 13 from:
  R. Alves Batista, A. Saveliev, E. M. de Gouveia Dal Pino. 
  Mon. Not. R. Astron. Soc. 489 (2019) 3836-3849.
  arXiv:1904.13345
which is based on:
  Broderick et al. 
  Astrophys. J. 752 (2012) 22. 
  arXiv:1106.5494
This is an upper limit based on analytical estimates of cascade development assuming inverse Compton,
and considering the plasma cooling dominated by the kinetic oblique mode.
 */

inline double beamDensityToLuminosity(double n, double Ee, double z, double a = 4.5) {
	return 1e38 * (n / 3.7e-16) * (crpropa::TeV / Ee) * pow((1. + z) / 2., 4. - 3. * a);
}

inline double beamLuminosityToDensity(double L, double Ee, double z, double a = 4.5) {
	return 3.7e-16 * (L / 1e38) * (Ee / crpropa::TeV) * pow((1. + z) / 2., 3. * a - 4.);
}

} // namespace grplinst

#endif // GRPLINST_FLOW_H