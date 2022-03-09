#include "grplinst/Flow.h"


namespace grplinst {


Flow::Flow(crpropa::Vector3d origin) {
}

Flow::~Flow() {
}

void Flow::setOrigin(crpropa::Vector3d centre) {
	origin = centre;
}

crpropa::Vector3d Flow::getOrigin() const {
	return origin;
}



FlowHomogeneous::FlowHomogeneous(double n, crpropa::Vector3d centre) {
	setOrigin(centre);
	setDensity(n);
}

FlowHomogeneous::~FlowHomogeneous() {
}

void FlowHomogeneous::setDensity(double n) {
	density = n;
}

double FlowHomogeneous::getDensity() const {
	return density;
}

double FlowHomogeneous::getValue(crpropa::Vector3d position, double redshift) const {
	return density * crpropa::pow_integer<3>(1 + redshift);
}



FlowJet1D::FlowJet1D(double n, crpropa::Vector3d centre) {
	setOrigin(centre);
	setDensity(n);
	initTable();
}

FlowJet1D::~FlowJet1D() {
}

void FlowJet1D::setDensity(double n) {
	// normalisation
	density = n;
}

double FlowJet1D::getDensity() const {
	return density;
}

double FlowJet1D::getValue(crpropa::Vector3d position, double redshift) const {
	double n = crpropa::interpolate((position - origin).getR(), _d, _n);
	// n = pow(10, n);
	return n * crpropa::pow_integer<3>(1 + redshift);
}

void FlowJet1D::initTable() {
	// the first element (0) is going to be inferred from extrapolation
	double d_[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
	double w_[] = {0, 2.65e-21, 4.17e-22, 1.19e-22, 4.36e-23, 1.87e-23, 8.94e-24, 4.67e-24, 2.63e-24, 1.58e-24};
	// double w_[] = {0, 2.65e-10, 4.17e-11, 1.19e-15, 4.36e-18, 1.87e-19, 8.94e-22, 4.67e-24, 2.63e-25, 1.58e-30};

	// // simple extrapolation
	// w_[0] = exp((log(w_[2] / w_[1]) / (d_[2] - d_[1])) * (d_[1] - d_[0]) - log(w_[1]));
	w_[0] = w_[1] - ((w_[2] - w_[1]) / (d_[2] - d_[1])) * (d_[1] - d_[0]);

	_d.insert(_d.begin(), std::begin(d_), std::end(d_));
	_n.insert(_n.begin(), std::begin(w_), std::end(w_));

	for (size_t i = 0; i < _d.size(); i++) {
		_d[i] *= crpropa::Mpc;
		_n[i] /= _n[0];
		_n[i] *= density;
		// _n[i] = log10(_n[i]);
	}
}

} // namespace grplinst