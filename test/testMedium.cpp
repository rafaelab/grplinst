#include <gtest/gtest.h>
#include <crpropa/Units.h>
#include "grplinst/Medium.h"

using namespace grplinst;
using namespace crpropa;


// ─── MediumTemperatureHomogeneous ────────────────────────────────────────────

TEST(MediumTemperatureHomogeneous, Constructor) {
	EXPECT_NO_THROW(MediumTemperatureHomogeneous(1e4));
}

TEST(MediumTemperatureHomogeneous, GetSetTemperatureValue) {
	MediumTemperatureHomogeneous med(1e4);
	EXPECT_DOUBLE_EQ(med.getTemperatureValue(), 1e4);
	med.setTemperatureValue(2e4);
	EXPECT_DOUBLE_EQ(med.getTemperatureValue(), 2e4);
}

TEST(MediumTemperatureHomogeneous, GetTemperatureAtZeroRedshift) {
	MediumTemperatureHomogeneous med(1e4);
	double T = med.getTemperature(Vector3d(0, 0, 0), 0.);
	EXPECT_DOUBLE_EQ(T, 1e4);
}

TEST(MediumTemperatureHomogeneous, GetTemperatureScalesWithRedshift) {
	MediumTemperatureHomogeneous med(1e4);
	double T = med.getTemperature(Vector3d(0, 0, 0), 1.);
	EXPECT_DOUBLE_EQ(T, 2e4);  // T * (1 + z) = 1e4 * 2
}

TEST(MediumTemperatureHomogeneous, GetTemperatureIsPositionIndependent) {
	MediumTemperatureHomogeneous med(1e4);
	double T1 = med.getTemperature(Vector3d(0, 0, 0), 0.);
	double T2 = med.getTemperature(Vector3d(100 * Mpc, 50 * Mpc, 0), 0.);
	EXPECT_DOUBLE_EQ(T1, T2);
}

TEST(MediumTemperatureHomogeneous, GetVelocityForElectronIsPositive) {
	MediumTemperatureHomogeneous med(1e4);
	double v = med.getVelocity(11, Vector3d(0, 0, 0), 0.);
	EXPECT_GT(v, 0.);
}

TEST(MediumTemperatureHomogeneous, GetVelocityScalesWithSqrtTemperature) {
	// v ~ sqrt(k_B * T / m), so doubling T multiplies v by sqrt(2)
	MediumTemperatureHomogeneous med1(1e4);
	MediumTemperatureHomogeneous med2(4e4);
	double v1 = med1.getVelocity(11, Vector3d(0, 0, 0), 0.);
	double v2 = med2.getVelocity(11, Vector3d(0, 0, 0), 0.);
	EXPECT_NEAR(v2 / v1, 2.0, 1e-10);  // sqrt(4) = 2
}

TEST(MediumTemperatureHomogeneous, GetVelocityDependsOnRedshift) {
	// v depends on T(z) = T0 * (1 + z), so v ~ sqrt(1 + z)
	MediumTemperatureHomogeneous med(1e4);
	double v0 = med.getVelocity(11, Vector3d(0, 0, 0), 0.);
	double v1 = med.getVelocity(11, Vector3d(0, 0, 0), 1.);
	EXPECT_NEAR(v1 / v0, std::sqrt(2.0), 1e-10);  // sqrt((1+1)/(1+0))
}


// ─── MediumDensityHomogeneous ─────────────────────────────────────────────────

TEST(MediumDensityHomogeneous, Constructor) {
	EXPECT_NO_THROW(MediumDensityHomogeneous(1e-1));
}

TEST(MediumDensityHomogeneous, GetSetDensityValue) {
	MediumDensityHomogeneous med(1e-1);
	EXPECT_DOUBLE_EQ(med.getDensityValue(), 1e-1);
	med.setDensityValue(5e-2);
	EXPECT_DOUBLE_EQ(med.getDensityValue(), 5e-2);
}

TEST(MediumDensityHomogeneous, GetDensityAtZeroRedshift) {
	MediumDensityHomogeneous med(1e-1);
	double n = med.getDensity(Vector3d(0, 0, 0), 0.);
	EXPECT_DOUBLE_EQ(n, 1e-1);
}

TEST(MediumDensityHomogeneous, GetDensityScalesWithRedshift) {
	MediumDensityHomogeneous med(1e-1);
	// n(z) = n0 * (1 + z)^3, at z=1: 0.1 * 2^3 = 0.8
	double n = med.getDensity(Vector3d(0, 0, 0), 1.);
	EXPECT_DOUBLE_EQ(n, 0.8);
}

TEST(MediumDensityHomogeneous, GetDensityIsPositionIndependent) {
	MediumDensityHomogeneous med(1e-1);
	double n1 = med.getDensity(Vector3d(0, 0, 0), 0.);
	double n2 = med.getDensity(Vector3d(100 * Mpc, 50 * Mpc, 0), 0.);
	EXPECT_DOUBLE_EQ(n1, n2);
}

TEST(MediumDensityHomogeneous, GetDensityIsPositive) {
	MediumDensityHomogeneous med(1e-1);
	EXPECT_GT(med.getDensity(Vector3d(0, 0, 0), 0.), 0.);
	EXPECT_GT(med.getDensity(Vector3d(0, 0, 0), 1.), 0.);
}
