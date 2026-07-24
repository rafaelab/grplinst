#include <gtest/gtest.h>
#include <crpropa/Units.h>
#include "grplinst/Flow.h"


using namespace grplinst;
using namespace crpropa;


///////////////////////////////////////////////////////////////////////////////
//                             FlowHomogeneous                               //
///////////////////////////////////////////////////////////////////////////////

TEST(FlowHomogeneous, DefaultConstructor) {
	EXPECT_NO_THROW(FlowHomogeneous());
}

TEST(FlowHomogeneous, ParameterizedConstructor) {
	double L = 1e38;
	Vector3d origin(1 * Mpc, 0, 0);
	FlowHomogeneous flow(L, origin);
	EXPECT_DOUBLE_EQ(flow.getLuminosity(), L);
	EXPECT_EQ(flow.getOrigin(), origin);
}

TEST(FlowHomogeneous, SetGetLuminosity) {
	FlowHomogeneous flow;
	flow.setLuminosity(2e40);
	EXPECT_DOUBLE_EQ(flow.getLuminosity(), 2e40);
}

TEST(FlowHomogeneous, SetGetOrigin) {
	FlowHomogeneous flow;
	Vector3d origin(5 * Mpc, 3 * Mpc, 0);
	flow.setOrigin(origin);
	EXPECT_EQ(flow.getOrigin(), origin);
}

TEST(FlowHomogeneous, GetDensityIsPositive) {
	FlowHomogeneous flow(1e38);
	double d = flow.getDensity(1 * TeV, Vector3d(0, 0, 0), 0.);
	EXPECT_GT(d, 0.);
}

TEST(FlowHomogeneous, GetDensityScalesWithEnergy) {
	FlowHomogeneous flow(1e38);
	Vector3d pos(0, 0, 0);
	double d1 = flow.getDensity(1 * TeV, pos, 0.);
	double d2 = flow.getDensity(10 * TeV, pos, 0.);
	EXPECT_NE(d1, d2);
}

TEST(FlowHomogeneous, GetDensityChangesWithRedshift) {
	FlowHomogeneous flow(1e38);
	Vector3d pos(0, 0, 0);
	double d0 = flow.getDensity(1 * TeV, pos, 0.0);
	double dz = flow.getDensity(1 * TeV, pos, 0.5);
	EXPECT_NE(d0, dz);
}

TEST(FlowHomogeneous, GetDensityScalesWithLuminosity) {
	FlowHomogeneous flow1(1e38);
	FlowHomogeneous flow2(2e38);
	Vector3d pos(0, 0, 0);
	double d1 = flow1.getDensity(1 * TeV, pos, 0.);
	double d2 = flow2.getDensity(1 * TeV, pos, 0.);
	EXPECT_NEAR(d2 / d1, 2.0, 1e-10);
}


///////////////////////////////////////////////////////////////////////////////
//                                FlowJet1D                                  //
///////////////////////////////////////////////////////////////////////////////

TEST(FlowJet1D, DefaultConstructor) {
	EXPECT_NO_THROW(FlowJet1D());
}

TEST(FlowJet1D, ParameterizedConstructor) {
	std::vector<double> dist = {1 * Mpc, 2 * Mpc, 4 * Mpc};
	std::vector<double> dens = {1e-3, 5e-4, 1e-4};
	FlowJet1D flow(dist, dens, 1e38);
	EXPECT_DOUBLE_EQ(flow.getLuminosity(), 1e38);
	EXPECT_EQ(flow.getDistanceProfile().size(), 3u);
	EXPECT_EQ(flow.getDensityProfile().size(), 3u);
}

TEST(FlowJet1D, MismatchedVectorsThrows) {
	std::vector<double> dist = {1.0, 2.0, 3.0};
	std::vector<double> dens = {1e-3, 5e-4};  // one fewer element
	EXPECT_THROW(FlowJet1D(dist, dens), std::invalid_argument);
}

TEST(FlowJet1D, SetGetDistanceProfile) {
	FlowJet1D flow;
	std::vector<double> dist = {0.0, 1.0, 2.0};
	flow.setDistanceProfile(dist);
	EXPECT_EQ(flow.getDistanceProfile(), dist);
}

TEST(FlowJet1D, SetGetDensityProfile) {
	FlowJet1D flow;
	std::vector<double> dens = {1e-3, 5e-4, 1e-5};
	flow.setDensityProfile(dens);
	EXPECT_EQ(flow.getDensityProfile(), dens);
}

TEST(FlowJet1D, GetDensityLinearInterpolation) {
	// Linear profile from 4 down to 1, sampled at integer distances
	std::vector<double> dist = {0.0, 1.0, 2.0, 3.0};
	std::vector<double> dens = {4.0, 3.0, 2.0, 1.0};
	FlowJet1D flow(dist, dens, 1.0, Vector3d(0, 0, 0), /*interpolateLog=*/false);

	// At |pos| = 1.5: density should interpolate to 2.5
	double d = flow.getDensity(1 * TeV, Vector3d(1.5, 0, 0), 0.);
	EXPECT_NEAR(d, 2.5, 1e-10);
}

TEST(FlowJet1D, GetDensityRedshiftScaling) {
	std::vector<double> dist = {0.0, 10.0};
	std::vector<double> dens = {2.0, 1.0};
	FlowJet1D flow(dist, dens, 1.0, Vector3d(0, 0, 0), /*interpolateLog=*/false);

	double d0 = flow.getDensity(1 * TeV, Vector3d(5.0, 0, 0), 0.);
	double dz = flow.getDensity(1 * TeV, Vector3d(5.0, 0, 0), 1.);
	// density scales as (1 + z)^3: at z=1, factor = 2^3 = 8
	EXPECT_NEAR(dz / d0, 8.0, 1e-10);
}

TEST(FlowJet1D, GetDensityUsesOriginOffset) {
	std::vector<double> dist = {0.0, 10.0};
	std::vector<double> dens = {2.0, 1.0};
	Vector3d origin(3.0, 0, 0);
	FlowJet1D flow(dist, dens, 1.0, origin, false);

	// Position (8, 0, 0): distance from origin = 5 → density = 1.5
	double d = flow.getDensity(1 * TeV, Vector3d(8.0, 0, 0), 0.);
	EXPECT_NEAR(d, 1.5, 1e-10);
}


///////////////////////////////////////////////////////////////////////////////
//                       FlowJet1D - Miniati 2013 model                      //
///////////////////////////////////////////////////////////////////////////////

TEST(FlowMiniati2013, ReturnsNonNull) {
	auto flow = createFlowMiniati2013(1e38);
	EXPECT_NE(flow.get(), nullptr);
}

TEST(FlowMiniati2013, ProfileHasSixteenPoints) {
	auto flow = createFlowMiniati2013(1e38);
	FlowJet1D* jet = dynamic_cast<FlowJet1D*>(flow.get());
	ASSERT_NE(jet, nullptr);
	EXPECT_EQ(jet->getDistanceProfile().size(), 16u);
	EXPECT_EQ(jet->getDensityProfile().size(), 16u);
}

TEST(FlowMiniati2013, DensityProfileScalesWithLuminosity) {
	auto flow1 = createFlowMiniati2013(1e38);
	auto flow2 = createFlowMiniati2013(2e38);

	FlowJet1D* jet1 = dynamic_cast<FlowJet1D*>(flow1.get());
	FlowJet1D* jet2 = dynamic_cast<FlowJet1D*>(flow2.get());
	ASSERT_NE(jet1, nullptr);
	ASSERT_NE(jet2, nullptr);

	// density values scale linearly with luminosity
	EXPECT_NEAR(jet2->getDensityProfile()[0] / jet1->getDensityProfile()[0], 2.0, 1e-10);
}

TEST(FlowMiniati2013, LuminosityIsSet) {
	auto flow = createFlowMiniati2013(1e38);
	EXPECT_DOUBLE_EQ(flow->getLuminosity(), 1e38);
}
