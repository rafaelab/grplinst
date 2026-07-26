#include <gtest/gtest.h>
#include <crpropa/Candidate.h>
#include <crpropa/Units.h>
#include "grplinst/PlasmaInstability.h"
#include "grplinst/Medium.h"
#include "grplinst/Flow.h"

using namespace grplinst;
using namespace crpropa;


namespace {

crpropa::Candidate makeElectron(double E = 1 * crpropa::TeV, double z = 0.) {
	crpropa::Candidate c;
	c.current.setId(11);
	c.current.setEnergy(E);
	c.current.setPosition(crpropa::Vector3d(0, 0, 0));
	c.setRedshift(z);
	return c;
}

crpropa::ref_ptr<FlowHomogeneous> makeFlow(double L = 1e38) {
	return new FlowHomogeneous(L, crpropa::Vector3d(0, 0, 0));
}

crpropa::ref_ptr<MediumDensityHomogeneous> makeDensity(double n = 0.1) {
	return new MediumDensityHomogeneous(n);
}

crpropa::ref_ptr<MediumTemperatureHomogeneous> makeTemperature(double T = 1e4) {
	return new MediumTemperatureHomogeneous(T);
}

} // anonymous namespace


// ─── PlasmaInstability base ──────────────────────────────────────────────────

TEST(PlasmaInstability, SetGetFlowProperties) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature());
	auto flow = makeFlow(2e38);
	pi.setFlowProperties(flow);
	EXPECT_EQ(pi.getFlowProperties().get(), flow.get());
}

TEST(PlasmaInstability, SetGetMediumDensity) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature());
	auto dens = makeDensity(0.5);
	pi.setMediumDensity(dens);
	EXPECT_EQ(pi.getMediumDensity().get(), dens.get());
}

TEST(PlasmaInstability, SetGetMediumTemperature) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature());
	auto temp = makeTemperature(2e4);
	pi.setMediumTemperature(temp);
	EXPECT_EQ(pi.getMediumTemperature().get(), temp.get());
}

TEST(PlasmaInstability, SetGetLimit) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature());
	pi.setLimit(0.05);
	EXPECT_DOUBLE_EQ(pi.getLimit(), 0.05);
}

TEST(PlasmaInstability, EfficiencyClampedToZeroWhenNegative) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature(), /*efficiency=*/-0.5);
	EXPECT_DOUBLE_EQ(pi.getEfficiencyFactor(), 0.);
}

TEST(PlasmaInstability, EfficiencyClampedToOneWhenExceedsOne) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature(), /*efficiency=*/1.5);
	EXPECT_DOUBLE_EQ(pi.getEfficiencyFactor(), 1.);
}

TEST(PlasmaInstability, EfficiencySetCorrectlyInRange) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature(), /*efficiency=*/0.7);
	EXPECT_DOUBLE_EQ(pi.getEfficiencyFactor(), 0.7);
}

TEST(PlasmaInstability, ComputeEnergyLossPerLengthIsPositive) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature());
	auto c = makeElectron();
	EXPECT_GT(pi.computeEnergyLossPerLength(c), 0.);
}

TEST(PlasmaInstability, ProcessReducesElectronEnergy) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature());
	crpropa::ref_ptr<crpropa::Candidate> c = new crpropa::Candidate();
	c->current.setId(11);
	c->current.setEnergy(1 * TeV);
	c->current.setPosition(Vector3d(0, 0, 0));
	c->setRedshift(0.);
	c->setCurrentStep(1 * kpc);
	pi.process(c);
	EXPECT_LE(c->current.getEnergy(), 1 * TeV);
}

TEST(PlasmaInstability, ProcessIgnoresNonElectrons) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature());
	crpropa::ref_ptr<crpropa::Candidate> c = new crpropa::Candidate();
	c->current.setId(22);  // photon
	double E0 = 1 * TeV;
	c->current.setEnergy(E0);
	c->current.setPosition(Vector3d(0, 0, 0));
	c->setRedshift(0.);
	c->setCurrentStep(1 * kpc);
	pi.process(c);
	EXPECT_DOUBLE_EQ(c->current.getEnergy(), E0);
}


// ─── energyLossTime for each model ──────────────────────────────────────────

TEST(PlasmaInstabilityBroderick2012, EnergyLossTimeIsPositive) {
	PlasmaInstabilityBroderick2012 pi(makeFlow(), makeDensity(), makeTemperature());
	auto c = makeElectron();
	EXPECT_GT(pi.energyLossTime(c), 0.);
}

TEST(PlasmaInstabilityBroderick2012, EnergyLossTimeDecreasesWithBeamDensity) {
	// Higher beam density → shorter loss time (both regimes)
	PlasmaInstabilityBroderick2012 pi1(makeFlow(1e38), makeDensity(), makeTemperature());
	PlasmaInstabilityBroderick2012 pi2(makeFlow(1e40), makeDensity(), makeTemperature());
	auto c = makeElectron();
	EXPECT_LT(pi2.energyLossTime(c), pi1.energyLossTime(c));
}

TEST(PlasmaInstabilitySchlickeiser2012, EnergyLossTimeIsPositive) {
	PlasmaInstabilitySchlickeiser2012 pi(makeFlow(), makeDensity(), makeTemperature());
	auto c = makeElectron();
	EXPECT_GT(pi.energyLossTime(c), 0.);
}

TEST(PlasmaInstabilitySironi2014, EnergyLossTimeIsPositive) {
	PlasmaInstabilitySironi2014 pi(makeFlow(), makeDensity(), makeTemperature());
	auto c = makeElectron();
	EXPECT_GT(pi.energyLossTime(c), 0.);
}

TEST(PlasmaInstabilityVafin2018, EnergyLossTimeIsPositive) {
	PlasmaInstabilityVafin2018 pi(makeFlow(), makeDensity(), makeTemperature());
	auto c = makeElectron();
	EXPECT_GT(pi.energyLossTime(c), 0.);
}

TEST(PlasmaInstabilityVafin2018, EnergyLossTimeDecreasesWithTemperature) {
	// tau ~ 1/T, so higher temperature → shorter loss time
	PlasmaInstabilityVafin2018 pi1(makeFlow(), makeDensity(), makeTemperature(1e4));
	PlasmaInstabilityVafin2018 pi2(makeFlow(), makeDensity(), makeTemperature(2e4));
	auto c = makeElectron();
	EXPECT_NEAR(pi1.energyLossTime(c) / pi2.energyLossTime(c), 2.0, 1e-10);
}

TEST(PlasmaInstabilityBret2010TwoStream, EnergyLossTimeIsPositive) {
	PlasmaInstabilityBret2010TwoStream pi(makeFlow(), makeDensity(), makeTemperature());
	auto c = makeElectron();
	EXPECT_GT(pi.energyLossTime(c), 0.);
}

TEST(PlasmaInstabilityBret2010Filamentation, EnergyLossTimeIsPositive) {
	PlasmaInstabilityBret2010Filamentation pi(makeFlow(), makeDensity(), makeTemperature());
	auto c = makeElectron();
	EXPECT_GT(pi.energyLossTime(c), 0.);
}

TEST(PlasmaInstabilityShalaby2020, EnergyLossTimeIsPositive) {
	PlasmaInstabilityShalaby2020 pi(makeFlow(), makeDensity(), makeTemperature());
	auto c = makeElectron();
	EXPECT_GT(pi.energyLossTime(c), 0.);
}


// ─── PlasmaInstabilityMiniati2013 ────────────────────────────────────────────

TEST(PlasmaInstabilityMiniati2013, DefaultConstructor) {
	EXPECT_NO_THROW(PlasmaInstabilityMiniati2013());
}

TEST(PlasmaInstabilityMiniati2013, ParameterizedConstructor) {
	EXPECT_NO_THROW(PlasmaInstabilityMiniati2013(1e38, makeDensity(), makeTemperature()));
}

TEST(PlasmaInstabilityMiniati2013, GettersAfterConstruction) {
	auto dens = makeDensity();
	auto temp = makeTemperature();
	PlasmaInstabilityMiniati2013 pi(1e38, dens, temp);
	EXPECT_EQ(pi.getMediumDensity().get(), dens.get());
	EXPECT_EQ(pi.getMediumTemperature().get(), temp.get());
	EXPECT_NE(pi.getFlowProperties().get(), nullptr);
}

TEST(PlasmaInstabilityMiniati2013, EnergyLossTimeDoesNotThrow) {
	PlasmaInstabilityMiniati2013 pi(1e38, makeDensity(), makeTemperature());
	// Place the candidate at a non-zero position to avoid log10(0)
	crpropa::Candidate c;
	c.current.setId(11);
	c.current.setEnergy(1 * TeV);
	c.current.setPosition(Vector3d(10 * Mpc, 0, 0));
	c.setRedshift(0.);
	EXPECT_NO_THROW(pi.energyLossTime(c));
}


// ─── Helper functions ────────────────────────────────────────────────────────

TEST(PlasmaFrequency, IsPositiveForElectron) {
	double wp = plasmaFrequency(1e-1, 11);
	EXPECT_GT(wp, 0.);
}

TEST(PlasmaFrequency, ScalesWithSqrtDensity) {
	double wp1 = plasmaFrequency(1e-1, 11);
	double wp4 = plasmaFrequency(4e-1, 11);
	// omega_p ~ sqrt(n), so quadrupling n doubles omega_p
	EXPECT_NEAR(wp4 / wp1, 2.0, 1e-10);
}

TEST(MaximumLinearGrowthFrequency, IsPositive) {
	double gamma = maximumLinearGrowthFrequency(1e-16, 1e-1, 1000.);
	EXPECT_GT(gamma, 0.);
}

TEST(MaximumLinearGrowthFrequency, LinearInBeamDensity) {
	double g1 = maximumLinearGrowthFrequency(1e-16, 1e-1, 1000.);
	double g2 = maximumLinearGrowthFrequency(2e-16, 1e-1, 1000.);
	EXPECT_NEAR(g2 / g1, 2.0, 1e-10);
}

TEST(MaximumLinearGrowthFrequency, InverseInLorentzFactor) {
	double g1 = maximumLinearGrowthFrequency(1e-16, 1e-1, 1000.);
	double g2 = maximumLinearGrowthFrequency(1e-16, 1e-1, 2000.);
	EXPECT_NEAR(g1 / g2, 2.0, 1e-10);
}
