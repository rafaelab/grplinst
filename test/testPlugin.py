import sys
import math
import unittest

from crpropa import *
from grplinst import *


def makeElectron(E = 1 * TeV, z = 0.):
    c = Candidate()
    c.current.setId(11)
    c.current.setEnergy(E)
    c.current.setPosition(Vector3d(0, 0, 0))
    c.setRedshift(z)
    return c


class TestFlowHomogeneous(unittest.TestCase):

    def testDefaultConstructor(self):
        FlowHomogeneous()

    def testParameterizedConstructor(self):
        flow = FlowHomogeneous(1e38, Vector3d(0, 0, 0))
        self.assertAlmostEqual(flow.getLuminosity(), 1e38)

    def testSetGetLuminosity(self):
        flow = FlowHomogeneous()
        flow.setLuminosity(2e40)
        self.assertAlmostEqual(flow.getLuminosity(), 2e40)

    def testSetGetOrigin(self):
        flow = FlowHomogeneous()
        origin = Vector3d(5 * Mpc, 0, 0)
        flow.setOrigin(origin)
        self.assertEqual(flow.getOrigin(), origin)

    def testGetDensityPositive(self):
        flow = FlowHomogeneous(1e38)
        d = flow.getDensity(1 * TeV, Vector3d(0, 0, 0), 0.)
        self.assertGreater(d, 0.)

    def testGetDensityScalesWithLuminosity(self):
        flow1 = FlowHomogeneous(1e38)
        flow2 = FlowHomogeneous(2e38)
        d1 = flow1.getDensity(1 * TeV, Vector3d(0, 0, 0), 0.)
        d2 = flow2.getDensity(1 * TeV, Vector3d(0, 0, 0), 0.)
        self.assertAlmostEqual(d2 / d1, 2.0, places = 10)

    def testGetDensityChangesWithRedshift(self):
        flow = FlowHomogeneous(1e38)
        d0 = flow.getDensity(1 * TeV, Vector3d(0, 0, 0), 0.0)
        dz = flow.getDensity(1 * TeV, Vector3d(0, 0, 0), 0.5)
        self.assertNotAlmostEqual(d0, dz)


class TestFlowJet1D(unittest.TestCase):

    def testDefaultConstructor(self):
        FlowJet1D()

    def testParameterizedConstructor(self):
        dist = [1 * Mpc, 2 * Mpc, 4 * Mpc]
        dens = [1e-3, 5e-4, 1e-4]
        flow = FlowJet1D(dist, dens, 1e38)
        self.assertAlmostEqual(flow.getLuminosity(), 1e38)
        self.assertEqual(len(flow.getDistanceProfile()), 3)
        self.assertEqual(len(flow.getDensityProfile()), 3)

    def testMismatchedVectorsRaises(self):
        with self.assertRaises(Exception):
            FlowJet1D([1., 2., 3.], [1e-3, 5e-4])

    def testSetGetProfiles(self):
        flow = FlowJet1D()
        dist = [0.0, 1.0, 2.0]
        dens = [3.0, 2.0, 1.0]
        flow.setDistanceProfile(dist)
        flow.setDensityProfile(dens)
        self.assertEqual(list(flow.getDistanceProfile()), dist)
        self.assertEqual(list(flow.getDensityProfile()), dens)

    def testLinearInterpolation(self):
        dist = [0.0, 1.0, 2.0, 3.0]
        dens = [4.0, 3.0, 2.0, 1.0]
        flow = FlowJet1D(dist, dens, 1.0, Vector3d(0, 0, 0), False)
        d = flow.getDensity(1 * TeV, Vector3d(1.5, 0, 0), 0.)
        self.assertAlmostEqual(d, 2.5, places = 10)

    def testRedshiftScaling(self):
        dist = [0.0, 10.0]
        dens = [2.0, 1.0]
        flow = FlowJet1D(dist, dens, 1.0, Vector3d(0, 0, 0), False)
        d0 = flow.getDensity(1 * TeV, Vector3d(5.0, 0, 0), 0.0)
        dz = flow.getDensity(1 * TeV, Vector3d(5.0, 0, 0), 1.0)
        self.assertAlmostEqual(dz / d0, 8.0, places = 10)


class TestMediumTemperatureHomogeneous(unittest.TestCase):

    def testConstructor(self):
        MediumTemperatureHomogeneous(1e4)

    def testGetSetTemperatureValue(self):
        med = MediumTemperatureHomogeneous(1e4)
        self.assertAlmostEqual(med.getTemperatureValue(), 1e4)
        med.setTemperatureValue(2e4)
        self.assertAlmostEqual(med.getTemperatureValue(), 2e4)

    def testGetTemperatureAtZeroRedshift(self):
        med = MediumTemperatureHomogeneous(1e4)
        self.assertAlmostEqual(med.getTemperature(Vector3d(0, 0, 0), 0.), 1e4)

    def testGetTemperatureScalesWithRedshift(self):
        med = MediumTemperatureHomogeneous(1e4)
        self.assertAlmostEqual(med.getTemperature(Vector3d(0, 0, 0), 1.), 2e4)

    def testGetTemperatureIsPositionIndependent(self):
        med = MediumTemperatureHomogeneous(1e4)
        T1 = med.getTemperature(Vector3d(0, 0, 0), 0.)
        T2 = med.getTemperature(Vector3d(100 * Mpc, 0, 0), 0.)
        self.assertAlmostEqual(T1, T2)

    def testGetVelocityForElectronIsPositive(self):
        med = MediumTemperatureHomogeneous(1e4)
        v = med.getVelocity(11, Vector3d(0, 0, 0), 0.)
        self.assertGreater(v, 0.)

    def testGetVelocityScalesWithSqrtTemperature(self):
        med1 = MediumTemperatureHomogeneous(1e4)
        med2 = MediumTemperatureHomogeneous(4e4)
        v1 = med1.getVelocity(11, Vector3d(0, 0, 0), 0.)
        v2 = med2.getVelocity(11, Vector3d(0, 0, 0), 0.)
        self.assertAlmostEqual(v2 / v1, 2.0, places = 10)  # sqrt(4) = 2


class TestMediumDensityHomogeneous(unittest.TestCase):

    def testConstructor(self):
        MediumDensityHomogeneous(1e-1)

    def testGetSetDensityValue(self):
        med = MediumDensityHomogeneous(1e-1)
        self.assertAlmostEqual(med.getDensityValue(), 1e-1)
        med.setDensityValue(5e-2)
        self.assertAlmostEqual(med.getDensityValue(), 5e-2)

    def testGetDensityAtZeroRedshift(self):
        med = MediumDensityHomogeneous(1e-1)
        self.assertAlmostEqual(med.getDensity(Vector3d(0, 0, 0), 0.), 1e-1)

    def testGetDensityScalesWithRedshift(self):
        med = MediumDensityHomogeneous(1e-1)
        n = med.getDensity(Vector3d(0, 0, 0), 1.)
        self.assertAlmostEqual(n, 0.8)  # 0.1 * (1+1)^3 = 0.8

    def testGetDensityIsPositionIndependent(self):
        med = MediumDensityHomogeneous(1e-1)
        n1 = med.getDensity(Vector3d(0, 0, 0), 0.)
        n2 = med.getDensity(Vector3d(100 * Mpc, 0, 0), 0.)
        self.assertAlmostEqual(n1, n2)


class TestPlasmaInstability(unittest.TestCase):

    def setUp(self):
        self.flow = FlowHomogeneous(1e38, Vector3d(0, 0, 0))
        self.density = MediumDensityHomogeneous(0.1)
        self.temperature = MediumTemperatureHomogeneous(1e4)

    def testEfficiencyClampedWhenNegative(self):
        pi = PlasmaInstabilityBroderick2012(self.flow, self.density, self.temperature, -0.5)
        self.assertAlmostEqual(pi.getEfficiencyFactor(), 0.)

    def testEfficiencyClampedWhenAboveOne(self):
        pi = PlasmaInstabilityBroderick2012(self.flow, self.density, self.temperature, 1.5)
        self.assertAlmostEqual(pi.getEfficiencyFactor(), 1.)

    def testEfficiencySetCorrectlyInRange(self):
        pi = PlasmaInstabilityBroderick2012(self.flow, self.density, self.temperature, 0.7)
        self.assertAlmostEqual(pi.getEfficiencyFactor(), 0.7)

    def testSetGetLimit(self):
        pi = PlasmaInstabilityBroderick2012(self.flow, self.density, self.temperature)
        pi.setLimit(0.05)
        self.assertAlmostEqual(pi.getLimit(), 0.05)

    def testBroderick2012EnergyLossTimePositive(self):
        pi = PlasmaInstabilityBroderick2012(self.flow, self.density, self.temperature)
        self.assertGreater(pi.energyLossTime(make_electron()), 0.)

    def testSchlickeiser2012EnergyLossTimePositive(self):
        pi = PlasmaInstabilitySchlickeiser2012(self.flow, self.density, self.temperature)
        self.assertGreater(pi.energyLossTime(make_electron()), 0.)

    def testSironi2014EnergyLossTimePositive(self):
        pi = PlasmaInstabilitySironi2014(self.flow, self.density, self.temperature)
        self.assertGreater(pi.energyLossTime(make_electron()), 0.)

    def testVafin2018EnergyLossTimeInverselyProportionalToTemperature(self):
        pi1 = PlasmaInstabilityVafin2018(self.flow, self.density, MediumTemperatureHomogeneous(1e4))
        pi2 = PlasmaInstabilityVafin2018(self.flow, self.density, MediumTemperatureHomogeneous(2e4))
        c = make_electron()
        self.assertAlmostEqual(pi1.energyLossTime(c) / pi2.energyLossTime(c), 2.0, places = 10)

    def testBret2010TwoStreamEnergyLossTimePositive(self):
        pi = PlasmaInstabilityBret2010TwoStream(self.flow, self.density, self.temperature)
        self.assertGreater(pi.energyLossTime(make_electron()), 0.)

    def testBret2010FilamentationEnergyLossTimePositive(self):
        pi = PlasmaInstabilityBret2010Filamentation(self.flow, self.density, self.temperature)
        self.assertGreater(pi.energyLossTime(make_electron()), 0.)

    def testShalaby2020EnergyLossTimePositive(self):
        pi = PlasmaInstabilityShalaby2020(self.flow, self.density, self.temperature)
        self.assertGreater(pi.energyLossTime(make_electron()), 0.)

    def testMiniati2013DefaultConstructor(self):
        PlasmaInstabilityMiniati2013()

    def testMiniati2013ParameterizedConstructor(self):
        pi = PlasmaInstabilityMiniati2013(1e38, self.density, self.temperature)
        self.assertIsNotNone(pi.getFlowProperties())
        self.assertEqual(pi.getMediumDensity(), self.density)
        self.assertEqual(pi.getMediumTemperature(), self.temperature)

    def testMiniati2013EnergyLossTimeDoesNotThrow(self):
        pi = PlasmaInstabilityMiniati2013(1e38, self.density, self.temperature)
        c = Candidate()
        c.current.setId(11)
        c.current.setEnergy(1 * TeV)
        c.current.setPosition(Vector3d(10 * Mpc, 0, 0))
        c.setRedshift(0.)
        pi.energyLossTime(c)  # should not raise

    def testComputeEnergyLossPerLengthPositive(self):
        pi = PlasmaInstabilityBroderick2012(self.flow, self.density, self.temperature)
        self.assertGreater(pi.computeEnergyLossPerLength(make_electron()), 0.)

    def testProcessReducesElectronEnergy(self):
        pi = PlasmaInstabilityBroderick2012(self.flow, self.density, self.temperature)
        c = Candidate()
        c.current.setId(11)
        c.current.setEnergy(1 * TeV)
        c.current.setPosition(Vector3d(0, 0, 0))
        c.setRedshift(0.)
        c.setCurrentStep(1 * kpc)
        pi.process(c)
        self.assertLessEqual(c.current.getEnergy(), 1 * TeV)

    def testProcessIgnoresPhotons(self):
        pi = PlasmaInstabilityBroderick2012(self.flow, self.density, self.temperature)
        c = Candidate()
        c.current.setId(22)  # photon
        E0 = 1 * TeV
        c.current.setEnergy(E0)
        c.current.setPosition(Vector3d(0, 0, 0))
        c.setRedshift(0.)
        c.setCurrentStep(1 * kpc)
        pi.process(c)
        self.assertAlmostEqual(c.current.getEnergy(), E0)


class TestHelperFunctions(unittest.TestCase):

    def testPlasmaFrequencyPositive(self):
        wp = plasmaFrequency(0.1, 11)
        self.assertGreater(wp, 0.)

    def testPlasmaFrequencyScalesWithSqrtDensity(self):
        wp1 = plasmaFrequency(1e-1, 11)
        wp4 = plasmaFrequency(4e-1, 11)
        self.assertAlmostEqual(wp4 / wp1, 2.0, places = 10)

    def testMaxLinearGrowthFrequencyPositive(self):
        g = maximumLinearGrowthFrequency(1e-16, 1e-1, 1000.)
        self.assertGreater(g, 0.)

    def testMaxLinearGrowthFrequencyLinearInBeamDensity(self):
        g1 = maximumLinearGrowthFrequency(1e-16, 1e-1, 1000.)
        g2 = maximumLinearGrowthFrequency(2e-16, 1e-1, 1000.)
        self.assertAlmostEqual(g2 / g1, 2.0, places = 10)

    def testMaxLinearGrowthFrequencyInverseInLorentzFactor(self):
        g1 = maximumLinearGrowthFrequency(1e-16, 1e-1, 1000.)
        g2 = maximumLinearGrowthFrequency(1e-16, 1e-1, 2000.)
        self.assertAlmostEqual(g1 / g2, 2.0, places = 10)


if __name__ == '__main__':
    unittest.main()
