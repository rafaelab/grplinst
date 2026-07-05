
import unittest
import numpy as np

from crpropa import *
from grplinst import *


###############################################################################
#                               test all classes                              #
###############################################################################

class TestImport(unittest.TestCase):

    def testAllPublicClassesAvailable(self):
        names = [
            'FlowHomogeneous', 
            'FlowJet1D',
            'MediumTemperatureHomogeneous', 
            'MediumDensityHomogeneous',
            'PlasmaInstabilityBroderick2012', 
            'PlasmaInstabilitySchlickeiser2012',
            'PlasmaInstabilitySironi2014', 
            'PlasmaInstabilityVafin2018',
            'PlasmaInstabilityBret2010TwoStream', 
            'PlasmaInstabilityBret2010Filamentation',
            'PlasmaInstabilityShalaby2020', 
            'PlasmaInstabilityMiniati2013',
            'plasmaFrequency', 
            'maximumLinearGrowthFrequency',
            'createFlowMiniati2013',
        ]
        import grplinst
        for name in names:
            self.assertTrue(hasattr(grplinst, name), f"grplinst.{name} not found")



###############################################################################
#                      vector ↔ Python list conversions                       #
###############################################################################

class TestVectorConversions(unittest.TestCase):
    """
    Verify that std::vector<double> SWIG typemaps accept Python lists and
    return objects that behave like sequences.
    """

    def testSetProfileFromList(self):
        dist = [0.0, 1.0, 2.0]
        dens = [3.0, 2.0, 1.0]
        flow = FlowJet1D()
        flow.setDistanceProfile(dist)
        flow.setDensityProfile(dens)
        self.assertEqual(list(flow.getDistanceProfile()), dist)
        self.assertEqual(list(flow.getDensityProfile()), dens)

    def testSetProfileFromNumpyArray(self):
        dist = np.array([0.0, 1.0, 2.0, 3.0])
        dens = np.array([4.0, 3.0, 2.0, 1.0])
        flow = FlowJet1D(dist, dens, 1.0, Vector3d(0, 0, 0), False)
        self.assertEqual(len(flow.getDistanceProfile()), 4)
        self.assertEqual(len(flow.getDensityProfile()), 4)

    def testConstructorFromNumpyArrayMatchesList(self):
        distList = [0.0, 1.0, 2.0]
        densList = [3.0, 2.0, 1.0]
        distNP = np.array(distList)
        densNP = np.array(densList)
        flowList = FlowJet1D(distList, densList, 1.0, Vector3d(0, 0, 0), False)
        flowNP = FlowJet1D(distNP, densNP, 1.0, Vector3d(0, 0, 0), False)
        self.assertEqual(
            list(flowList.getDistanceProfile()),
            list(flowNP.getDistanceProfile()),
        )



###############################################################################
#                   test directors: sub-classing of virtuals                  #
###############################################################################

class TestDirectorMediumDensity(unittest.TestCase):
    """
    SWIG directors let Python override pure-virtual C++ methods.
    A subclass of MediumDensity that returns a fixed value should be usable
    wherever a MediumDensity ref_ptr is accepted.
    """

    def testPythonSubclassIsAccepted(self):
        class ConstantDensity(MediumDensity):
            def getDensity(self, position, redshift = 0.):
                return 42.0

        dens = ConstantDensity()
        temp = MediumTemperatureHomogeneous(1e4)
        flow = FlowHomogeneous(1e38)

        # if the director wiring is broken this constructor call will raise
        pi = PlasmaInstabilityBroderick2012(flow, dens, temp)
        self.assertIsNotNone(pi)

    def testPythonSubclassValueIsForwardedToCpp(self):
        class ConstantDensity(MediumDensity):
            def getDensity(self, position, redshift = 0.):
                return 99.0

        dens = ConstantDensity()
        flow = FlowHomogeneous(1e38)
        temp = MediumTemperatureHomogeneous(1e4)

        # getMediumDensity returns the same object back through ref_ptr
        pi = PlasmaInstabilityBroderick2012(flow, dens, temp)
        retrieved = pi.getMediumDensity()
        self.assertAlmostEqual(retrieved.getDensity(Vector3d(0, 0, 0)), 99.0)


class TestDirectorMediumTemperature(unittest.TestCase):

    def testPythonSubclassIsAccepted(self):
        class ConstantTemperature(MediumTemperature):
            def getTemperature(self, position, redshift = 0.):
                return 1e6

        temp = ConstantTemperature()
        flow = FlowHomogeneous(1e38)
        dens = MediumDensityHomogeneous(0.1)
        pi = PlasmaInstabilityVafin2018(flow, dens, temp)
        self.assertIsNotNone(pi)

    def testPythonSubclassValueIsForwardedToCpp(self):
        class ConstantTemperature(MediumTemperature):
            def getTemperature(self, position, redshift = 0.):
                return 5e5

        temp = ConstantTemperature()
        flow = FlowHomogeneous(1e38)
        dens = MediumDensityHomogeneous(0.1)
        pi = PlasmaInstabilityVafin2018(flow, dens, temp)
        self.assertAlmostEqual(pi.getMediumTemperature().getTemperature(Vector3d(0, 0, 0)), 5e5)


class TestDirectorFlow(unittest.TestCase):

    def testPythonSubclassIsAccepted(self):
        class ConstantFlow(Flow):
            def getDensity(self, energy, position, redshift = 0.):
                return 1e-16

        flow = ConstantFlow()
        dens = MediumDensityHomogeneous(0.1)
        temp = MediumTemperatureHomogeneous(1e4)
        pi = PlasmaInstabilityBroderick2012(flow, dens, temp)
        self.assertIsNotNone(pi)

    def testPythonSubclassValueIsForwardedToCpp(self):
        class ConstantFlow(Flow):
            def getDensity(self, energy, position, redshift = 0.):
                return 7e-20

        flow = ConstantFlow()
        dens = MediumDensityHomogeneous(0.1)
        temp = MediumTemperatureHomogeneous(1e4)
        pi = PlasmaInstabilityBroderick2012(flow, dens, temp)
        retrieved = pi.getFlowProperties()
        self.assertAlmostEqual(retrieved.getDensity(1 * TeV, Vector3d(0, 0, 0)), 7e-20)


if __name__ == '__main__':
    unittest.main()
