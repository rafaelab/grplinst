# grplinst
Module for the CRPropa code to calculate the effect of plasma instabilities on the development of electromagnetic cascades


## Installation Instruction

To install *grplinst*, you will need to have CRPropa 3 installed. 
Go to https://github.com/CRPropa/CRPropa3/ and follow the instructions.

Now proceed with the installation of this module.

1. Download the latest version of the code.
```
git clone https://github.com/rafaelab/grplinst.git
```

2. Navigate into the downloaded folder and create a folder called "build/".

3. Install the code with CMake and then Make:

```
cmake ..
make
```

4. If the code was properly compiled, you are ready to go!
You can now include it as a CRPropa module as follows:

```
from crpropa import *
import grplinst as plinst


z = 0.14 # redshift of the source
nIGM = 1e-13 # density of the intergalactic medium in m^-3
L = 1e38 # luminosity of the beam in W
T = 1e5 # temperature of the intergalactic medium

obs = Observer()
obs.add(ObserverPoint())
output = TextOutput('test.txt', Output.Event1D)
output.setEnergyScale(eV)
obs.onDetection(output)

sim = ModuleList()
sim.add(SimplePropagation())
sim.add(Redshift())
sim.add(EMInverseComptonScattering(CMB, True))
sim.add(EMInverseComptonScattering(EBL, True))
sim.add(EMPairProduction(CMB, True))
sim.add(EMPairProduction(EBL, True))
sim.add(plinst.PlasmaInstability(nIGM, T, L, instabilityModel))
sim.add(MaximumTrajectoryLength(4000 * Mpc))
sim.add(MinimumEnergy(1e9 * eV))
sim.add(obs)

source = Source()
source.add(SourcePowerLawSpectrum(1e9 * eV, 1e14 * eV, -1))
source.add(SourceParticleType(22))
source.add(SourcePosition(Vector3d(redshift2ComovingDistance(z), 0, 0)))
source.add(SourceRedshift1D())
source.add(SourceDirection(Vector3d(-1, 0, 0)))


sim.setShowProgress(True)
sim.run(source, 10000, True)
```







