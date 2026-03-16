import sys
sys.path.append('../build')

from crpropa import *
from grplinst import *


# background photon fields
ebl = IRB_Saldana21()
cmb = CMB()

# source parameters
z = 0.0308  # Mrk 421
nIGM = 1e-1 # m^-3
L = 1e38 # W
T = 1e4 # K

# plasma instability configuration
temperature = MediumTemperatureHomogeneous(T)
density = MediumDensityHomogeneous(nIGM)
beam = FlowHomogeneous(L, Vector3d(0, 0, 0))
plinst = PlasmaInstabilityBroderick2012(beam, density, temperature)


photons = electrons = True
ppEBL = EMPairProduction(ebl, electrons)
ppCMB = EMPairProduction(cmb, electrons)
icEBL = EMInverseComptonScattering(ebl, photons)
icCMB = EMInverseComptonScattering(cmb, photons)
# plinst = PlasmaInstability(model, beam, density, temperature)
redshift = Redshift()
processes = [redshift, ppEBL, ppCMB, icEBL, icCMB, plinst]

maxTrajectory = MaximumTrajectoryLength(4000 * Mpc)
minEnergy = MinimumEnergy(1e9 * eV)
breakCondition = [maxTrajectory, minEnergy]

source = Source()
source.add(SourcePowerLawSpectrum(1e9 * eV, 1e13 * eV, -1))
source.add(SourceParticleType(22))
source.add(SourcePosition(Vector3d(redshift2ComovingDistance(z), 0, 0)))
source.add(SourceRedshift1D())
source.add(SourceDirection(Vector3d(-1, 0, 0)))

output = TextOutput('test.txt', Output.Event1D)
output.setEnergyScale(eV)
output.set(output.WeightColumn, True)

observer = Observer()
observer.add(Observer1D())
observer.onDetection(output)

sim = ModuleList()
sim.add(SimplePropagation(1e-3 * kpc, 10 * Mpc))
for p in processes:
	sim.add(p)
for bc in breakCondition:
	sim.add(bc)
sim.add(observer)

sim.setShowProgress(True)
sim.run(source, 10000, True)
