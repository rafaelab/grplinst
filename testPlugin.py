from crpropa import *
import grplinst as plinst


EBL = IRB_Gilmore12
z = 0.14 # 1ES 0209+200
nIGM = 1e-7 # cm^-3
L = 1e38 # W
T = 1e5 # K
model = 'C'

obs = Observer()
obs.add(ObserverPoint())
output = TextOutput('test.txt', Output.Event1D)
output.setEnergyScale(eV)
output.set(output.WeightColumn, True)
obs.onDetection(output)

sim = ModuleList()
sim.add(SimplePropagation(1e-10 * kpc, 10 * kpc))
sim.add(Redshift())
sim.add(EMPairProduction(CMB, True))
sim.add(EMPairProduction(EBL, True))
sim.add(EMInverseComptonScattering(CMB, True))
sim.add(EMInverseComptonScattering(EBL, True))
sim.add(plinst.PlasmaInstability(nIGM, T, L, model))
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
