
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np

pdb = PDBFile("3.0m_box.pdb")
PDBFile.writeFile(pdb.topology, pdb.positions, open("3.0m_box.pdb", 'w'))

ff = ForceField("ff_our_nonbonded.xml") # this will create a NonbondedForce

modeller = Modeller(pdb.topology, pdb.positions)
# Add  water and cations
modeller.addSolvent(ff,model='spce',boxSize=(60.000000,60.000000,60.000000)*angstroms,positiveIon='Na+')

# Create the OpenMM system
print('Creating OpenMM System')

system = ff.createSystem(
    modeller.topology, nonbondedMethod=PME,
    nonbondedCutoff=1.280000*nanometers, constraints=AllBonds, rigidWater=True, ewaldErrorTolerance=0.0005)

# Create the integrator to do Langevin dynamics
integrator = LangevinIntegrator(
                        298.15*kelvin,       # Temperature of heat bath
                        1.0/picoseconds,  # Friction coefficient
                        2.0*femtoseconds, # Time step
)
integrator.setConstraintTolerance(0.00001)

# NPT ensemble
barostat = MonteCarloBarostat(1.0*bar, 298.15*kelvin, 25) 
system.addForce(barostat)
# Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = Platform.getPlatformByName('CUDA')

prop = dict(CudaPrecision='mixed',CudaDeviceIndex='0,1')

# Create the Simulation object
sim = Simulation(modeller.topology, system, integrator, platform)

print(platform.getPropertyValue(sim.context,'CudaDeviceIndex'))

# Set the particle positions
sim.context.setPositions(modeller.positions)
# Minimize the energy
print('Minimizing energy')
sim.minimizeEnergy(tolerance=1*kilojoule/mole, maxIterations=1000)
LocalEnergyMinimizer.minimize(sim.context,tolerance=1*kilojoule/mole,maxIterations=1000)

sim.context.setVelocitiesToTemperature(298.15*kelvin)

sim.reporters.append(DCDReporter('out.dcd', 1000))

sim.reporters.append(StateDataReporter(open('out_file', 'w'), 1000, step=True,
      potentialEnergy=True, totalEnergy=True, temperature=True, density=True,
          progress=True, remainingTime=True, speed=True, separator='	', totalSteps = 20000000))

print('Running Production...')
sim.step(20000000)
with open('out.chk', 'wb') as f:
      f.write(sim.context.createCheckpoint())
print('Saving pdb file')
positions = sim.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(sim.topology, positions, open('out.pdb', 'w'))
print('Done!')