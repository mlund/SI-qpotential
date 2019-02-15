
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np

elec_to_kJmol = (
    constants.elementary_charge**2 *
    AVOGADRO_CONSTANT_NA / (4*np.pi*1.0*8.854187817e-12
                            * constants.farad/constants.meter)).value_in_unit(kilojoule_per_mole*nanometer)

def findForce(system, forcetype, add=True):
# Finds a specific force in the system force list - added if not found
    for force in system.getForces():
        if isinstance(force, forcetype):
            return force
    if add==True:
        system.addForce(forcetype())
        return findForce(system, forcetype)
    return None

# return a q-potential system (TO 

pdb = PDBFile("1.0m_box.pdb")
PDBFile.writeFile(pdb.topology, pdb.positions, open("1.0m_box.pdb", 'w'))

ff = ForceField("ff_our_nonbonded.xml") # this will create a NonbondedForce

modeller = Modeller(pdb.topology, pdb.positions)
# Add  water and cations
modeller.addSolvent(ff,model='spce',boxSize=(60.000000,60.000000,60.000000)*angstroms,positiveIon='Na+')

ff = ForceField("ff_our_nonbonded_custom.xml") # this will create a CustomNonbondedForce

# Create the OpenMM system
print('Creating OpenMM System')

system = ff.createSystem(
    modeller.topology, nonbondedMethod=CutoffPeriodic,
    nonbondedCutoff=1.280000*nanometers, constraints=AllBonds, rigidWater=True)

def qPochhammerSymbol( Rc, moments ):
    if isinstance(Rc, Quantity):
        Rc = Rc / nanometer # strip unit
    qP = 1.0
    r = np.linspace(0, Rc, 5000)
    for i in range( moments ):
        qP *= (1 - (r/Rc)**(i+1) )
    return qP

qP = Continuous1DFunction( qPochhammerSymbol(1.280000*nanometers, 3), 0*nanometers, 1.280000*nanometers)

nonbonded = findForce(system, CustomNonbondedForce)
nonbonded.addTabulatedFunction( 'qP', qP )         # 'qP(r)' can now be used in energy function
nonbonded.addGlobalParameter( 'f', elec_to_kJmol ) # convert to kJ/mol

nonbonded.setEnergyFunction(
    'f * charge1 * charge2 * qP(r)/r'     ' + 4 * epsilon * ( (sigma/r)^12 - (sigma/r)^6 )'     ' ; sigma = 0.5 * ( sigma1+sigma2 ); epsilon = sqrt( epsilon1*epsilon2 )'
)

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