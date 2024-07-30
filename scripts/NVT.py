import sys
sys.path.append("source")
import argparse
import os

#data imports
from openmm import XmlSerializer
from openmm.app import PDBFile
from openmm.app import StateDataReporter, CheckpointReporter, PDBFile, PDBReporter


from openmm import Platform, unit
from openmm.app import Simulation, Modeller
from openmm import LangevinIntegrator

from alive_progress import alive_bar



parser = argparse.ArgumentParser(description="NVT Equilibriation")


parser.add_argument("-o", "--output", default='data', help="Output directory for simulation files.")
parser.add_argument("-s", "--steps", type=int, default=1e5, help="Number of steps")
parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps")
parser.add_argument("-f", "--friction-coeff", type=float, default=1, help="Friction coefficient (ps)")
parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
parser.add_argument("-p", "--platform", default="CUDA", help="Platform to perform computations on. [CUDA, CPU, Reference]")

def simulate(args):
    
    with open(os.path.join(args.output, 'systems', 'system_0.xml')) as input:
        system = XmlSerializer.deserialize(input.read())
        
    pdbfile = PDBFile(os.path.join(args.output, 'topologies', 'topology_0.pdb'))

    modeller = Modeller(pdbfile.topology, pdbfile.positions)

    platform = Platform.getPlatformByName(args.platform)

    integrator = LangevinIntegrator(args.temperature*unit.kelvin,
                                    1/unit.picoseconds,
                                    args.step_size*unit.femtoseconds)

    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)


    simulation.reporters.append(PDBReporter(os.path.join(args.output, 'reports', 'NVT.pdb'), args.interval))
    simulation.reporters.append(StateDataReporter(
            os.path.join(args.output, 'reports', 'NVT_report.csv'),
            args.interval,
            step=True,
            potentialEnergy=True,
            totalEnergy=True,
            density=True,
            temperature=True,
            progress=True,
            totalSteps=args.steps,
            separator='\t'))
    simulation.reporters.append(CheckpointReporter(os.path.join(args.output, 'checkpoints', 'NVT.chk'), args.interval))

    print("Minimizing energy.")
    simulation.minimizeEnergy()
    print("Done.")

    simulation.context.setVelocitiesToTemperature(0*unit.kelvin)

    print('Warming up the system...')
    temperature = 0*unit.kelvin
    temperature_end = 300.0*unit.kelvin
    integrator.setTemperature(temperature)

    mdsteps = 1e5
    with alive_bar(100, force_tty=True) as bar:
        for i in range(100):
            temperature += temperature_end/100
            integrator.setTemperature(temperature)
            simulation.step(int(mdsteps/100))
            bar()
    print("Done.")

    with open(os.path.join(args.output, 'systems', 'system_1.xml'), 'w') as output:
        output.write(XmlSerializer.serialize(system))


    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(os.path.join(args.output, 'states', 'state_1.xml'), 'w') as output:
        output.write(XmlSerializer.serialize(state))

    with open(os.path.join(args.output, 'topologies', 'topology_1.pdb'), 'w') as output:
        PDBFile.writeFile(modeller.topology, state.getPositions(), output)

if __name__ == '__main__':
    args = parser.parse_args()
    simulate(args)