import sys
sys.path.append("source")
import os
import argparse
from openmmtools import integrators
from openmm import Platform, unit, CustomBondForce
from openmm.app import Simulation, PDBReporter, CheckpointReporter, PDBFile, StateDataReporter,  DCDReporter, Modeller
from openmm import Platform, unit, XmlSerializer, MonteCarloBarostat, LangevinIntegrator, NoseHooverIntegrator
from alive_progress import alive_bar
import numpy as np
#

from simulation_utils import load_system_state, save_system_state


parser = argparse.ArgumentParser(description="NPT Equilibriation")


parser.add_argument("-o", "--output", default='data', help="Output directory for simulation files.")
parser.add_argument("-s", "--steps", type=int, default=2e5, help="Number of steps")
parser.add_argument("-z", "--step-size", type=float, default=0.1, help="Step size (ps")
parser.add_argument("-i", "--interval", type=int, default=10000, help="Reporting interval")
parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
parser.add_argument("-p", "--platform", default="CUDA", help="Platform to perform computations on. [CUDA, CPU, Reference]")



def main(args):

    system, pdbfile, state = load_system_state(args.output, 2)
    print(state)
    combined_modeller = Modeller(pdbfile.topology, pdbfile.positions)

    print("Production Run")
    platform =Platform.getPlatformByName(args.platform)
    integrator = NoseHooverIntegrator(300*unit.kelvin, 1/unit.picoseconds, args.step_size*unit.femtoseconds)
    simulation = Simulation(combined_modeller.topology, system, integrator, platform)
    simulation.context.setState(state)

    simulation.reporters.append(PDBReporter(os.path.join(args.output, 'reports', 'production.pdb'),  args.interval)) 
    simulation.reporters.append(StateDataReporter(
                                                    os.path.join(args.output, 'reports', 'production_report.csv'),
                                                    args.interval,
                                                    step=True,
                                                    potentialEnergy=True,
                                                    totalEnergy=True,
                                                    temperature=True,
                                                    density=True,
                                                    progress=True,
                                                    totalSteps=args.steps,
                                                    separator='\t'
                                                    ))
    simulation.reporters.append(CheckpointReporter(os.path.join(args.output, 'checkpoints', 'NPT.chk'), args.interval))
        
    mdsteps = 2e4
        
    with alive_bar(100, force_tty=True) as bar:
            simulation.step(mdsteps/100)
            bar()

    save_system_state(system, simulation, combined_modeller.topology, 3, args.output)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)