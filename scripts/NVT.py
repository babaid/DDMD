import sys
sys.path.append("source")
import os

import argparse
from openmm import Platform, unit, CustomBondForce
from openmm.app import StateDataReporter, CheckpointReporter
from openmmtools import integrators
from openmm.app import Simulation, PDBReporter, DCDReporter, Modeller
from openmm import MonteCarloBarostat, LangevinIntegrator, NoseHooverIntegrator
from alive_progress import alive_bar
import numpy as np
from simulation_utils import save_system_state, load_system_state

parser = argparse.ArgumentParser(description="NVT Equilibriation")


parser.add_argument("-o", "--output", default='data', help="Output directory for simulation files.")
parser.add_argument("-s", "--steps", type=int, default=1e5, help="Number of steps")
parser.add_argument("-z", "--step-size", type=float, default=2, help="Step size (fs)")
parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
parser.add_argument("-p", "--platform", default="CUDA", help="Platform to perform computations on. [CUDA, CPU, Reference]")


def main(args):
    #load system
    system, pdbfile, _ = load_system_state(args.output, 0)

    combined_modeller = Modeller(pdbfile.topology, pdbfile.positions)

    
    defvals = {} 
    for force in system.getForces():
        if isinstance(force, CustomBondForce):
            defvals[force.getGlobalParameterName(0)] = 0

    platform = Platform.getPlatformByName(args.platform)
    
    integrator = LangevinIntegrator(300*unit.kelvin, 1/unit.picoseconds, args.step_size*unit.femtoseconds)

    simulation = Simulation(combined_modeller.topology, system, integrator, platform)
    simulation.context.setPositions(combined_modeller.positions)

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

    for j, (name, val) in enumerate(defvals.items()):
        print(name)
        simulation.context.setParameter(name, 0)
        simulation.context.setParameter(f'Toggle_{name}', 0)


    print("Minimizing energy.")
    simulation.minimizeEnergy()
    print("Done.")


        
    simulation.context.setVelocitiesToTemperature(0*unit.kelvin)
    print('Warming up the system...')
    temperature = 0*unit.kelvin
    temperature_end = 300.0*unit.kelvin
    integrator.setTemperature(temperature)

    K = 0.0

    mdsteps = args.steps
    with alive_bar(100, force_tty=True) as bar:
        for i in range(100):
            temperature += temperature_end/100
            K+=1.0
            for j, (name, val) in enumerate(defvals.items()):
                simulation.context.setParameter(name, K)
            integrator.setTemperature(temperature)
            simulation.step(int(mdsteps/100))
            bar()
    print("Done.")
    
    save_system_state(system, simulation, combined_modeller, 1, args.output)

if __name__ == '__main__':
    args = parser.parse_args()
    print("Performing NVT equilibriation.")
    main(args)