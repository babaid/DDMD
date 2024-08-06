import sys
sys.path.append("source")
import os
import argparse
from openmmtools import integrators
from openmm.app import Simulation, PDBReporter, CheckpointReporter, PDBFile, StateDataReporter,  DCDReporter, Modeller
from openmm import Platform, unit, XmlSerializer, MonteCarloBarostat, LangevinIntegrator, NoseHooverIntegrator
from alive_progress import alive_bar
import numpy as np


from simulation_utils import save_system_state, load_system_state



parser = argparse.ArgumentParser(description="NPT Equilibriation")


parser.add_argument("-o", "--output", default='data', help="Output directory for simulation files.")
parser.add_argument("-s", "--steps", type=int, default=2e5, help="Number of steps")
parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps")
parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
parser.add_argument("-p", "--platform", default="CUDA", help="Platform to perform computations on. [CUDA, CPU, Reference]")

def main(args):
    
    system, pdbfile, state = load_system_state(args.output, 1)
    combined_modeller = Modeller(pdbfile.topology, pdbfile.positions)

    print("NPT Equilibriation")
    platform =Platform.getPlatformByName("CUDA")
    integrator = NoseHooverIntegrator(300*unit.kelvin, 1/unit.picoseconds, 1.0*unit.femtoseconds)
    simulation = Simulation(combined_modeller.topology, system, integrator, platform)
    simulation.context.setState(state)

    simulation.reporters.append(PDBReporter(os.path.join(args.output, 'reports', 'NPT.pdb'),  args.interval))# atomSubset=subset))   
    simulation.reporters.append(StateDataReporter(
                                                    os.path.join(args.output, 'reports', 'NPT_report.csv'),
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

        

    #barostat
    system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin, 25))
    simulation.context.reinitialize(True)


    mdsteps = 2e5
    with alive_bar(100, force_tty=True) as bar:
        for i in range(100):
            simulation.step(mdsteps/100)
            #simulation.context.setParameter('kres', 100-i*0.98)
            simulation.context.setParameter('kbb', 100-i*0.98)
            bar()

    #simulation.context.setParameter('kres',0)
    simulation.context.setParameter('kbb', 0)

    save_system_state(system, simulation, combined_modeller, 2, args.output)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)