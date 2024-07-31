import argparse
import sys
sys.path.append("source")

from openmmtools import integrators
from openmm import Platform, unit, CustomBondForce
from openmm.app import Simulation, PDBReporter, CheckpointReporter, PDBFile, StateDataReporter,  DCDReporter, Modeller
from openmm import Platform, unit, XmlSerializer, MonteCarloBarostat, LangevinIntegrator, NoseHooverIntegrator
from alive_progress import alive_bar
import numpy as np
import argparse
import os


parser = argparse.ArgumentParser(description="NVT Equilibriation")


parser.add_argument("-o", "--output", default='data', help="Output directory for simulation files.")
parser.add_argument("-s", "--steps", type=int, default=1e5, help="Number of steps")
parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps")
parser.add_argument("-l", "--ligand", tpye=str)
parser.add_argument("-p", "--platform", default="CUDA", help="Platform to perform computations on. [CUDA, CPU, Reference]")


def main(args):

    with open(os.path.join(args.output, "systems", "system_3.xml")) as input:
        system = XmlSerializer.deserialize(input.read())
        
    pdbfile = PDBFile(os.path.join(args.output, "topologies", "topology_2.pdb"))

    modeller = Modeller(pdbfile.topology, pdbfile.positions)

    with open(os.path.join(args.output, "states", "state_2.xml")) as input:
        state = XmlSerializer.deserialize(input.read())




    print("Production Run")
    platform =Platform.getPlatformByName("CUDA")
    integrator = NoseHooverIntegrator(300*unit.kelvin, 1/unit.picoseconds, 1*unit.femtoseconds)
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setState(state)

    simulation.reporters.append(PDBReporter('production.pdb', 1000)) 
    simulation.reporters.append(StateDataReporter(
                                                    "production_report.csv",
                                                    1000,
                                                    step=True,
                                                    potentialEnergy=True,
                                                    totalEnergy=True,
                                                    temperature=True,
                                                    density=True,
                                                    progress=True,
                                                    totalSteps=2e6,
                                                    separator='\t'
                                                    ))

    simulation.reporters.append(CheckpointReporter("production.chk", 1e4))
        
    mdsteps = 2e6



    simulation.step(mdsteps)



    state = simulation.context.getState(getPositions=True, getVelocities=True)


    with open('PROD_state.xml', 'w') as output:
        output.write(XmlSerializer.serialize(state))

    with open('topology_prod.pdb', 'w') as output:
        PDBFile.writeFile(pdbfile.topology, state.getPositions(), output)
if __name__ == '__main__':
    
    args = parser.parse_args()
    main(args)
