import sys
sys.path.append("source")

from openmmtools import integrators
from openmm import Platform, unit, CustomBondForce
from openmm.app import Simulation, PDBReporter, CheckpointReporter, PDBFile, StateDataReporter,  DCDReporter, Modeller
from openmm import Platform, unit, XmlSerializer, MonteCarloBarostat, LangevinIntegrator, NoseHooverIntegrator
from alive_progress import alive_bar
import numpy as np





def main():

    with open('system_NPT.xml') as input:
        system = XmlSerializer.deserialize(input.read())
        
    pdbfile = PDBFile('topology_npt.pdb')

    modeller = Modeller(pdbfile.topology, pdbfile.positions)

    with open('NPT_state.xml') as input:
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
