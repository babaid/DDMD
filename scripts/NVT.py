import sys
sys.path.append("source")
import argparse

from openff.toolkit import Molecule
from openmm import XmlSerializer
from openmm import Platform, unit, CustomBondForce
from openmm.app import StateDataReporter, CheckpointReporter, PDBFile
from openmm.app import Simulation, PDBReporter, Modeller
from openmm import MonteCarloBarostat, LangevinIntegrator
from alive_progress import alive_bar
import numpy as np
from rdkit import Chem

from openmm import app
from openmm import unit, Vec3

from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from openmm.app import ForceField
from openmm.app import PDBFile

from topology_tools import move_compound_by_vector

pdb = PDBFile('data/1pin_docking.pdb')
rdkit_mol = Chem.MolFromMol2File("data/docked_meb.mol2")

npos = []
import numpy as np
for pos in pdb.positions:
    l = pos -  pdb.positions.mean()
    npos.append(l.value_in_unit(unit.nanometer))
npos = np.array(npos)

molecule = Molecule.from_rdkit(rdkit_mol, hydrogens_are_explicit=False)
molecule_topology = molecule.to_topology()


modeller = Modeller(pdb.topology, pdb.positions)
modeller.add(molecule_topology.to_openmm(),
             molecule_topology.get_positions().to_openmm())

positions = move_compound_by_vector(modeller.topology, modeller.positions, Vec3(3.0, 3.0, 3.0)*unit.nanometers)
modeller.positions = positions

gaff = GAFFTemplateGenerator(molecules=molecule)
forcefield = ForceField('amber/protein.ff14SB.xml',
                        'amber/tip3p_standard.xml',
                        'amber/tip3p_HFE_multivalent.xml')
forcefield.registerTemplateGenerator(gaff.generator)

modeller.addHydrogens(forcefield)
modeller.addSolvent(forcefield, model="tip3p", boxSize=Vec3(6.0, 6.0, 6.0)*unit.nanometers)

system = forcefield.createSystem(modeller.topology,
                                 nonbondedCutoff=1.1*unit.nanometers,
                                 switchDistance=0.9*unit.nanometers,
                                 constraints=app.HBonds,
                                 hydrogenMass=4.0*unit.amu,
                                 rigidWater=True, nonbondedMethod=app.PME)


platform = Platform.getPlatformByName('CPU')

integrator = LangevinIntegrator(300*unit.kelvin,
                                1/unit.picoseconds,
                                2*unit.femtoseconds)

simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)


simulation.reporters.append(PDBReporter('NVT.pdb', 1))
simulation.reporters.append(StateDataReporter(
        "NVT_report.csv",
        1000,
        step=True,
        potentialEnergy=True,
        totalEnergy=True,
        density=True,
        temperature=True,
        progress=True,
        totalSteps=1e5,
        separator='\t'))
simulation.reporters.append(CheckpointReporter("NVT.chk", 1e4))

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

state = simulation.context.getState(getPositions=True, getVelocities=True)
with open('NVT_state.xml', 'w') as output:
    output.write(XmlSerializer.serialize(state))

with open('topology.pdb', 'w') as output:
    PDBFile.writeFile(modeller.topology, state.po, output)
