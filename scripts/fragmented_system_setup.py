
import sys
sys.path.append("source")
import os 
import argparse
from openff.toolkit import Molecule
from openmm import XmlSerializer
from openmm.app import  PDBFile,  Modeller, ForceField
from openmm import CustomCentroidBondForce, unit, Vec3
from alive_progress import alive_bar
import numpy as np
from rdkit import Chem
from openmm import app
from openmmforcefields.generators import GAFFTemplateGenerator


from simulation_utils import add_backbone_posres

from topology_tools import move_compound_by_vector

from utils import Mol2MolSupplier

import parmed as pmd

parser = argparse.ArgumentParser(description="System Setup")

parser.add_argument("-p", "--protein", required=True, help="Protein PDB file")
parser.add_argument("-l", "--ligand", required=True, help="Ligand molfile")
parser.add_argument("-o", "--output", default='data', help="Output directory for simulation files.")
parser.add_argument("--padding", type=float, default=1, help="Padding for solvent box (A)")
parser.add_argument("--water-model", default="tip3p",
                    choices=["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"],
                    help="Water model for solvation")
parser.add_argument("--positive-ion", default="Na+", help="Positive ion for solvation")
parser.add_argument("--negative-ion", default="Cl-", help="Negative ion for solvation")
parser.add_argument("--ionic-strength", type=float, default="0", help="Ionic strength for solvation")
parser.add_argument("--no-neutralize", action='store_true', help="Don't add ions to neutralize")
parser.add_argument("--protein-force-field", default='amber/ff14SB.xml', help="Protein force field")
parser.add_argument("--ligand-force-field", default='gaff-2.11', help="Ligand force field")
parser.add_argument("--water-force-field", default='amber/tip3p_standard.xml', help="Ligand force field")
parser.add_argument("--export-to-gmx", default=True, type=bool, help="Should the system be exported to gromacs.")


def main(args):

    pdb = PDBFile(args.protein)
    mols = Mol2MolSupplier(args.ligand, sanitize=False)

    for rdmol in mols:
        #Chem.AddHs(rdmol)
        Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    
    npos = []
    for pos in pdb.positions:
        l = pos -  pdb.positions.mean()
        npos.append(l.value_in_unit(unit.nanometer))
    npos = np.array(npos)

    molecules = [ Molecule.from_rdkit(rdmol,  allow_undefined_stereo=True) for rdmol in mols]
    molecules_topologies  = []

    for molecule in molecules:
        #molecule.generate_conformers()
        molecule.generate_unique_atom_names()
        molecules_topologies.append(molecule.to_topology())

    protein_chain_ids = []
    for chain in pdb.topology.chains():
        protein_chain_ids.append(chain.id)


    modeller = Modeller(pdb.topology, pdb.positions)

    abc = ["K", "L", "M", "N", "O"]
    for i, topology in enumerate(molecules_topologies):
        
        ot = topology.to_openmm()
        for chain in ot.chains():
            chain.id = abc[i]

        modeller.add(ot,
                topology.get_positions().to_openmm())
    


    gaff = GAFFTemplateGenerator(molecules=molecules)
    forcefield = ForceField(args.protein_force_field,
                            args.water_force_field,
                            'amber/tip3p_HFE_multivalent.xml')
    forcefield.registerTemplateGenerator(gaff.generator)

    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, model=args.water_model, padding = args.padding*unit.nanometers)

    system = forcefield.createSystem(modeller.topology,
                                    nonbondedCutoff=1.1*unit.nanometers,
                                    switchDistance=0.9*unit.nanometers,
                                    constraints=app.HBonds,
                                    hydrogenMass=1.5*unit.amu,
                                    rigidWater=True, nonbondedMethod=app.PME)


    protd = []
    for chain in modeller.topology.chains():
        if chain.id in protein_chain_ids:
            for atom in chain.atoms():
                protd.append(atom.index)

    # WE NEED TO RESTRAIN THE LIGAND TO STAY IN THE BINDING POCKET.
    restrain_lst = []
    for chain in modeller.topology.chains():
        if chain.id.startswith("X"):
            for atom in chain.atoms():
                if atom.name.startswith("C"):
                    restrain_lst.append(atom.index)

    """"
    molecule_position = molecule_topology.get_positions().to_openmm().value_in_unit(unit.nanometer)
    move_vec = npos.mean(axis=0) - molecule_position.mean(axis=0)

    force = CustomCentroidBondForce(2, '0.5*kres*step(distance(g1, g2)-r0)*(distance(g1, g2)-r0)^2')
    force.addGlobalParameter('r0', 1.3*np.sqrt(move_vec[0]**2+ move_vec[1]**2 + move_vec[2]**2)*unit.nanometers)
    force.addGlobalParameter('kres', 100*unit.kilojoule_per_mole/unit.nanometers**2)
    force.addGroup(restrain_lst, [1.0 for el in restrain_lst])
    force.addGroup(protd, [1.0 for el in protd])
    force.addBond([0, 1])
    system.addForce(force)
    """



    add_backbone_posres(system, modeller, 100)

    with open(os.path.join(args.output, 'systems',  'system_0.xml'), 'w') as output:
        output.write(XmlSerializer.serialize(system))
    
    with open(os.path.join(args.output, 'topologies', 'topology_0.pdb'), 'w') as output:
        PDBFile.writeFile(modeller.topology, modeller.positions, output)
    if args.export_to_gmx:
        structure = pmd.openmm.load_topology(modeller.topology, system, modeller.positions)
        structure.save(os.path.join(args.output, "systems", 'system.top'), format='gromacs', overwrite=True)
        structure.save(os.path.join(args.output, "systems", 'system.gro'), format='gromacs', overwrite=True)
if __name__ == '__main__':
    args = parser.parse_args()
    
    subdirs = ["states", "systems", "topologies", "reports", "checkpoints",  "prediction"]

    if not os.path.isdir(args.output):
        print("Creating output directory.")
        os.mkdir(args.output)
        for subdir in subdirs:
            os.mkdir(os.path.join(args.output, subdir))
    else:
        print("Output directory seems to already exist. This ma cause issues, take care.")
        for subdir in subdirs:
            if not os.path.isdir(os.path.join(args.output, subdir)):
                os.mkdir(os.path.join(args.output, subdir))

    print("Starting System Setup.")
    main(args)
