import os
from openmm import CustomExternalForce, unit
from pymol import cmd
from openmm import XmlSerializer
from openmm.app import PDBFile
def add_backbone_posres(system, modeller, restraint_force):
  force = CustomExternalForce("kbb*periodicdistance(x, y, z, x0, y0, z0)^2")
  force_amount = restraint_force * unit.kilocalories_per_mole/unit.angstroms**2
  force.addGlobalParameter("kbb", force_amount)
  force.addPerParticleParameter("x0")
  force.addPerParticleParameter("y0")
  force.addPerParticleParameter("z0")
  for i, (atom_crd, atom) in enumerate(zip(modeller.positions, modeller.topology.atoms())):
    if atom.name in  ('CA', 'C', 'N'):
      force.addParticle(i, atom_crd.value_in_unit(unit.nanometers))
  #posres_sys = deepcopy(system)
  system.addForce(force)


def seperate_complex(path, output):
  cmd.load(path)
  cmd.select('protein', 'polymer')
  cmd.select('ligand', 'organic')
  cmd.save(os.path.join(output, 'protein.pdb'), 'protein')
  cmd.save(os.path.join(output, 'ligand.sdf'), 'ligand')

def save_system_state(system, simulation, pdb_top, id: int, output_root_dir: str):
  state = None
  if simulation is not None:
      state = simulation.context.getState(getPositions=True, getVelocities=True)
      with open(os.path.join(output_root_dir, 'states', f'state_{id}.xml'), 'w') as output:
          output.write(XmlSerializer.serialize(state))
  if system:
      with open(os.path.join(output_root_dir, 'systems', f'system_{id}.xml'), 'w') as output:
          output.write(XmlSerializer.serialize(system))
  if pdb_top is not None:
      with open(os.path.join(output_root_dir, 'topologies', f'topology_{id}.pdb'), 'w') as output:
          if state is not None:
              PDBFile.writeFile(pdb_top.topology, state.getPositions(), output)
          else:
              PDBFile.writeFile(pdb_top.topology, pdb_top.positions, output)



def load_system_state(root_dir, id):
    with open(os.path.join(root_dir, "systems", f"system_{id}.xml")) as input:
        system = XmlSerializer.deserialize(input.read())

    pdbfile = PDBFile(os.path.join(root_dir, "topologies", f"topology_{id}.pdb"))
    state = None
    if id != 0:
        with open(os.path.join(root_dir, 'states', f"state_{id}.xml")) as input:
            state = XmlSerializer.deserialize(input.read())
    return system, pdbfile, state 