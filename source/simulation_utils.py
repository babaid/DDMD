import os
from openmm import CustomExternalForce, unit
from pymol import cmd

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
  
