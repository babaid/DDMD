from openmm.app import Topology

def move_compound_by_vector(topology: Topology, positions, pos):
        for chain in topology.chains():
            for atom in chain.atoms():
                positions[atom.index] += pos
        return positions