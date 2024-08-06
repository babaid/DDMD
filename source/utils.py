from biopandas.pdb import PandasPdb
import shutil

from rdkit import Chem
import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from rdkit import Chem


def expand_coordinates(coordinates, r, num_points=10):
            expanded_coords = []
            for coord in coordinates:
                for _ in range(num_points):
                    random_direction = np.random.randn(3)
                    random_direction /= np.linalg.norm(random_direction)
                    new_point = coord + r * random_direction
                    expanded_coords.append(new_point)
            return np.array(expanded_coords)
    
def find_pocket(pockets_dir, ligand_mol2_file):
    mol = Chem.MolFromMol2File(ligand_mol2_file)
    
    #COM of a molecule
    commol = np.array([0., 0., 0.])
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        commol+=np.array([positions.x, positions.y, positions.z])
    commol /= 10*mol.GetNumAtoms()

    
    pockets_dir = "tmp_complex_out/pockets/"
    path_pockets = [str(i) for i in Path(pockets_dir).iterdir()]
    path_pockets_pdb = [str(i) for i in Path(pockets_dir).iterdir() if PurePath(i).suffix == '.pdb']
    path_pockets_pqr = [str(i) for i in Path(pockets_dir).iterdir() if PurePath(i).suffix == '.pqr']

    pocket_points = {}
    for pock in path_pockets_pqr:
        g = re.findall('(?:pocket)(\d+)(?:_\w+)\.(\w+)', pock)
        i = g[0][0]
        suff = g[0][1]
        
        fname = pock.split(".")[1] + "_pqr.pdb"
        shutil.copyfile(pock, fname)
        
        p = PandasPdb().read_pdb(fname).df["ATOM"]
        coords = np.array([p["x_coord"], p["y_coord"], p["z_coord"]])
        coords = coords.reshape(com.shape[1], 3)
        pocket_points[i] = coords
    for k, v in pocket_points.items():
        coordinates = np.array(coms['7'])
        expanded_coords = expand_coordinates(coordinates, 2.5, 100)
        hull_delaunay = Delaunay(expanded_coords)
 
        if hull_delaunay.find_simplex(commol) >= 0:
            return k
        




def Mol2MolSupplier(file_path, sanitize=False):
    mols = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        start_indices = np.array([i for i, line in enumerate(lines) if '@<TRIPOS>MOLECULE' in line])
        end_indices = start_indices + np.array(list(start_indices[1:]) + [len(lines)])
        print(f"There are {len(start_indices)} molecule fragments in the mol2 file.")
        for start, end in zip(start_indices, end_indices):
            print(start, end)
            mol_block = ''.join(lines[start:end])
            mol = Chem.MolFromMol2Block(mol_block, sanitize=False)
            mols.append(mol)
            print(f"Molecule has {mol.GetNumAtoms()} atoms.")
    
    return mols














    
        