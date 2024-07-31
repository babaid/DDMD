#!/bin/bash

python scripts/system_setup.py -p data/2f0z_fixed.pdb -l data/2f0z_ligand.mol2 -o data
python scripts/NVT.py -o data
python scripts/NPT.py -o data

