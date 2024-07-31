import os
import sys
import argparse


parser = argparse.ArgumentParser(description="NVT Equilibriation")


parser.add_argument("-o", "--output", default='data', help="Output directory for simulation files.")
parser.add_argument("--N", type=int, default=0, help="Number of step")
parser.add_argument("-r", "--results", type=str, default='lig_res', help="Step size (ps")
parser.add_argument("-p", "--protein", type=str,  help="Protein")
parser.add_argument("-l", "--ligand", type=str,  help="Ligand")




def generate(args):


if __name__ =='__main__':
    args = parser.parse_args()
    generate(args)