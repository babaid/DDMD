# Introduction

Recently deep generative models are used more and more in computer aided drug discovery. Typical models used are Diffusion Models, VAE's and different methods utilizing GNN's. In this project I will try to utilize different generative models. Here I try to focus on pocket-specific models, especially molecular linkers and fragment based models in combination with a molecular dynamics simulation. The overall goal is an improvement over a classical Docking-Simulation-Screening pipeline.

A list of the models I am currently exploring:

- [DeepICL](https://github.com/ACE-KAIST/DeepICL/tree/master): in theory it is the perfect model for this task as it works with interactions between a starting ligand and a pocket. In reality it doesnt seem to work with any protein, I am yet to figure out why. It also has some unclarities in design which I try to fix currently, without breaking anything.

- [autofragdiff](https://github.com/ghorbanimahdi73/autofragdiff): Very promising diffusion based appraoch that also works with the binding pocket directly.



Similar to how [Cosolvent Molecular Dynamics](https://pubs.acs.org/doi/10.1021/acs.jmedchem.6b00399) allow for finding states that usually arent available at docking (e.g. binding pockets forming later on in a simulation), this proposed pipeline should also allow for similar sideffects.

# Simulation Pipeline

The simulation will have following steps:

1. System Setup
2. NVT Equilibriation
3. NPT Equilibriation
4. Production Run for N steps
5. Generating a new ligand based on the current topology:
6. Repeat steps 4-5 for k times.

For now I will work with protein-ligand complex f20z used in the DeepICL example, as the binding pocket is quite well defined and yet the protein isnt too big. For simplification I deleted most of the ligand, leaving only the ring behind.

# Simulation Stages

For clarity, each stage, i.e.
    System Setup, 
    NVT, 
    NPT, 
    Production and 
    Ligand Generation 
is in a seperate script. Each step saves a state and a topology of the last simulation step and there are also reports, pdb-reports and checkpoints for safety and analysis. Each of these stages has a corresponding folder. There is also a subfolder for the ligand generation process. 

This allows continuation and dynamic fixing of issues that come up during simulation, i.e. one doesnt have to start over from 0 if ligand generation fails, or can even start multiple pipelines starting from the same basis state, e.g. starting after NPT equilibriation.

# Issues

A possible issue is the seeting of velocities of the generated new atoms. Currently I think the easiest is to start with 0 or the thermal velocity derived from the temperature.
