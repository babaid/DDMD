# Introduction

Recently deep generative models are used more and more in computer aided drug discovery. Typical models used are Diffusion Models, VAE's and different methods utilizing GNN's. In this project I will try to utilize [DeepICL](https://github.com/ACE-KAIST/DeepICL/tree/master), to combine a generative model with a Molecular dynamics simulation. The overally goal is an improvement over a classical Docking-Simulation-Screening pipeline.

Similar how [Cosolvent Molecular Dynamics](https://pubs.acs.org/doi/10.1021/acs.jmedchem.6b00399) allow for finding states that usually arent available at docking (e.g. binding pockets forming later on in a simulation), this proposed pipeline should also allow for conveniences.

# Simulation Pipeline

The simulation will have following steps:

1. System Setup with constraints 
2. NVT Equilibriation
3. NPT Equilibriation
4. Production Run for N steps
5. Generating a new ligand based on the current topology
6. Repeat steps 4-5 k times.

