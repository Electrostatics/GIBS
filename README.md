# GIBS: Grand Canonical Monte Carlo (GCMC) simulation program for simulating ion-biomolecule interactions with molecular solvent models

The ionic environment of biomolecules strongly influences their structure, conformational stability, and 
inter-molecular interactions. GIBS is a Grand canonical Monte Carlo (GCMC) simulation program for computing the 
thermodynamic properties of ionic solutions and their distributions around biomolecules. GIBS implements 
algorithms that automate the excess chemical potential calculations for a given target salt concentration. GIBS 
uses a cavity-bias algorithm to achieve high sampling acceptance rates for inserting ions and solvent hard 
spheres when simulating dense ionic systems. In the current version, 
ion-ion 
interactions are 
de-
scribed using Coulomb, hard-sphere, or Lennard-Jones (L-J) potentials; solvent-ion
interactions are described using hard-sphere, L-J and attractive square-well potentials;
and, solvent-solvent interactions are described using hard-sphere repulsions.

interactions are described 
using Coulomb, hard-sphere, or Lennard-Jones (L-J) potentials; solvent-ion interactions are described using 
hard-sphere, L-J and attractive square-well potentials; and solvent-solvent interactions are described using 
hard-sphere repulsions. GIBS can be used as a platform for evaluating new implicit solvent and coarse-grained 
models for predicting the thermodynamics properties of ionic solutions. GIBS is written in C++ and is available 
freely for the community to use as an educational and research tool.

The GIBS program was written by Dr. Dennis G. Thomas in collaboration with Dr. Nathan A. Baker, at Pacific 
Northwest National Laboratory. The program was developed as part of projects funded by the National Institutes 
of Health through R01 Grant Nos. GM076121-04S1 and GM099450.

## Features

1. Automate the excess chemical potential calculations for bulk electrolyte solutions.
2. Fast and efficient GCMC sampling of ion distributions in bulk electrolyte solutions and around fixed molecular solutes 
3. Ion-Ion interactions can be specified `Coulomb`, `hard-sphere`, `Lennard-Jones` potentials.
4. Ion-Solvent interactions can be specified using `hard-sphere`, `Lennard-Jones`, 'attractive square well` potentials.
5. Solvent-Solvent interactions can be specified using `hard-sphere` repulsions.
6. Solvent representation: dielectric continuum (primitive model) or hard spheres (solvent primitive model 
(SPM)) 
7. Ion representation: hard spheres


## Download Instructions and Code Compilation

### File Folders:

1. `src` : This folder contains the source code files.

2. `inputfiles` : This folder contains templates of input files for the simulations.

3. `tutorials` : This folder contains example applications of GIBS.


### Code compilation and use


## Citing GIBS

To acknowledge the use of GIBS, please cite:

```


```

