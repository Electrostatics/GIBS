# GIBS: Grand Canonical Monte Carlo (GCMC) simulation program for simulating ion-biomolecule interactions with molecular solvent models

The ionic environment of biomolecules strongly influences their structure, conformational stability, and 
inter-molecular interactions. GIBS is a Grand canonical Monte Carlo (GCMC) simulation program for computing the 
thermodynamic properties of ionic solutions and their distributions around biomolecules. GIBS implements 
algorithms that automate the excess chemical potential calculations for a given target salt concentration. GIBS 
uses a cavity-bias algorithm to achieve high sampling acceptance rates for inserting ions and solvent hard 
spheres when simulating dense ionic systems. In the current version of GIBS, ion-ion interactions are described 
using Coulomb, hard-sphere, or Lennard-Jones (L-J) potentials; solvent-ion interactions are described using 
hard-sphere, L-J and attractive square-well potentials; and solvent-solvent interactions are described using 
hard-sphere repulsions. GIBS can be used as a platform for evaluating new implicit solvent and coarse-grained 
models for predicting the thermodynamics properties of ionic solutions. GIBS is written in C++ and is available 
freely for the community to use as an educational and research tool.

The GIBS program was written by Dr. Dennis G. Thomas in collaboration with Dr. Nathan A. Baker, at Pacific 
Northwest National Laboratory. The program was developed as part of projects funded by the National Institutes 
of Health through R01 Grant Nos. GM076121-04S1 and GM099450.


## Download Instructions and Code Compilation


## Citing GIBS

Please cite the following paper for using GIBS:

```
Dennis G. Thomas and Nathan A. Baker, GIBS: Grand Canonical Monte Carlo (GCMC) simulation program for simulating 
ion-biomolecule interactions. XYZ 2017...

```

