![GIBS](/images/gibs_icon.png)
# GIBS: Grand Canonical Monte Carlo (GCMC) simulation program for computing ion distributions around biomolecules with hard sphere solvent models

GIBS is a Grand canonical Monte Carlo (GCMC) simulation program for computing the 
thermodynamic properties of ionic solutions and their distributions around 
biomolecules. GIBS implements algorithms that automate the excess chemical potential 
calculations for a given target salt concentration. GIBS uses a cavity-bias 
algorithm to achieve high sampling acceptance rates for inserting ions and solvent 
hard spheres when simulating dense ionic systems. In the current version, 
ion-ion interactions can be described using Coulomb, hard-sphere, or Lennard-Jones (L-J)
potentials; solvent-ion interactions can be described using hard-sphere, L-J and 
attractive square-well potentials; and, solvent-solvent interactions are described 
using hard-sphere repulsions. GIBS can be used as a platform to evaluate new implicit 
solvent and coarse-grained models for predicting the thermodynamics properties of 
ionic solutions. GIBS is written in C++ and is available freely for the community to 
use as an educational and as a research tool.

The GIBS program was written by Dr. Dennis G. Thomas in collaboration with Dr. Nathan A. Baker, at the Pacific Northwest National Laboratory. The program was developed as part of projects funded by the National Institutes of Health through R01 Grant Nos. GM076121-04S1 and GM099450.

## Main Features

1. Automated excess chemical potential calculations for bulk electrolyte solutions.
2.	Fast and efficient GCMC sampling of ion distributions in bulk electrolyte solutions and around fixed molecular solutes.
3.	Models for Ion-Ion interactions using `Coulomb`, `hard-sphere`, `Lennard-Jones` potentials.
4.	Models for Ion-Solvent interactions using `hard-sphere`, `Lennard-Jones`, `attractive square well` potentials.
5.	Models for Solvent-Solvent interactions using `hard-sphere` repulsions.
6.	Solvent representation as dielectric continuum (primitive model) or as hard spheres (solvent primitive model or SPM).
7.	Ion representation as charged hard spheres.

## Current Version

The current version is [version 1](/release_version_1) 

## Citing GIBS

A publication has been planned.

```

```
## Disclaimer
Disclaimer
This material was prepared as an account of work sponsored by an agency of the 
United States Government.  Neither the United States Government nor the United 
States Department of Energy, nor Battelle, nor any of their employees, nor any 
jurisdiction or organization that has cooperated in the development of these 
materials, makes any warranty, express or implied, or assumes any legal liability 
or responsibility for the accuracy, completeness, or usefulness or any 
information, apparatus, product, software, or process disclosed, or represents 
that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service by trade 
name, trademark, manufacturer, or otherwise does not necessarily constitute or 
imply its endorsement, recommendation, or favoring by the United States 
Government or any agency thereof, or Battelle Memorial Institute. The views and 
opinions of authors expressed herein do not necessarily state or reflect those of 
the United States Government or any agency thereof.
PACIFIC NORTHWEST NATIONAL LABORATORY
operated by
BATTELLE
for the
UNITED STATES DEPARTMENT OF ENERGY
under Contract DE-AC05-76RL01830

