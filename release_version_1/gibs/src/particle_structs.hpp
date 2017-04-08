/*
 * particle_structs.hpp
 *
 *  Created on: Aug 7, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file particle_structs.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains declarations of data structures for particles in the simulation box
 */

#ifndef PARTICLE_STRUCTS_HPP_
#define PARTICLE_STRUCTS_HPP_

/**
 * @brief Data structure for each particle in the simulation box
 */
struct particle_struct{

	/*Floating point variable to store charge on the particle (in e units)*/
	double charge;
	/*Floating point variable to store the particle radius (in Angstrom)*/
	double radius;
	/* Floating point value to store the LJ potential well depth of
 	 	 the particle type (kcal/mol)*/
	double ljeps;
	/* Floating point variable to store the LJ collision diameter of the
	 * particle type(kcal/mol)*/
	double ljsig;
	/* Floating point variables to store the x,y,z coordinate of the particle's position*/
	double x,y,z;
	/*Integer variable to store A non-zero positive integer assigned to the particle to identify its type.
	 If solvent is represented as hard spheres, then its ptype = nparticletypes */
	int ptype;
	/* Integer variables to store the x,y,z coordinate cell indices of the cavity grid
	 * containing the center of the particle.*/
	int gridcell_ix,gridcell_iy,gridcell_iz;
	/*Integer variable to store the 1D cavity grid cell index of the particle */
	int gridcell_1Dindex;

	/*Integer variable to store the LJ 1D cell index containing the center of the particle*/
	int lj_gridcell_1Dindex;
	/* String variable to store a label that identifies the type of particle (e.g., Na, Cl, Water, etc.) */
	std::string label;

};
typedef struct particle_struct struct_particle;

#endif /* PARTICLE_STRUCTS_HPP_ */
