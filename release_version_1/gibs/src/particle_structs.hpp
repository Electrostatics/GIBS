////////////////////////////////////////////////////////////////////////////////
// Copyright © 2017, Battelle Memorial Institute
// All rights reserved.
// 1.	Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any person or 
// 	entity lawfully obtaining a copy of this software and associated documentation files 
// 	(hereinafter “the Software”) to redistribute and use the Software in source and binary 
// 	forms, with or without modification.  Such person or entity may use, copy, modify, 
// 	merge, publish, distribute, sublicense, and/or sell copies of the Software, and may 
// 	permit others to do so, subject to the following conditions:
// •	Redistributions of source code must retain the above copyright notice, this list of 
// 	conditions and the following disclaimers. 
// •	Redistributions in binary form must reproduce the above copyright notice, this list of 
// 	conditions and the following disclaimer in the documentation and/or other materials 
// 	provided with the distribution. 
// •	Other than as used herein, neither the name Battelle Memorial Institute or Battelle may 
// 	be used in any form whatsoever without the express written consent of Battelle.  
// 2.	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
// 	CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
// 	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
// 	MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// 	DISCLAIMED. IN NO EVENT SHALL BATTELLE OR CONTRIBUTORS BE 
// 	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
// 	OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// 	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
// 	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
// 	AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
// 	LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
// 	IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
// 	THE POSSIBILITY OF SUCH DAMAGE.
////////////////////////////////////////////////////////////////////////////////
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
