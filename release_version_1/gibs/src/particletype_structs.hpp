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
 * particletype_structs.hpp
 *
 *  Created on: Aug 7, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file particletype_structs.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains declarations of data structures for particle types in the simulation box
 */

#ifndef PARTICLETYPE_STRUCTS_HPP_
#define PARTICLETYPE_STRUCTS_HPP_


struct ljcell_cavityindices_struct{

	int num;
	std::vector<int> cavity_1Dcellindices;
};

typedef struct ljcell_cavityindices_struct struct_ljcell_cavityindices;


struct particletype_struct{

	/*Boolean variable to indicate whether there is a Lennard-Jones potential interaction
	 * between the two particle types of a particle pair type*/
	bool lj_switch;

	/*A string label that identifies the type of particle (e.g., Na, Cl, Water, etc.)*/
	std::string label;

	/* A floating point variable to store and update the number of each particle type present
	 * in the simulation box at a given simulation step.*/
	double num;
	/* Integer variable to store the number of cavity grid cells available for
	 * placing the center of a particle type hard sphere*/
	int num_cavities;
	/* Floating point variable to store the charge of the particle type (in e units) */
	double charge;
	/* Floating point variable to store the radius of the particle type (in Angstrom) */
	double radius;

	/*Floating point variable to store the LJ potential well depth of the
	 * particle type (kcal/mol) */
	double ljeps;
	/* Floating pointe variable to store the LJ collision diameter of the
	 * particle type(Angstrom) */
	double ljsig;

	/*Integer variable to store the non-zero positive integer assigned to the
	 * particle to identify its type. If solvent is represented as hard
	 * spheres, then its ptype = the number of particle types*/
	int ptype;

	/*Floating point variable to store the bulk concentration of the
	 * particle type in molar (M) units*/
	double conc;

	/* Floating point variable to store the molar mass (g/mol) of
	 * the particle type */

	double molmass;

	/*Floating point variable to store the excess chemical potential of the
	 * particle type (kcal/mol) */
	double mu_ex;
	/*Floating point variable to store the total chemical potential of the
	 * particle type (kcal/mol) from previous iteration */
	double mu_old;

	/*Floating point variable to store the total chemical potential of the
	 * particle type (kcal/mol) in the current iteration*/
	double mu_new;

	/*Floating point variable to store the average number of the
	 * particle type after each iteration of GCMC chemical potential simulation */
	double num_avg;

	/*Floating point variable to store the average concentration (in Molar units) of the
		 * particle type after each iteration of GCMC chemical potential simulation */
	double conc_avg;

	/*Array of integers to store and update the number of cavity grid cells
	 * available for the particle type in each cavity grid segment*/
	int *num_cavities_in_segment;

	double epsi;
	double epsim1;
	double Bip1;
	double Bi ;

	double targ_dens;
	double targ_mu;
	double targ_num;
	double mu_ex_m1;
	double mu_i;


	std::tr1::unordered_map<int, struct_ljcell_cavityindices > ljcell_cavityindices;

};
typedef struct particletype_struct struct_particletype;


#endif /* PARTICLETYPE_STRUCTS_HPP_ */
