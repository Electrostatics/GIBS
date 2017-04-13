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
 * solute_structs.hpp
 *
 *  Created on: Aug 12, 2015
 *      Author:  Dennis G. Thomas
 *
 * @file solute_structs.hpp
 * @author Dennis G. Thomas
 *
 * @brief Contains declarations of data structures for solute atoms and electrostatic
 * 			potential
 *
 */

#ifndef SOLUTE_STRUCTS_HPP_
#define SOLUTE_STRUCTS_HPP_

#include "generic.hpp"
#include "constants.hpp"

/**
 * @brief Data structure of an atom
 */
struct atom_struct {

// Integer variable to store the atom index
	int atom_number;
	// Integer variable to store the residue index (if residues are present)
	int residue_number;
	// Floating point variable to store the x, y, and z coordinates of the atom
	double x, y, z;
	// Floating point variable to store the radius of the atom in Angstrom units
	double radius;
	// Floating point variable to store the charge of the atom in e units
	double charge;
	/*Floating point variable to store the well-depth parameter  (in kcal/mol)
	of the atom in the Lennard-Jones potential*/
	double lj_epsilon;
	// String variable to store the atom type name
	std::string atom_name;
	// String variable to store the residue type name
	std::string residue_name;
	// String variable to store the chain ID of the atom
	std::string chain_id;
};

typedef struct atom_struct struct_atom;

/**
 * @brief Data structure for an all-atom solute model
 *
 *
 */

struct solute_allatom_struct{

	// Integer variable to store the number of atoms of the solute molecule
	int num_atoms;
	// String variable to store the name of the solute molecule
	std::string molecule_name;
	// Array of structs containing information about the atoms in the solute
	struct_atom atoms[MAX_NUM_SOLUTE_ATOMS];
	// Floating point variable to store the total charge on the solute in e units
	double total_charge;
	// Floating point variables to store the minimum and maximum value of x among the (x,y,z) atom coordinates
	double min_x, max_x;
	// Floating point variables  to store the minimum and maximum value of y among the (x,y,z) atom coordinates
	double min_y, max_y;
	// Floating point variables to store the minimum and maximum value of z among the (x,y,z) atom coordinates
	double min_z, max_z;
};
typedef struct solute_allatom_struct struct_solute_allatom;


/**
 * @brief Data structure for solute's electrostatic potential
 *
 *
 */
struct solute_allatom_potentialgrid_struct{

	// Integer variables to store the number of grid points in x, y, and z dimensions
	int nx,ny,nz;
	// Floating point variables  to store the grid spacings along x, y, and z dimensions
	double hx,hy,hz;
	/* String variable to store the name of the file (in .dx format)
	containing the potential values at each grid point */
	std::string dxmap_filename;
	// 3-D array variable to store the potential values
	double ***data;
	// 1-D array variable of floating point values to store the grid point positions along x-axis
	double *xgrid;
	// 1-D array variable of floating point values to store the grid point positions along y-axis
	double *ygrid;
	// 1-D array variable of floating point values to store the grid point positions along z-axis
	double *zgrid;

};
typedef struct solute_allatom_potentialgrid_struct struct_allatom_potentialgrid;


#endif /* SOLUTE_STRUCTS_HPP_ */
