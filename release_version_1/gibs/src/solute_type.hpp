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
 * solute_type.hpp
 *
 *  Created on: Aug 12, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @brief Class definitions and function prototypes for solute model
 *
 *  @file solute_type.hpp
 *  @author Dennis G. Thomas
 *
 */

#ifndef SOLUTE_TYPE_HPP_
#define SOLUTE_TYPE_HPP_

#include "input_parameters.hpp"
#include "enum_types.hpp"
#include "solute_structs.hpp"
#include "box.hpp"


/**
 * @brief Class for storing and accessing the parameters of the solute model
 *
 * Use defined type CSoluteModel_t for class declarations
 */

class CSoluteModel{

	enum_solute_models solute_model;
	struct_solute_allatom *solute_allatom;
	struct_allatom_potentialgrid *potential_grid;

public:

	// constructor
	CSoluteModel(const CInputParameters_t &parameters,const CBox_t &box);

	// destructor
	~CSoluteModel(){
		switch(solute_model){

		case none:

			break;

		case all_atom:
			delete solute_allatom;
			solute_allatom = NULL;


			for (int ind1 = 0; ind1 < potential_grid->nx; ind1++) {

				for (int ind2 = 0; ind2 < potential_grid->ny; ind2++) {
					delete[] potential_grid->data[ind1][ind2];
				}
				delete[] potential_grid->data[ind1];
			}
			delete[] potential_grid->data;
			potential_grid->data = NULL;

			delete potential_grid;
			potential_grid = NULL;

			break;

		case cylindrical_polyion:

			break;

		case spherical_macroion:

			break;

		}

	}

	// read and set methods
	void readSoluteAllAtomCoordsFromFile(std::string file_name,
			enum_molecule_coordinate_file_type file_type); // for solute all atom model

	void readSolutePotentialFromFile(const CBox_t &box,double temp);

	// get methods
	int getSoluteAllAtom_NumAtoms() const {return solute_allatom->num_atoms;}
	std::string getSoluteAllAtom_MoleculeName() const {return solute_allatom->molecule_name;}
	double getSoluteAllAtom_TotalCharge() const {return solute_allatom->total_charge;}
	double getSoluteAllAtom_MinX() const {return solute_allatom->min_x;}
	double getSoluteAllAtom_MinY() const {return solute_allatom->min_y;}
	double getSoluteAllAtom_MinZ() const {return solute_allatom->min_z;}
	double getSoluteAllAtom_MaxX() const {return solute_allatom->max_x;}
	double getSoluteAllAtom_MaxY() const {return solute_allatom->max_y;}
	double getSoluteAllAtom_MaxZ() const {return solute_allatom->max_z;}
	double getSoluteAllAtom_AtomCoordX(int i) const {return solute_allatom->atoms[i].x;}
	double getSoluteAllAtom_AtomCoordY(int i) const {return solute_allatom->atoms[i].y;}
	double getSoluteAllAtom_AtomCoordZ(int i) const {return solute_allatom->atoms[i].z;}
	double getSoluteAllAtom_AtomRadius(int i) const {return solute_allatom->atoms[i].radius;}

	double getSoluteAllAtom_Potential(double x, double y, double z,const CBox_t box) const;

};
typedef class CSoluteModel CSoluteModel_t;

#endif /* SOLUTE_TYPE_HPP_ */
