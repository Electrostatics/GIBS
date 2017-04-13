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
/**
 * solvent_type.hpp
 *
 *  Created on: Apr 13, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file solvent_type.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Class declarations and function prototypes/definitions
 *  		for solvent type
 */

#ifndef SOLVENT_TYPE_HPP_
#define SOLVENT_TYPE_HPP_

#include "generic.hpp"
#include "constants.hpp"
#include "input_parameters.hpp"
#include "convert_functions.hpp"

/**
 * @brief Class declaration for storing solvent parameters
 *
 */
class CSolvent{

	/* Floating point variable to store the total charge (in e units) on
	 * solvent molecule, if represented as hard spheres
	 */
	double charge;
	/* Floating point variable to store the radius (in Angstrom) of
	 * solvent hard sphere
	 */
	double radius;
	/*
	 * String variable to store the label used for
	 * identifying a particle as solvent
	 */
	std::string label;
	/*
	 * Floating point variable to store the solvent packing fraction
	 */
	double packing_fraction;
	/*
	 * Floating point variable to store the solvent dielectric constant
	 */
	double dielectric;
	/*
	 * Integer variable to store the target number of solvent hard
	 * spheres in SPM simulations
	 */
	int target_num;
	/*
	 * Floating point variable to store the solvent concentration (M)
	 */
	double conc;
	/*
	 * Floating point variable to store the excess chemical potential
	 * of solvent hard spheres.
	 */
	double mu_ex;

public:
	// constructor
	CSolvent(double box_vol);

	// destructor
	~CSolvent(){}

	void setParameters(const CInputParameters_t &parameters,double box_vol);

	double getCharge() const {return charge;}
	double getRadius()const {return radius;}
	std::string getLabel() const {return label;}
	double getPackingFraction()const {return packing_fraction;}
	double getDielectricConstant() const {return dielectric;}
	double getTargetNumber() const {return target_num;}
	double getExcessChemicalPotential()const {return mu_ex;}
	double getConcentration() const {return conc;}// in Molar units
};





typedef class CSolvent CSolvent_t;

#endif /* SOLVENT_TYPE_HPP_ */
