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
 * solvent_type.cpp
 *
 *  Created on: Aug 5, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file solvent_type.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Definitions for member functions of class CSolvent
 *
 */

#include "solvent_type.hpp"

/*
 * @brief CSolvent constructor
 *
 * Sets default values for the solvent parameters
 *
 * @param[in] box_volume(double) Box volume in Angstrom cube units
 *
 */
CSolvent::CSolvent(double box_vol){
	charge = 0.0; //
	radius = 1.4; // radius in Angstrom
	label = "Water";
	packing_fraction = 0.001;
	dielectric = 78.5;
	/*target_num = (int) (box_vol*packing_fraction*3/
			(PI*4*radius*radius*radius));*/
	target_num = (int) (box_vol*packing_fraction/
				(pow(2.0*radius,3)));
	conc = numToConcMolarity(target_num,box_vol);
	mu_ex = 0.0;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Sets the solvent parameters with values specified by the user in
 * 			in the input parameter file
 *
 * Total charge is always set to zero
 *
 * @param[in]	parameters(CInputParameters_t)	Input Parameters specified by user
 * @param[in]	box_vol(double)		Box volume in Angstrom cube
 *
 * @return void
 */
void CSolvent::setParameters(const CInputParameters_t &parameters,double box_vol){

	radius = parameters.solvent_radius;
	label = parameters.solvent_label;
	packing_fraction = parameters.solvent_packing_fraction;
	dielectric = parameters.solvent_dielectric;
	charge = 0;
	/*target_num = (int) (box_vol*packing_fraction*6/
				(PI*8*pow(radius,3)));*/
	target_num = (int) (box_vol*packing_fraction/
				(pow(2.0*radius,3)));
	conc = numToConcMolarity(target_num,box_vol);
	mu_ex = parameters.solvent_mu_ex;

}
