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
 * ion_type.hpp
 *
 *  Created on: Apr 13, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file ion_type.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains class declaration and function prototypes/definitions for ion types
 */

#ifndef ION_TYPE_HPP_
#define ION_TYPE_HPP_


#include "generic.hpp"
#include "input_parameters.hpp"
#include "iontype_structs.hpp"
#include "convert_functions.hpp"

/**
 * @brief Class for storing parameters and chemical
 * potential values of each ion type
 *
 * Use defined type CIonType_t for class declarations
 */
class CIonType{

	// Vector variable of type struct_iontype
	std::vector<struct_iontype> iontypes;

public:
	// constructor
	CIonType(int n){iontypes.resize(n);}

	// destructor
	~CIonType(){}

/**
 * @brief Function to add parameters for each ion type
 *
 */
void addIonType(int i,const CInputParameters_t &parameters,double box_vol);

/**
 * @brief Function to update the number of an ion species during simualation
 *
 */
void updateIonTypeNum(int i,int num){iontypes[i].num = num;}

// get methods
struct_iontype getIonType(int i) const {return iontypes[i];}

/**
 * @brief Get the ion concentration
 *
 * @param[in] i	The index (integer) corresponding to the position of
 * 					an iontype in the iontypes vector
 *
 * @return The concentration (double) of the ion species in Molar units
 */
double getIonConc(int i) const{return iontypes[i].conc;}

};

typedef class CIonType CIonType_t;

#endif /* ION_TYPE_HPP_ */
