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
 * ion_type.cpp
 *
 *  Created on: Aug 5, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file ion_type.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains definitions for member functions of class CIonType
 *
 */


#include "ion_type.hpp"

/**
 * @brief Adds parameter values for each ion type
 *
 * @param[in]		i(int) 		Index of an element of iontypes struct vector
 * @param[in,out]	parameters(CInputParameters_t)	Input parameters
 * @param[in]		box_volume(double)	Box volume in Angstrom cube unit
 *
 * @return void
 */
void CIonType::addIonType(int i,const CInputParameters_t &parameters,double box_vol){
	struct_iontype iont;

	iont.radius = parameters.iontypes[i].radius;
			iont.conc = parameters.iontypes[i].conc;
			iont.num = 0;//parameters.iontypes[i].num;
			iont.charge = parameters.iontypes[i].charge;
			iont.itype = parameters.iontypes[i].itype;
			iont.mu_ex = parameters.iontypes[i].mu_ex;
			iont.label = parameters.iontypes[i].label;
			//iont.lj_epsilon = parameters.iontypes[i].lj_epsilon;
			//iont.lj_sigma = parameters.iontypes[i].lj_sigma;


			iont.mu_old = kTLnConc(iont.conc, parameters.temperature) + iont.mu_ex;
						iont.mu_new = iont.mu_old;
						iont.nion_avg = 0;//concMolarityToNum(iont.conc, box_vol);
						iont.cion_avg = 0;//iont.conc;

	this->iontypes[i]=iont;
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
