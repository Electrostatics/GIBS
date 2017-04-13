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
 * mwater.hpp
 *
 *  Created on: Mar 1, 2016
 *      Author:  Dennis G. Thomas
 *
 *  @file mwater.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Class definitions and parameters of Monoatomic water (MWater) model
 *
 *
 * Reference for MWater model:
 * 	Valeria Molinero* and Emily B. Moore, "Water Modeled As an Intermediate Element
 * 	between Carbon and Silicon", J. Phys. Chem. B, 2009, 113 (13), pp 4008–4016
 */

#ifndef MWATER_HPP_
#define MWATER_HPP_


#include "generic.hpp"


class CMWater{

public:
	double A;
	double B;
	int p;
	int q;
	double gamma;
	double a;
	double costheta0;
	double lamda;
	double epsilon;	// kcal/mol
	double sigma;	// Angstrom
	double asigma; // asigma = a x sigma
	double AepsilonB; // AepsilonB = A x epsilon x B
	double Aepsilon; // Aepsilon = A x epsilon
	double lamdaepsilon; // lamdaepsilon = lamda x epsilon
	double gammasigma; // gammasigma  = gamma x sigma

	std::vector<int> solvent_indices;

	CMWater(){

		this->A = 7.049556277;
		this->B = 0.6022245584;
		this->p = 4;
		this->q = 0;
		this->gamma = 1.2;
		this->a = 1.8;
		this->costheta0 =-0.333333333333; // theta0 = 109.47; costheta0 = -1/3.
		this->lamda = 23.15;
		this->epsilon = 6.189;   // units: kcal/mol
		this->sigma = 2.3925;
		this->asigma = a * sigma;
		this->AepsilonB = A * epsilon * B;
		this->Aepsilon = A * epsilon;
		this->lamdaepsilon = lamda * epsilon;
		this->gammasigma = gamma * sigma;

		this->solvent_indices.resize(1000);


	}

	~CMWater(){}

};


#endif /* MWATER_HPP_ */
