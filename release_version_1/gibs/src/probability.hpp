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
 * probability.hpp
 *
 *  Created on: Jul 15, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file probability.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains function prototypes for computing accepting probabilities
 *  		and associated quantities (e.g., Boltzmann factor)
 */

#ifndef PROBABILITY_HPP_
#define PROBABILITY_HPP_

#include "generic.hpp"

/**
 * @brief Function to compute the Boltzmann factor
 */
double boltzmannFactor(double ener_chg,double temp);

/**
 *
 * @brief Function to compute the acceptance probability for inserting
 * 			a particle
 */
double insertProb(double cavity_prob,int N, double delU, double mu_chpot,
		double temp, double vol);

/**
 * @brief Function to compute the acceptance probability for deleting a
 * 			particle
 */
double deleteProb(double cavity_prob, int N, double delU, double mu_chpot,
		double temp, double vol);
/**
 * @brief Function to calculate the acceptance probability for displacing
 * 			a particle
 */
double dispProb(double ener_chg, double temp);

#endif /* PROBABILITY_HPP_ */
