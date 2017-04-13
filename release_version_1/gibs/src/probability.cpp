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
 * probability.cpp
 *
 *  Created on: Jul 15, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file probability.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains function definitions for calculating acceptance
 *  		probabilities and associated quantities
 *
 */

#include "probability.hpp"

/**
 * @brief Computes acceptance probability of displacing a particle
 *
 * @param[in]	ener_chg (double)	Change in energy in kcal/mol
 * @param[in]	temp (double)		Temperature in Kelvin
 *
 * @return The displacement probability (double)
 *
 */
double dispProb(double ener_chg, double temp) {
	double prob;

	prob = boltzmannFactor(ener_chg, temp);

	return (prob);
}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
/**
 * @brief Computes the acceptance probability of deleting a particle
 *
 * @param[in]	cavity_prob(double)	Probability of finding a cavity in a configuration of N-1 particles
 * @param[in]	N(int)		Number of particles before deletion
 * @param[in]	delU(double)	Change in energy in kcal/mol
 * @param[in]	mu_chpot(double)	Chemical potential of the particle in kcal/mol
 * @param[in]	temp(double)	Temperature in Kelvin unit
 * @param[in]	vol(double)		Box volume  in Angstrom cube
 *
 * @return The deletion probability (double)
 *
 */
double deleteProb(double cavity_prob, int N, double delU, double mu_chpot,
		double temp, double vol) {
	double prob, bfactor, ener;

	/* Acceptance probability for deleting a particle:
	 * [N!/(N-1)! ] *[1/f_N-1] *exp(-delU-kTln<n>-mu_ex)/kT
	 * f_N-1: probability of finding a cavity in the configuration of N-1 particles (cavity_prob)
	 *  N: number of particles before deletion
	 *  <n>:  concentration * volume (average number of particles in bulk)
	 */

	ener = delU + mu_chpot;
	bfactor = boltzmannFactor(ener, temp); // exp(-ener/(kB*temp))

	prob = bfactor * N / (6.022045 * 0.0001 * vol * cavity_prob);

	if (prob > 1.0) {
		prob = 1.0;
	}

	return (prob);
}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/**
 * @brief Computes the acceptance probability of inserting a particle
 *
 * @param[in]	cavity_prob(double)	Probability of selecting a cavity in a configuration of N particles
 * @param[in]	N(int)	Number of particles before insertion
 * @param[in]	delU(double)	Change in energy in kcal/mol
 * @param[in]	mu_chpot(double)	Chemical potential of the particle in kcal/mol
 * @param[in]	temp(double)	Temperature in Kelvin unit
 * @param[in]	vol(double)		Box volume in Angstrom cube
 *
 * @return The insertion probability (double)
 *
 */

double insertProb(double cavity_prob, int N, double delU, double mu_chpot,
		double temp, double vol) {

	double prob, bfactor, ener;

	/* Acceptance probability for inserting a particle:
	 * [N!/(N+1)! ] *[f_N] *exp(-delU+kTln<n>+mu_ex)/kT
	 * f_N: probability of finding a cavity in the configuration of N particles (cavity_prob)
	 *  N: number of particles before insertion
	 *  <n>:  concentration * volume (average number of particles in bulk)
	 */


	ener = delU - mu_chpot;
	bfactor = boltzmannFactor(ener, temp); // exp(-ener/(kB*temp))

	prob = bfactor * (vol * 6.022045 * 0.0001)*cavity_prob / (N + 1);


	if (prob > 1.0) {
		prob = 1.0;
	}

	return (prob);
}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
/**
 * @brief Computes the Boltzmann factor
 *
 ^ @param[in]	ener_chg(double)	Change in energy  in kcal/mol
 * @param[in]	temp(double)		Temperature in Kelvin
 *
 * @return The Boltzmann factor (double)
 *
 */
double boltzmannFactor(double ener_chg, double temp) {
	double bzfactor, beta, expnt;

	// beta = 1.0/(k_B*avognum*convJtoCal*0.001*temp) = 1.0/(0.001985862*temp)

	//beta = 1.0/(k_B*avognum*convJtoCal*0.001*temp);
	// beta = 1.0/(k_B*avognum*convJtoCal*0.001*temp) = 1.0/(0.00198719137*temp)
	beta = 503.222797 / temp;
	// kBT in kcal/mol

	//cout <<"[probability.cpp: boltzmannFactor] beta (1/(kcal/mol)) = " << beta << endl;
	expnt = ener_chg * beta;

	//cout <<"[probability.cpp: boltzmannFactor] expnt = " << expnt << endl;

	bzfactor = exp(-expnt);
	if (-expnt >= pow(10, 45)) {
		expnt = -pow(10, 45);
		bzfactor = exp(-expnt);
	}

	if (-expnt <= -pow(10, 45)) {
		expnt = pow(10, 45);
		bzfactor = 0;
	}

	return (bzfactor);
}

