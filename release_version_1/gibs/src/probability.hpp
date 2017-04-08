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
