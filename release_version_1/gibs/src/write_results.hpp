/*
 * write_results.hpp
 *
 *  Created on: Jul 9, 2015
 *      Author: Dennis G. Thomas
 *
 *
 *  @file write_results.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Prototypes of functions defined in write_results.cpp
 *
 */

#ifndef WRITE_RESULTS_HPP_
#define WRITE_RESULTS_HPP_

#include "generic.hpp"
#include "particle_type.hpp"

/**
 * @brief Function to write out the results of each particle type count to a file, at each simulation step
 */
void writeResults_ParticleTypeCount(int istep1, int step1max,const CParticleType_t &particle_types);

/**
 * @brief Function to write out the acceptance rates to a file, at each simulation step
 */
void writeResults_AcceptanceRate(int istep1, int nsteps, int accept_count, int attempt_count, const char* fname);

/**
 * @brief Function to write out the total energy of the system to a file at each simulation step
 */
void writeResults_SystemEnergy(int istep1, int step1max, double ener) ;


#endif /* WRITE_RESULTS_HPP_ */
