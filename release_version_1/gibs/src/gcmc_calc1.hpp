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
 * gcmc_calc1.hpp
 *
 *  Created on: Jul 15, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file gcmc_calc1.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Prototypes of functions defined in gcmc_calc1.cpp
 *
 */

#ifndef GCMC_CALC1_HPP_
#define GCMC_CALC1_HPP_

#include "generic.hpp"
#include "probability.hpp"
#include "energy.hpp"
#include "input_parameters.hpp"
#include "particle_state.hpp"
#include "box.hpp"
#include "cavity_grid.hpp"
#include "rdf_methods.hpp"
#include "particle_type.hpp"
#include "particlepair_type.hpp"
#include "write_results.hpp"
#include "convert_functions.hpp"
#include "solute_type.hpp"

#include "mwater.hpp"

/**
 * @brief Function to run the steps of the GCMC simulation
 */
void gcmcCalc1(const CInputParameters_t &parameters, const CBox_t &box,
		CRdf_t &rdf, CParticleType_t &particle_types,
		CParticlePairType_t &particlepair_types, CMWater &mwater,
		CParticleState_t &particle_state, CCavityGrid_t &cavity_grid,
		const CSoluteModel_t &solute, double &energy_old);

/**
 * @brief Function to randomly insert a particle type
 */
void insertParticle_1(const CInputParameters_t &parameters, const CBox_t &box,
		CParticleState &particle_state, CParticleType_t &particle_types,
		CParticlePairType_t &particlepair_types, CCavityGrid &cavity_grid,
		 CMWater &mwater,const CSoluteModel_t &solute, double &energy_old,
		int &num_accept_insert, double rand_cavity, double randpos_x,
		double randpos_y, double randpos_z, double rand_ins,
		double rand_ptcltype,bool calc_statistics,int &state_count, CRdf_t &rdf);

/**
 * @brief Function to randomly select and delete a particle type
 */
void deleteParticle_1(const CInputParameters_t &parameters, const CBox_t &box,
		CParticleState_t &particle_state, CParticleType_t &particle_types,
		CParticlePairType_t &particlepair_types, CCavityGrid_t &cavity_grid,
		 CMWater &mwater,const CSoluteModel_t &solute, double &energy_old,
		int &num_accept_delete, double rand_countindex, double rand_del,
		double rand_ptcltype,bool calc_statistics,int &state_count, CRdf_t &rdf);

/**
 * @brief Function to randomly select and displace a particle type
 */
void displaceParticle_1(const CInputParameters_t &parameters, const CBox_t &box,
		CParticleState_t &particle_state, CParticleType_t &particle_types,
		CParticlePairType_t &particlepair_types, CCavityGrid_t &cavity_grid,
		 CMWater &mwater,const CSoluteModel_t &solute, double &energy_old, int &num_accept_disp,
		double rand_cavity, double randpos_x, double randpos_y,
		double randpos_z, double rand_disp, double rand_countindex,
		double rand_ptcltype,bool calc_statistics,int &state_count,CRdf &rdf);



#endif /* GCMC_CALC1_HPP_ */
