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
