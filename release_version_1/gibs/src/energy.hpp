/*
 * energy.hpp
 *
 *  Created on: Jul 15, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file energy.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief 	Prototypes for functions defined in energy.cpp
 */

#ifndef ENERGY_HPP_
#define ENERGY_HPP_

#include "generic.hpp"
#include "particle_type.hpp"
#include "particle_state.hpp"
#include "particlepair_type.hpp"
#include "cavity_grid.hpp"
#include "box.hpp"
#include "splines.hpp"

#include "solute_type.hpp"
#include "mwater.hpp"

/**
 * @brief Function to calculate the total interaction energy between particles of the same type
 */
double calcParticlePairEnergy_SameType(bool state_change,bool use_spm,bool use_mwater,const CParticleType_t &particle_types,
		const CParticleState_t &particle_state, const CParticlePairType_t &particlepair_types, CMWater &mwater,
		const CCavityGrid_t &cavity_grid,const CBox_t &box,
		double solvdiel, int ptclindex, int ptcltype);

/**
 * @brief Function to calculate the total interaction energy between particles of different types
 */
double calcParticlePairEnergy_DifferentType(bool state_change,bool use_spm,bool solvent_ion_att,const
		CParticleType_t &particle_types,const CParticleState_t	&particle_state, const
		CParticlePairType_t &particlepair_types,const CCavityGrid_t &cavity_grid,const CBox_t &box,double solvdiel, int ptclindex,int ptcltype);


/**
 * @brief Function to calculate the total interaction energy between particles
 */
double calcParticlePairEnergy(bool state_change,bool use_spm,bool use_mwater,bool solvent_ion_att,
		double solvdiel,const CParticleState_t &particle_state,const CParticlePairType_t &particlepair_types,
		const CParticleType_t &particle_types, CMWater &mwater,const CCavityGrid_t & cavity_grid, const CBox_t &box,int ptclindex,int ptcltype);

/**
 * @brief Function to calculate the interaction energy between solute and particles
 */
double calcSoluteParticleEnergy(bool state_change,bool use_spm,const CParticleState_t &particle_state,
		const CParticleType_t &particle_types,const CSoluteModel_t &solute,const CBox_t &box,int ptclindex);


#endif /* ENERGY_HPP_ */
