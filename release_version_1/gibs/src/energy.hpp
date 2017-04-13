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
