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
 * write_results.cpp
 *
 *  Created on: Jul 9, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file write_results.cpp
 *  @author Dennis G. Thomas
 *
 * @brief Write functions defined for writing results to a file
 *
 */


#include "write_results.hpp"

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Writes the results of each particle type count to a file
 *
 * @param[in]	istep1(int)	Simulation step number
 * @param[in]	step1max(int)	Maximum or total number of simulation steps
 * @param[in]	particle_types(CParticleType_t) Parameters and data of all particle types
 *
 * @return void
 */
void writeResults_ParticleTypeCount(int istep1, int step1max,const CParticleType_t &particle_types) {
	int nparticle_types = particle_types.getNumParticleTypes();

	std::ofstream writefile[nparticle_types];

	int ind;

	std::string filename;
	const char* fname = "";

	for (ind = 0; ind < nparticle_types; ind++) {


		filename = "./outputfiles/particlecount_";

		filename.append(particle_types.getParticleType(ind).label);
		filename.append(".out");
		fname = filename.c_str();
		if(istep1>-1){
		writefile[ind].open(fname, std::ios::app);
		}else{
			writefile[ind].open(fname);
		}
		writefile[ind] << istep1 + 1 << '\t' << particle_types.getNum(ind) << '\t' << particle_types.getNumCavities(ind) <<
				'\t' << particle_types.getMuEx(ind) << '\t' << particle_types.getBip1(ind) << '\t' << particle_types.getEpsi(ind)
				<< '\t' << particle_types.getMuOld(ind) << std::endl;


		if (istep1 == step1max - 1) {
			writefile[ind].close();
		}
	}

}
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Writes the acceptance rates to a file
 *
 * @param[in]	istep1(int) 		Simulation step number
 * @param[in]	nsteps(int)			Total number of simulation steps
 * @param[in]	accept_count(int)	Number of times a move was accepted
 * @param[in]	attempt_count(int)	Number of times a move was attempted
 * @param[in]	fname(char*)		File name to write out the acceptance rates
 *
 * @return void
 */

void writeResults_AcceptanceRate(int istep1, int nsteps, int accept_count, int attempt_count, const char* fname) {
	std::ofstream writefile;

	double rate = (double) accept_count/attempt_count;

	if(istep1>1){
	writefile.open(fname, std::ios::app);
	}else{
		writefile.open(fname);
	}
	writefile << istep1  + 1<< '\t' << attempt_count << '\t' << accept_count << '\t' << rate << std::endl;

	if (istep1 == nsteps - 1) {
		writefile.close();
	}

}

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Writes out the total energy of the system to a file, at a given simulation step
 *
 * @param[in]	istep1(int)	Simulation step number
 * @param[in]	step1max(int)	Maximum or total number of simulation steps
 * @param[in]	ener(double)		Energy of the system in kcal/mol
 *
 * @return void
 *
 */
void writeResults_SystemEnergy(int istep1, int step1max, double ener) {
	std::ofstream writefile;

	std::string filename;
	const char* fname = "";

	filename = "./outputfiles/energy.out";
	fname = filename.c_str();

	if(istep1 > -1){
	writefile.open(fname, std::ios::app);
	}else{
		writefile.open(fname);
	}
	writefile << istep1 + 1 << '\t' << ener << std::endl;

	if (istep1 == step1max - 1) {
		writefile.close();
	}
}
