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
 * writeStateToFile.cpp
 *
 *  Created on: Jul 1, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file writeStateToFile.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief CParticleState class member function definitions to write out the coordinates and
 *  		type of particles in the simulation box
 *
 */

#include "particle_state.hpp"
/**
 * @brief  Writes out the coordinates, charge, and particle type label of all the particles in
 * 			the simulation box, at a given simulation step
 *
 * @param[in]	istep(int)	Simulation step number
 * @param[in]	file_index(int)	Index number to append to the file
 * @param[in]	box(CBox)		Siimulation box parameters
 * @param[in]	energy(double)	Energy  of the system in kcal/mol
 * @param[in]	filename(string)	Name of the file to write out the state of the system]
 *
 * @return void
 */
void CParticleState::writeToFile(int istep, int file_index,const CBox &box,double energy,std::string filename) {
	std::ofstream writefile;
	std::string sindex;
	std::stringstream sout;

	const char* fname = "";

	int ind;

	double dist;

	double xij, yij, zij;

	sout << file_index + 1;
	sindex = sout.str();

	filename.append(sindex);
	fname = filename.c_str();

	writefile.open(fname);
	writefile << "system_state_at_step = " << istep + 1 << std::endl;
	writefile << "number_of_particles_in_the_system = " << npart << std::endl;
	writefile << "Energy	" << energy << std::endl;
	writefile << "distance_from_box_center" << '\t' << "x" << '\t' << "y"
			<< '\t' << "z" << '\t' << "radius" << '\t' << "charge_(e_units)"
			<< '\t' << "particle_type_index" << '\t' << "particle_type_label"
			<< std::endl;
	for (ind = 0; ind < npart; ind++) {
		xij = particles[ind].x - box.getBoxXCenter();
		yij = particles[ind].y - box.getBoxYCenter();
		zij = particles[ind].z - box.getBoxZCenter();
		dist = sqrt(xij * xij + yij * yij + zij * zij);
		writefile << dist << '\t' << particles[ind].x << '\t' << particles[ind].y << '\t'
				<< particles[ind].z << '\t' << particles[ind].radius << '\t'
				<< particles[ind].charge << '\t' << particles[ind].ptype << '\t'
				<< particles[ind].label << std::endl;
	}
	writefile.close();
}


