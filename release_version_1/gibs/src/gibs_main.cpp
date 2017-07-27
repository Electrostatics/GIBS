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
 * gibs_main.cpp
 *
 *  Created on: Apr 7, 2015
 *      Author:  Dennis G. Thomas
 *
 *	@file gibs_main.cpp
 *  @author Dennis G. Thomas
 *  @version 1.0 12/08/2016
 *  @date 12/08/2016
 *
 *  @brief This file has the main() function of the GIBS simulation program.
 *
 *  GIBS stands for [G]rand Canonical Monte Carlo (GCMC) simulation program for
 *  computing [I]on distributions around [B]iomolecules with hard sphere [S]olvent models.
 *
 *  GIBS can be used to perform the following types of simulations:
 *  1: [GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION] GCMC simulation for
 *  computing the excess chemical potential and mean activity coefficients
 *  of individual ions and salt.
 *	2: [GCMC_WITH_SOLUTE_ALLATOM_MODEL] GCMC simulation for computing ion
 *	distributions around a fixed solute represented in atomic detail.
 *
 */


#include "generic.hpp"
#include "RunGCMCSimulation1.hpp"

/**
 * @brief	The main program
 *
 * @return  zero(integer)
 */
int main(){

	double time_elapsed; // elapsed time in hours
	time_t start,end;	// start and end time
	struct tm * timeinfo;
	time (&start);
	timeinfo = localtime(&start);

	std::ios::sync_with_stdio(false);

	std::cout << "Using GIBS software, version 1.0" << std::endl;
	std::cout << "Simulation started on (local time): " << asctime(timeinfo) << std::endl;

	runSimulation_1(); // driver program

	time (&end);

	timeinfo = localtime(&end);
	std::cout << "Simulation ended on (local time): " << asctime(timeinfo) << std::endl;

	time_elapsed = difftime(end,start);
	time_elapsed = time_elapsed/3600.0;

	std::cout << "Total hours taken for the simulation to complete = " << time_elapsed << std::endl;
	std::cout << std::endl;
	std::cout << "Thank you for using GIBS!" << std::endl;

	return 0;
}

