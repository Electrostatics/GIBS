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
 *  computing [I]on distributions around [B]iomolecules with molecular [S]olvent models.
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
 * @return  zero(int)
 */
int main(){

	double time_elapsed; // elapsed time in hours
	time_t start,end;	// start and end time
	struct tm * timeinfo;
	time (&start);
	timeinfo = localtime(&start);

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

