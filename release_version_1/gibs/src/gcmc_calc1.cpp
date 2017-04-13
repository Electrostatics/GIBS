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
 * gcmc_calc1.cpp
 *
 *  Created on: Jul 15, 2015
 *      Author: Dennis G. Thomas
 *
 * @file gcmc_calc1.cpp
 * @author Dennis G. Thomas
 *
 * @brief Functions implementing the various steps of the GCMC algorithm
 *
 */

#include "gcmc_calc1.hpp"

/*
 * @brief Performs GCMC simulation
 *
 * @param[in]		parameters(CInputParameters_t)		Input parameters
 * @param[in]		box(CBox_t)					Simulation box parameters
 * @param[in,out]	rdf(CRdf_t)					Radial distribution functions and parameters
 * @param[in]		particle_types(CParticleType_t)		Parameters of all particle types
 * @param[in]		particlepair_types(CParticlePairType_t)	Parameters of all particle pair types
 * @param[in]		mwater(CMWater)			Parameters of MWater model
 * @param[in,out]	particle_state(CParticleState_t)	 Position and type data of all the particles in the system
 * @param[in]		cavity_grid(CavityGrid_t)		Cavity grid data
 * @param[in]		solute(CSoluteModel_t)		Parameters and data of the solute model
 * @param[in,out]	energy_old(double)		Current energy of the system (kcal/mol)
 *
 * @return void
 */
void gcmcCalc1(const CInputParameters_t &parameters, const CBox_t &box,
		CRdf_t &rdf, CParticleType_t &particle_types,
		CParticlePairType_t &particlepair_types,  CMWater &mwater,
		CParticleState_t &particle_state, CCavityGrid_t &cavity_grid,
		const CSoluteModel_t &solute, double &energy_old) {

	int record_step = (parameters.nsteps - parameters.neqlb_steps)/ parameters.num_record_states;
	int record_state_count = 0;
	int state_count = 0;
	int nparticle_types = particle_types.getNumParticleTypes();		// Number of particle types

	/*-------Calculate the energy of the initial state of the system-------*/

	/*for (int i =0; i< particle_state.getParticleNum(); i++){
		std::cout << " particle type = " << particle_state.getParticle(i).label << ", cavity cell index = " << particle_state.getParticle(i).gridcell_1Dindex
				<< " , LJ cell index = " << particle_state.getParticle(i).lj_gridcell_1Dindex << std::endl;
	}*/
	if (energy_old==0){
		energy_old = calcParticlePairEnergy(0, parameters.use_spm,parameters.use_mwater,
				parameters.solvent_ion_attraction, parameters.solvent_dielectric,
				particle_state, particlepair_types, particle_types, mwater,cavity_grid,box, 0, 0);


		if (parameters.solute_model == all_atom) {
			energy_old += calcSoluteParticleEnergy(false,parameters.use_spm,particle_state, particle_types, solute, box,
					0);

		}
		std::cout << "[gcmc_calc1.cpp (gcmcCalc1) energy (kcal/mol) = "
				<< energy_old << std::endl;

	}

	/* Integer variable to store the number of times an insertion attempt was accepted*/
	int num_accept_insert = 0;
	/* Integer variable to store the number of times a deletion attempt was accepted*/
	int num_accept_delete = 0;
	/* Integer variable to store the number of times a single
	 *  particle displacement attempt was accepted*/
	int num_accept_disp = 0;

	/* Integer variable to store the number of insertion attempts*/
	int num_insert_attempts = 0;
	/* Integer variable to store the number of deletion attempts*/
	int num_delete_attempts = 0;
	/* Integer variable to store the number of single particle displacement attempts*/
	int num_disp_attempts = 0;

	/* MTRand function to generate a  random floating point number [0,1) for
	 * randomly selecting one of the three GCMC moves: insertion, deletion, or displacement */
	MTRand rand_insdelmov;

	/* MTRand function to generate a random floating point number [0,1) for
	 * randomly selecting an insertion or a deletion attempt*/
	MTRand rand_insdel;

	/* MTRand function to generate a random floating point number [0,1) for
	 * randomly selecting a particle of a given type  */
	MTRand rand_countindex;

	/* MTRand function to generate a random floating point number [0,1) for
	 *  randomly selecting  an unoccupied grid cell within a cavity grid segment */
	MTRand rand_cavity;

	/* MTRand_open function to generate a random floating point number (0,1) for
	 * randomly selecting a x,y,z position for a particle.Grid insertion algorithm does
	 * not randomly select the particle's position, as it places the particle
	 * at the center of the randomly selected cavity grid cell.*/
	MTRand_open randpos;

	/* MTRand_closed function to generate a random floating point number [0,1] that is
	 * used to determine whether or not to accept a DISPLACEMENT attempt by comparing it to
	 * the acceptance probability  */
	MTRand_closed rand_disp;

	/*MTRand function to generate a random floating point number [0,1) for randomly
	 * selecting a particle type for insertion/deletion/displacement */
	MTRand rand_ptcltype;  // [0,1)

	/* MTRand_closed function to generate a random floating point number [0,1] that is
	 * used to determine whether or not to accept a DELETION attempt by comparing it to
	 * the acceptance probability  */
	MTRand_closed rand_del;

	/* MTRand_closed function to generate a random floating point number [0,1] that is
	 * used to determine whether or not to accept an INSERTION attempt by comparing it to
	 * the acceptance probability  */
	MTRand_closed rand_ins;

	/*Write out the initial count of each particle type.*/
	writeResults_ParticleTypeCount(-1, parameters.nsteps, particle_types);
	/* Write out the initial potential energy of the system.*/
	writeResults_SystemEnergy(-1, parameters.nsteps, energy_old);


	/*Boolean variable for calculating the statistics when the simulation step number exceed
	 * the number of steps for equilibration */
	bool calc_statistics = false;

	for (int istep1 = 0; istep1 < parameters.nsteps; istep1++) {

		if (istep1 >= parameters.neqlb_steps) {
			calc_statistics = true;
		}
		/*-------GCMC routine for calculating the chemical potential------- */

		/* Each simulation step consists of "parameters.gcmc_cycl" insertion/deletion cycles,
		 * and "parameters.disp_cycl" displacement steps. Each insertion/deletion cycle
		 * consists of an insertion step and a deletion step.
		 */
		if (parameters.simulation_type
				== GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION) {

			/*Insertion/deletion cycle */

			for (int istep2 = 0; istep2 < parameters.gcmc_cycl; istep2++) {

				/*single-particle insertion */
				double insert_rand_cavity = rand_cavity();
				double insert_randpos_x = 0.5; // randpos(); //0.5
				double insert_randpos_y = 0.5;//randpos(); // 0.5
				double insert_randpos_z = 0.5; //randpos(); //0.5
				double insert_rand_ins = rand_ins();
				double insert_rand_ptcltype = rand_ptcltype();

				insertParticle_1(parameters, box, particle_state,
						particle_types, particlepair_types, cavity_grid, mwater,solute,
						energy_old, num_accept_insert, insert_rand_cavity,
						insert_randpos_x, insert_randpos_y, insert_randpos_z,
						insert_rand_ins, insert_rand_ptcltype,calc_statistics,state_count,rdf);

				++num_insert_attempts;

				/*single-particle deletion*/

				double delete_rand_countindex = rand_countindex();
				double delete_rand_del = rand_del();
				double delete_rand_ptcltype = rand_ptcltype();

				deleteParticle_1(parameters, box, particle_state,
						particle_types, particlepair_types, cavity_grid, mwater, solute,
						energy_old, num_accept_delete, delete_rand_countindex,
						delete_rand_del, delete_rand_ptcltype,calc_statistics,state_count,rdf);

				++num_delete_attempts;
			}
			/*single-particle displacement*/
			for (int istep3 = 0; istep3 < parameters.disp_cycl; istep3++) {

				double disp_rand_cavity = rand_cavity();
				double disp_randpos_x = 0.5; //randpos();
				double disp_randpos_y = 0.5; //randpos();
				double disp_randpos_z = 0.5; //randpos();
				double disp_rand_disp = rand_disp();
				double disp_rand_countindex = rand_countindex();
				double disp_rand_ptcltype = rand_ptcltype();


				displaceParticle_1(parameters, box, particle_state,
						particle_types, particlepair_types, cavity_grid, mwater, solute,
						energy_old, num_accept_disp, disp_rand_cavity,
						disp_randpos_x, disp_randpos_y, disp_randpos_z,
						disp_rand_disp, disp_rand_countindex,
						disp_rand_ptcltype,calc_statistics,state_count,rdf);

				++num_disp_attempts;

			}
		}
		/*------- GCMC routine for computing the distribution of particles
		 * around and all-atom solute -------*/

		/* The  The parameters.gcmc_cycl and parameters.disp_cycl values are
		 * normally set such that 50-70% of the random moves (nseps) are
		 * single-particle displacements. In particular, the parameters.gcmc_cycl value is
		 * normally set equal to or less than the average number of total ions,
		 * expected to be present in the simulation box.  e.g..,
		 * if select_moveset > # particles in the box, then an insertion/deletion
		 * is attempted; otherwise, a single-particle displacement is attempted./
		 */

		else if (parameters.simulation_type
				== GCMC_WITH_SOLUTE_ALLATOM_MODEL) {

			int select_moveset = int(
					rand_insdelmov()
					* (parameters.gcmc_cycl + parameters.disp_cycl))+ 1;
			if (select_moveset > particle_state.getParticleNum()) { // select insertion / deletion
				// Insertion and deletion steps are randomly selected with equal probability
				if (rand_insdel() >= 0.5) {
					double insert_rand_cavity = rand_cavity();
					double insert_randpos_x = 0.5;//randpos();
					double insert_randpos_y = 0.5;//randpos();
					double insert_randpos_z = 0.5;//randpos();
					double insert_rand_ins = rand_ins();
					double insert_rand_ptcltype = rand_ptcltype();

					insertParticle_1(parameters, box, particle_state,
							particle_types, particlepair_types, cavity_grid, mwater,
							solute, energy_old, num_accept_insert,
							insert_rand_cavity, insert_randpos_x,
							insert_randpos_y, insert_randpos_z, insert_rand_ins,
							insert_rand_ptcltype,calc_statistics,state_count,rdf);

					++num_insert_attempts;
				} else {

					double delete_rand_countindex = rand_countindex();
					double delete_rand_del = rand_del();
					double delete_rand_ptcltype = rand_ptcltype();

					deleteParticle_1(parameters, box, particle_state,
							particle_types, particlepair_types, cavity_grid, mwater,
							solute, energy_old, num_accept_delete,
							delete_rand_countindex, delete_rand_del,
							delete_rand_ptcltype,calc_statistics,state_count,rdf);

					++num_delete_attempts;
				}
			} else { // single-particle displacement

				double disp_rand_cavity = rand_cavity();
				double disp_randpos_x = 0.5;//randpos();
				double disp_randpos_y = 0.5;//randpos();
				double disp_randpos_z = 0.5;//randpos();
				double disp_rand_disp = rand_disp();
				double disp_rand_countindex = rand_countindex();
				double disp_rand_ptcltype = rand_ptcltype();

				displaceParticle_1(parameters, box, particle_state,
						particle_types, particlepair_types, cavity_grid, mwater,solute,
						energy_old, num_accept_disp, disp_rand_cavity,
						disp_randpos_x, disp_randpos_y, disp_randpos_z,
						disp_rand_disp, disp_rand_countindex,
						disp_rand_ptcltype,calc_statistics,state_count,rdf);

				++num_disp_attempts;
			}
		}

		/*-------Record particle type count-------*/
		if (parameters.record_particletype_count == true) {
			if ((istep1 == 0)
					|| ((istep1 + 1) % parameters.record_every_n_steps == 0)) {
				writeResults_ParticleTypeCount(istep1, parameters.nsteps,
						particle_types);
			}
		}
		/*-------Record acceptance rates-------*/
		if (parameters.record_acceptance_rate == true) {
			if ((istep1 == 0)
					|| ((istep1 + 1) % parameters.record_every_n_steps == 0)) {
				writeResults_AcceptanceRate(istep1, parameters.nsteps,
						num_accept_insert, num_insert_attempts,
						"./outputfiles/insertion_acceptance_rate");

				writeResults_AcceptanceRate(istep1, parameters.nsteps,
						num_accept_delete, num_delete_attempts,
						"./outputfiles/deletion_acceptance_rate");

				writeResults_AcceptanceRate(istep1, parameters.nsteps,
						num_accept_disp, num_disp_attempts,
						"./outputfiles/displacement_acceptance_rate");

				num_accept_insert = 0;
				num_accept_delete = 0;
				num_accept_disp = 0;
				num_insert_attempts = 0;
				num_delete_attempts = 0;
				num_disp_attempts = 0;

			}
		}

		/*-------Record system energy -------*/
		if (parameters.record_system_energy == true) {
			if ((istep1 == 0)
					|| ((istep1 + 1) % parameters.record_every_n_steps == 0)) {
				writeResults_SystemEnergy(istep1, parameters.nsteps,
						energy_old);
			}
		}

		/*-------Record state-------*/
		if (istep1
				== -1 + parameters.neqlb_steps
				+ record_state_count * record_step) {

			particle_state.writeToFile(istep1, istep1, box,energy_old,
					"./outputfiles/state_");
			++record_state_count;
		}

		/* Obtain statistics for calculating the average number of each particle type
		 * and RDFs of each particle pair*/
		//if ((istep1 >= parameters.neqlb_steps) && ((istep1 - parameters.neqlb_steps) / 100 == state_count)) {
		/*	if (istep1 >= parameters.neqlb_steps) {
			++state_count;
			for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
				double num_avg = particle_types.getParticleType(ind1).num_avg;
				num_avg += particle_types.getParticleType(ind1).num;
				particle_types.updateNumAvg(ind1, num_avg);

			}

			if (parameters.calculate_particle_pair_rdf == true) {
				rdf.calcParticlePairBinCounts(particle_types, particle_state);
				rdf.incrementParticlePairBinAvgCounts();
			}

			if (parameters.calculate_solute_particle_rdf == true) {
				rdf.calcSoluteParticleBinCounts(particle_types, particle_state);
				rdf.incrementSoluteParticleBinAvgCounts();
			}
		}*/
	} //<----------end of istep1 for loop

	/* -------Calculate average number and concentration of each particle type ------- */
	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
		double num_avg = particle_types.getParticleType(ind1).num_avg
				/ double(state_count);
		particle_types.updateNumAvg(ind1, num_avg);
		double conc_avg = numToConcMolarity(num_avg, box.getBoxVol());
		particle_types.updateConcAvg(ind1, conc_avg);
		std::cout <<"[gcmc_calc1.cpp(gcmcCalc1)]: particle type " << particle_types.getLabel(ind1)
																						<< ", average number = " << num_avg << std::endl;
	}

	/*------- Compute particle pair RDF-------*/
	if (parameters.calculate_particle_pair_rdf == true) {

		/*Calculate average bin counts */
		rdf.finalizeParticlePairBinAvgCounts(state_count);
		/*Compute particle pair RDFs */
		rdf.calcParticlePairRDF(particle_types);

		/*Write out the radial distribution function of each particle pair type */
		rdf.writeParticlePairRDFtoFile();
	}


	/*------- Compute solute particle RDF -------*/
	if (parameters.calculate_solute_particle_rdf == true) {
		/*Calculate average bin counts */
		rdf.finalizeSoluteParticleBinAvgCounts(state_count);

		/*Compute solute particle pair RDFs */
		rdf.calcSoluteParticleRDF(particle_types);
		/*Write out the radial distribution function of each solute- particle type */
		rdf.writeSoluteParticleRDFtoFile();

	}

	std::cout << "[gcmc_calc1.cpp (gcmcCalc1)] Number of states averaged over = "
			<< state_count << std::endl;

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Inserts a randomly selected particle type in a randomly selected cavity site.
 *
 * @param[in]		parameters(CInputParameters_t)		Input parameters
 * @param[in]		box(CBox_t)				Simulation box parameters
 * @param[in,out]	particle_state(CParticleState_t)	 Position and type data of all the particles in the system
 * @param[in,out]	particle_types(CParticleType_t)		Parameters of all particle types
 * @param[in]		particlepair_types(CParticlePairType_t)	Parameters of all particle pair types
 * @param[in]		cavity_grid(CavityGrid_t)		Cavity grid data
 * @param[in]		mwater(CMWater)			Parameters of MWater model
 * @param[in]		solute(CSoluteModel_t)			Parameters and data of the solute model
 * @param[in,out]	energy_old(double)		Current energy of the system (kcal/mol)
 * @param[in,out]	num_accept_insert(int)	Number of times an insertion attempt was accepted
 * @param[in]		rand_cavity(double)		Random number [0,1) for randomly selecting  an unoccupied grid cell within a cavity grid segment
 * @param[in]		randpos_x(double)		Random number (0,1) for randomly selection the x position a in a grid cell (currently a fixed value of 0.5 is used)
 * @param[in]		randpos_y(double)		Random number (0,1) for randomly selection the y position a in a grid cell (currently a fixed value of 0.5 is used)
 * @param[in]		randpos_z(double)		Random number (0,1) for randomly selection the z position a in a grid cell (currently a fixed value of 0.5 is used)
 * @param[in]		rand_ins(double)		Random number [0,1] to test the acceptance criterion for accepting an INSERTION attempt
 * @param[in]		rand_ptcltype(double)	Random number (0,1] for selecting the type of particle for insertion
 * @param[in]		calc_statistics(bool)	True or False option (boolean) to determine when to start computing the statistics
 * @param[in,out]	state_count(int)		Number of states used for averaging the statistical quantities
 * @param[in,out]	rdf(CRdf_t)				Radial distribution functions and parameters
 *
 * @return void
 */
void insertParticle_1(const CInputParameters_t &parameters, const CBox_t &box,
		CParticleState &particle_state, CParticleType_t &particle_types,
		CParticlePairType_t &particlepair_types, CCavityGrid &cavity_grid,  CMWater &mwater,
		const CSoluteModel_t &solute, double &energy_old,
		int &num_accept_insert, double rand_cavity, double randpos_x,
		double randpos_y, double randpos_z, double rand_ins,
		double rand_ptcltype,bool calc_statistics,int &state_count, CRdf_t &rdf) {


	/*Get the number of particle types */
	int nparticle_types = particle_types.getNumParticleTypes();
	/*Get the total number of cavity grid cells */
	int num_cells = cavity_grid.getCavityGridNumCells();

	/*Randomly select the type of particle to be inserted. The value will always be
	 * less than the number of particle types (nparticle_types)*/
	int ptcltype_index = int(rand_ptcltype * nparticle_types);

	/* Get the number of cavity grid cells available for inserting the center of a new
	 * particle of the randomly selected particle type ( ptcltype_index)*/
	int num_cavities = particle_types.getNumCavities(ptcltype_index);
	//for (int ptcltype_index = 0; ptcltype_index < nparticle_types; ptcltype_index++){

	if (num_cavities>= 1) { // only attempt to insert the particle if there is a cavity grid cell available.

		//	std::cout << "[gcmc_calc1.cpp (insertParticle_1) attempting to insert particle of type-------->, "<< particle_types[ptcltype_index].label << std::endl;
		// randomly select a cavity  from the cavity list of the selected particle
		//	MTRand rand_cavity;

		/*Declare a random number generator to randomly select a cavity grid segment */
		MTRand rand_cavity_segment;
		/* Initialize a boolean variable to false, which will later be set to True, if
		 * a cavity grid cell is found for inserting the particle*/
		bool found_cavityindex= false;
		/* Get the total number of cavity grid segments */
		int num_cavity_segments = cavity_grid.getNumCavitySegments();
		/*Get the total number of cavity grid cells in the selected segment */
		int num_cells_in_segment =cavity_grid.getNumCellsInSegment();

		int cellindex_1D;
		/*------- Find a cavity grid cell for inserting the particle -------*/
		while(found_cavityindex == false ){
			/* Select a segment*/
			int cavity_segmentindex = int(rand_cavity_segment()*num_cavity_segments);
			/*Get the number of cavities available for the current particle type in the segment.*/
			int num_cavities_in_segment = particle_types.getNumCavitiesInSegment(ptcltype_index,
					cavity_segmentindex);
			/*Select a cavity out of the list of all cavities in the segment for the current particle type.*/
			int seg_cavitycount_index = int(rand_cavity*num_cavities_in_segment);

			/*Get the first cell index of the segment*/
			int seg_first_cellindex = cavity_grid.getFirstCellIndexOfSegment(cavity_segmentindex);

			/*Get the last cell index of the segment*/
			int seg_last_cellindex = seg_first_cellindex + num_cells_in_segment-1;

			if(num_cavities_in_segment > 0){
				int it1 = seg_first_cellindex;

				int count_index = 0;
				while((it1 <= seg_last_cellindex) && (count_index != seg_cavitycount_index+1) ){
					cellindex_1D = particle_state.getCellIndexByCavityCountIndex(
							ptcltype_index, it1);
					if(cellindex_1D != -1){
						++count_index;
					}
					++it1;
				}
				found_cavityindex = true;
			}else{
				found_cavityindex = false;
			}
		}


		// Randomly select the index of a cavity site
		/*		int cavity_countindex = int(rand_cavity * num_cavities);

		// From the cavity count index, get the grid cell index (cellindex_1D) of the selected cavity
		int cellindex_1D;
		int it1=0;

		int count_num = 0;
		while ((it1 < num_cells) && (count_num != cavity_countindex+1)){
			cellindex_1D = particle_state.getCellIndexByCavityCountIndex(
					ptcltype_index, it1);

			if(cellindex_1D != -1){
				++count_num;
			}
			++it1;
		}*/



		/*int cellindex_1D = particle_state.getCellIndexByCavityCountIndex(
				ptcltype_index, cavity_countindex);*/

		/*-------Assign particle type parameters to the new particle-------*/
		struct_particle ptcl;
		particle_types.assignParticleTypeParametersToParticle(ptcl,
				ptcltype_index);

		/*-------Set the grid cell indices and the position for the particle in
		 * the selected grid cell cavity.-------*/
		struct_cavitygrid_cellindex_coords cellindex3D_coords;
		/*Get 3D grid cell indices*/
		cellindex3D_coords = cavity_grid.getCellIndex3D(cellindex_1D);
		//MTRand_open randpos;

		ptcl.x = cavity_grid.getX(cellindex3D_coords.icx)+ randpos_x * cavity_grid.getHx();
		ptcl.y = cavity_grid.getY(cellindex3D_coords.icy)+ randpos_y * cavity_grid.getHy();
		ptcl.z = cavity_grid.getZ(cellindex3D_coords.icz)+ randpos_z * cavity_grid.getHz();

		ptcl.gridcell_ix = cellindex3D_coords.icx;
		ptcl.gridcell_iy = cellindex3D_coords.icy;
		ptcl.gridcell_iz = cellindex3D_coords.icz;
		ptcl.gridcell_1Dindex = cellindex_1D;

		/*If Lennard-Jones potential is used, then get the LJ grid cell index and add it to
		 * the LJ cavity list of the selected particle type */
		if(particle_types.getLJSwitch(ptcltype_index)==true){

			ptcl.lj_gridcell_1Dindex = cavity_grid.getLJCellindexByCavityCellIndex(cellindex_1D);
			particle_types.addLJCellCavityIndices(ptcltype_index,ptcl.lj_gridcell_1Dindex,cellindex_1D);

		}

		/* Check if there are any other particles that are within the HS cutoff distance.
		 * If yes, then do not proceed to reject this insertion attempt.
		 */
		/*int ptcl_present = particle_state.searchParticleWithinHSCutoff(ptcl,
				particle_types, particlepair_types, cavity_grid);*/

		/*If there is no particle present within  cut-off distance,
		 *  then attempt to insert the particle.
		 */
		/*if (ptcl_present == 0) {

			insert_attempted = true;*/

		/*------- Add particle to particle state vector -------*/

		// Increment the particle type count by one
		particle_types.incrementParticleTypeCountByOne(ptcltype_index);
		int numtype = particle_types.getNum(ptcltype_index);

		particle_state.addParticle(ptcl, ptcltype_index, numtype,
				cellindex_1D);


		/*------- Calculate the new energy and the difference between the new and old energies of the system -------*/
		/*	if(ptcl.label=="Water"){
				std::cout << "[gcmc_calc1.cpp(insertParticle_1)]: attempting to insert water" << std::endl;

			}*/
		double delE = calcParticlePairEnergy(1, parameters.use_spm,parameters.use_mwater,
				parameters.solvent_ion_attraction,
				parameters.solvent_dielectric, particle_state,
				particlepair_types, particle_types, mwater,cavity_grid, box,
				particle_state.getParticleNum() - 1, ptcltype_index + 1);

		if (parameters.solute_model == all_atom && ptcl.charge !=0) {
			delE += calcSoluteParticleEnergy(1, parameters.use_spm,particle_state,particle_types, solute, box,
					particle_state.getParticleNum() - 1);
		}
		double energy_new = energy_old + delE;

		/*------- Calculate the probability of finding a cavity -------*/
		double cavity_prob = particle_types.getNumCavities(ptcltype_index)/ double(num_cells);
		int n = numtype - 1; // number of particles of this type before insertion
		double mu_old = particle_types.getMuOld(ptcltype_index);
		double temp = parameters.temperature;

		/*------- Calculate the acceptance probability for inserting the new particle -------*/
		double probinsert = insertProb(cavity_prob, n, delE, mu_old, temp,
				box.getBoxVol());

		/*------- Apply the acceptance criterion for inserting the new particle ------- */
		// generate a random number
		//	MTRand_closed rand_ins;
		//double raninsert = rand_ins();

		//		if (raninsert <= probinsert) {
		if (rand_ins <= probinsert) {
			energy_old = energy_new; // create a copy of the new energy of the system

			/*Update cavity list of all particle types*/
			particle_state.updateCavityList(ptcl, particle_types,
					particlepair_types, cavity_grid, insertion);

			/*Increment the number of times an insertion attempt was accepted by 1 */
			++num_accept_insert;
			/*	if(ptcl.label=="Water"){
						std::cout << "[gcmc_calc1.cpp(insertParticle_1)]: insertion accepted, delE ="
								<< delE << ", probinsert = " << probinsert << std::endl;

					}*/


			/* Update excess chemical potential if the PID method is used*/
			if(parameters.use_pid == true && parameters.simulation_type == GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION){
				particle_types.updateChemPotentialPID(ptcltype_index,box.getBoxVol(),parameters.temperature);
				//	particle_types.updateChemPotentialPIDwSlothCorr(ptcltype_index,box.getBoxVol(),
				// parameters.temperature,parameters.solvent_dielectric,box.getBoxXLen());

				//	particle_types.updateChemPotentialPIDtargetNum(ptcltype_index,box.getBoxVol(),parameters.temperature);
			}/*else if(parameters.use_pid == true && parameters.simulation_type == GCMC_WITH_SOLUTE_ALLATOM_MODEL){
				particle_types.targetChemPotentialPID(ptcltype_index,box.getBoxVol(),parameters.temperature);
			}*/

			/*------- Get statistics for calculating average quantites and RDFs -------*/
			if (calc_statistics == true) {
				/*	Update the count of sampled states */
				++state_count;

				/* Count the number of each particle type in the state and update the count
				 * from all states sampled so far */
				for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
					double num_avg = particle_types.getParticleType(ind1).num_avg;
					num_avg += particle_types.getParticleType(ind1).num;
					particle_types.updateNumAvg(ind1, num_avg);

				}
				/* Calculate RDF of each particle pair type */
				if (parameters.calculate_particle_pair_rdf == true) {

					if(state_count==1){
						rdf.calcParticlePairBinCounts(particle_types, particle_state);
					}else{
						int ptclindex = particle_state.getParticleNum()-1;
						rdf.calcParticlePairBinCountsPerStateChange(insertion,particle_types,
								particlepair_types, particle_state,ptcl,ptclindex, ptcltype_index);
					}
					rdf.incrementParticlePairBinAvgCounts();

				}

				/* Calculate RDFs of each particle type from the solute */
				if (parameters.calculate_solute_particle_rdf == true) {
					if(state_count == 1){
						rdf.calcSoluteParticleBinCounts(particle_types, particle_state);
					}else{
						int ptclindex = particle_state.getParticleNum()-1;
						rdf.calcSoluteParticleBinCountsPerStateChange(insertion,particle_types,particle_state,ptcl, ptclindex,ptcltype_index);
					}

					rdf.incrementSoluteParticleBinAvgCounts();
				}
			}

		} else { /*------- Reject the insertion attempt ------- */
			particle_state.updateParticleIndicesByCount(ptcltype_index,
					numtype - 1, -1);
			particle_state.updateParticleIndicesByCellindex(ptcltype_index,
					cellindex_1D, -1);

			/*Reduce the number of particles of the selected particle type and the total number of particles  by 1*/
			particle_state.decrementParticleCountByOne();
			particle_types.decrementParticleTypeCountByOne(ptcltype_index);

			/*Remove cavity cellindex from LJ cell*/
			if(particle_types.getLJSwitch(ptcltype_index)==true){
				particle_types.removeCavityCellIndexFromLJCell(ptcltype_index,ptcl.lj_gridcell_1Dindex,cellindex_1D);
				//std::cout << "[gcmc_calc1.cpp(insert)]: rejecting inserted ion, " << ptcl.label
				//	<< ", LJ cell index, " << ptcl.lj_gridcell_1Dindex << ", cavity cell index, " << ptcl.gridcell_1Dindex << std::endl;

			}

		}

		/*} // end of if statement (ptcl_present == 0)
		else {
			insert_attempted = false;
		}*/

	}

	//} // for loop ptclindex

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Deletes a randomly selected particle type.
 *
 * @param[in]		parameters(CInputParameters_t)		Input parameters
 * @param[in]		box (CBox_t)				Simulation box parameters
 * @param[in,out]	particle_state(CParticleState_t)	 Position and type data of all the particles in the system
 * @param[in,out]	particle_types(CParticleType_t)	Parameters of all particle types
 * @param[in]		particlepair_types(CParticlePairType_t)	Parameters of all particle pair types
 * @param[in]		cavity_grid(CavityGrid_t)		Cavity grid data
 * @param[in]		mwater(CMWater)			Parameters of MWater model
 * @param[in]		solute(CSoluteModel_t)		Parameters and data of the solute model
 * @param[in,out]	energy_old(double)		Current energy of the system (kcal/mol)
 * @param[in,out]	num_accept_delete(int)	Number of times a deletion attempt was accepted
 * @param[in]		rand_countindex(double)		A random number [0,1) for randomly selecting  a particle of the selected particle type for deletion
 * @param[in]		rand_del(double)		A random number [0,1] to test the acceptance criterion for accepting a DELETION attempt
 * @param[in]		rand_ptcltype(double)	A random number (0,1] for selecting the type of particle for deletion
 * @param[in]		calc_statistics(bool)	True or False option (boolean) to determine when to start computing the statistics
 * @param[in,out]	state_count(int)		Number of states used for averaging the statistical quantities
 * @param[in,out]	rdf(CRdf_t)				Radial distribution functions and parameters
 *
 * @return void
 */
void deleteParticle_1(const CInputParameters_t &parameters, const CBox_t &box,
		CParticleState_t &particle_state, CParticleType_t &particle_types,
		CParticlePairType_t &particlepair_types, CCavityGrid_t &cavity_grid,  CMWater &mwater,
		const CSoluteModel_t &solute, double &energy_old,
		int &num_accept_delete, double rand_countindex, double rand_del,
		double rand_ptcltype,bool calc_statistics,int &state_count, CRdf_t &rdf) {
	/*
	 MTRand rand_countindex;
	 MTRand_closed rand_del;*/
	/* Get the number of particle types*/
	int nparticle_types = particle_types.getNumParticleTypes();
	/*Get the total number of cavity grid cells */
	int num_cells = cavity_grid.getCavityGridNumCells();

	// Randomly select the type of particle to be deleted
	//MTRand rand_ptcltype;  // [0,1)

	/*Randomly select the type of particle to be deleted. The value will always be
	 * less than the number of particle types (nparticle_types)*/
	int ptcltype_index = int(rand_ptcltype * nparticle_types); // ptcltype_index will always be less than nparticle_types

	//  attempt to delete a particle of each type
	//for (int ptcltype_index = 0; ptcltype_index < nparticle_types; ptcltype_index++){

	/*Get the number of particles of the selected type */
	int numtype = particle_types.getNum(ptcltype_index);

	/*------- Attempt to delete the particle if there is at least one particle of the selected type -------*/
	if (numtype > 1) {
		//	MTRand rand_countindex;

		/*Randomly select the particle to be deleted*/
		int ptclcount_index = int(rand_countindex * numtype);
		int ptclindex = particle_state.getParticleIndicesByCount(ptcltype_index,
				ptclcount_index);

		/*std::cout << "[gcmc_calc1.cpp(deleteParticle_1)]: selected particle at index " << ptclindex <<
		 ", for deletion." << std::endl;*/

		struct_particle ptcl = particle_state.getParticle(ptclindex); // the particle selected for deletion

		/*------- Calculate the new energy and the difference between the new and old energies of the system -------*/
		/*if(ptcl.label=="Water"){
							std::cout << "[gcmc_calc1.cpp(deleteParticle_1)]: attempting to delete water" << std::endl;

						}*/

		double delE = calcParticlePairEnergy(1, parameters.use_spm,parameters.use_mwater,
				parameters.solvent_ion_attraction,
				parameters.solvent_dielectric, particle_state,
				particlepair_types, particle_types, mwater,cavity_grid, box, ptclindex,
				ptcltype_index + 1);

		if (parameters.solute_model == all_atom && ptcl.charge !=0) {
			delE += calcSoluteParticleEnergy(1,parameters.use_spm, particle_state,particle_types, solute, box,
					ptclindex);
		}

		double energy_new = energy_old - delE;

		/* Update the cavity list of each particle type */
		particle_state.updateCavityList(ptcl, particle_types,
				particlepair_types, cavity_grid, deletion);

		/* Calculate the probability of finding a cavity for the selected particle type*/
		double cavity_prob = particle_types.getNumCavities(ptcltype_index)
																	/ double(num_cells);
		int n = particle_types.getNum(ptcltype_index); // number of particles of this type before deletion
		/*Get the current chemical potential of the particle type */
		double mu_old = particle_types.getMuOld(ptcltype_index);
		/*Get the temperature */
		double temp = parameters.temperature;

		/*------- Calculate the acceptance probability for deleting the new particle -------*/

		double probdelete = deleteProb(cavity_prob, n, -delE, mu_old, temp,
				box.getBoxVol());

		/*------- Apply the acceptance criterion for deleting the particle ------*/
		// generate a random number
		//	MTRand_closed rand_del;

		//double randelete = rand_del();

		//	std::cout << "[gcmc_calc1.cpp (deleteParticle_1) random number for deletion (rand_del), "<< randelete << std::endl;

		//		if (randelete <= probdelete) {
		if (rand_del <= probdelete) {
			energy_old = energy_new; // create a copy of the new energy of the system

			/*	if(ptcl.label=="Water"){
						std::cout << "[gcmc_calc1.cpp(deleteParticle_1)]: deleted water, delE = "
								<< delE << ", probdelete = " << probdelete <<std::endl;

					}*/


			/*------- Get statistics for calculating average quantites and RDFs -------*/
			if (calc_statistics == true) {

				/*	Update the count of sampled states */
				++state_count;

				/* Count the number of each particle type in the state and update the count
				 * from all states sampled so far */
				for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
					double num_avg = particle_types.getParticleType(ind1).num_avg;
					if(ind1 == ptcltype_index){
						int tmp = particle_types.getParticleType(ind1).num - 1;
						num_avg += tmp;
					}else{
						num_avg += particle_types.getParticleType(ind1).num;
					}
					particle_types.updateNumAvg(ind1, num_avg);

				}
				/* Calculate RDF of each particle pair type */
				if (parameters.calculate_particle_pair_rdf == true && state_count > 1) {


					rdf.calcParticlePairBinCountsPerStateChange(deletion,particle_types,
							particlepair_types, particle_state,ptcl,ptclindex, ptcltype_index);

					rdf.incrementParticlePairBinAvgCounts();

				}

				/* Calculate RDFs of each particle type from the solute */
				if (parameters.calculate_solute_particle_rdf == true && state_count > 1) {
					rdf.calcSoluteParticleBinCountsPerStateChange(deletion,particle_types,particle_state,ptcl, ptclindex,ptcltype_index);
					rdf.incrementSoluteParticleBinAvgCounts();
				}

			}

			/*Remove cavity cellindex from LJ cell*/
			if(particle_types.getLJSwitch(ptcltype_index)==true){
				particle_types.removeCavityCellIndexFromLJCell(ptcltype_index,ptcl.lj_gridcell_1Dindex,ptcl.gridcell_1Dindex);
				//	std::cout << "[gcmc_calc1.cpp(delete)]: deletion accepted for, " << ptcl.label
				//					<< ", LJ cell index, " << ptcl.lj_gridcell_1Dindex << ", cavity cell index, " << ptcl.gridcell_1Dindex << std::endl;

			}


			/*------- Delete particle data struct from particle state vector -------*/

			particle_state.deleteParticle(ptcl.gridcell_1Dindex,
					ptclcount_index, ptcltype_index, ptclindex, particle_types);

			/*Decrement the number of particles of the selected particle type  by 1 */
			particle_types.decrementParticleTypeCountByOne(ptcltype_index);

			/*Increment the number of times a deletion attempt was accepted by 1 */
			++num_accept_delete;

			/* Update excess chemical potential if the PID method is used*/
			if(parameters.use_pid == true && parameters.simulation_type == GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION){
				particle_types.updateChemPotentialPID(ptcltype_index,box.getBoxVol(),parameters.temperature);
				//particle_types.updateChemPotentialPIDwSlothCorr(ptcltype_index,box.getBoxVol(),
				//parameters.temperature,parameters.solvent_dielectric,box.getBoxXLen());

				//particle_types.updateChemPotentialPIDtargetNum(ptcltype_index,box.getBoxVol(),parameters.temperature);
			}/*else if(parameters.use_pid == true && parameters.simulation_type == GCMC_WITH_SOLUTE_ALLATOM_MODEL){
				particle_types.targetChemPotentialPID(ptcltype_index,box.getBoxVol(),parameters.temperature);
			}*/
			if(state_count==1){
				if(parameters.calculate_particle_pair_rdf==true){
					rdf.calcParticlePairBinCounts(particle_types, particle_state);
					rdf.incrementParticlePairBinAvgCounts();
				}
				if (parameters.calculate_solute_particle_rdf == true) {
					rdf.calcSoluteParticleBinCounts(particle_types, particle_state);
					rdf.incrementSoluteParticleBinAvgCounts();
				}
			}


		} else { /* Restore the cavity list*/
			// No changes required for particle indices vectors, and the number of particles
			particle_state.updateCavityList(ptcl, particle_types,
					particlepair_types, cavity_grid, insertion);

		}

	}
	//} // for loop ptclindex

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Displace a randomly selected particle type.
 *
 * @param[in]		parameters(CInputParameters_t)		Input parameters
 * @param[in]		box(CBox_t)			Simulation box parameters
 * @param[in,out]	particle_state(CParticleState_t)	 Position and type data of all the particles in the system
 * @param[in,out]	particle_types(CParticleType_t)	Parameters of all particle types
 * @param[in]		particlepair_types(CParticlePairType_t)	Parameters of all particle pair types
 * @param[in]		cavity_grid(CavityGrid_t)		Cavity grid data
 * @param[in]		mwater(CMWater)			Parameters of MWater model
 * @param[in]		solute(CSoluteModel_t)			Parameters and data of the solute model
 * @param[in,out]	energy_old(double)		Current energy of the system (kcal/mol)
 * @param[in,out]	num_accept_disp(int)	Number of times a deletion attempt was accepted
 * @param[in]		rand_cavity(double)		Random number [0,1) for randomly selecting  an unoccupied grid cell within a cavity grid segment
 * @param[in]		randpos_x(double)		Random number (0,1) for randomly selecting the x position in a grid cell (currently a fixed value of 0.5 is used)
 * @param[in]		randpos_y(double)		Random number (0,1) for randomly selecting the y position in a grid cell (currently a fixed value of 0.5 is used)
 * @param[in]		randpos_z(double)		Random number (0,1) for randomly selecting the z position in a grid cell (currently a fixed value of 0.5 is used)
 * @param[in]		rand_disp(double)		Random number [0,1] to test the acceptance criterion for accepting a displacement attempt
 * @param[in]		rand_countindex(double)	Random number [0,1) for randomly selecting  a particle of the selected particle type for displacement
 * @param[in]		rand_ptcltype(double)	Random number (0,1] for selecting the type of particle for displacement
 * @param[in]		calc_statistics(bool)	True or False option to determine when to start computing the statistics
 * @param[in,out]	state_count(int)		Number of states used for averaging the statistical quantities
 * @param[in,out]	rdf(CRdf_t)				Radial distribution functions and parameters
 *
 * @return void
 */
void displaceParticle_1(const CInputParameters_t &parameters, const CBox_t &box,
		CParticleState_t &particle_state, CParticleType_t &particle_types,
		CParticlePairType_t &particlepair_types, CCavityGrid_t &cavity_grid, CMWater &mwater,
		const CSoluteModel_t &solute, double &energy_old, int &num_accept_disp,
		double rand_cavity, double randpos_x, double randpos_y,
		double randpos_z, double rand_disp, double rand_countindex,
		double rand_ptcltype,bool calc_statistics,int &state_count,CRdf &rdf) {

	/*MTRand rand_countindex;
	 MTRand rand_cavity;
	 MTRand_open randpos;
	 MTRand_closed rand_disp;
	 */


	/*Get the number of particle types */
	int nparticle_types = particle_types.getNumParticleTypes();
	/*Get the total number of cavity grid cells */
	int num_cells = cavity_grid.getCavityGridNumCells();

	// Randomly select the type of particle to be deleted
	//	MTRand rand_ptcltype;  // [0,1)

	/*Randomly select the type of particle to be displaced. The value will always be
	 * less than the number of particle types (nparticle_types)*/

	int ptcltype_index = int(rand_ptcltype * nparticle_types);

	/*Get the number of particles of the selected type */
	int numtype = particle_types.getNum(ptcltype_index);

	/* Get the number of cavity grid cells available for inserting the center of the
	 * displaced particle of the randomly selected particle type ( ptcltype_index)*/
	int num_cavities = particle_types.getNumCavities(ptcltype_index);

	/*------- Attempt to displace the particle if there is at least one cavity grid cell available
	 * and if there is at least one particle of the selected type ------ */

	if (num_cavities >= 1 && numtype >= 1) {

		/*Randomly select the index of the particle to be deleted*/
		//MTRand rand_countindex;
		int ptclcount_index = int(rand_countindex * numtype);
		int ptclindex = particle_state.getParticleIndicesByCount(ptcltype_index,
				ptclcount_index);

		/*------- Randomly select a cavity from the cavity list of the selected particle -------*/

		/*Declare a random number generator to randomly select a cavity grid segment */

		MTRand rand_cavity_segment;
		/* Initialize a boolean variable to false, which will later be set to True, if
		 * a cavity grid cell is found for inserting the displaced particle*/
		bool found_cavityindex= false;
		/* Get the total number of cavity grid segments */
		int num_cavity_segments = cavity_grid.getNumCavitySegments();
		/*Get the total number of cavity grid cells in the selected segment */
		int num_cells_in_segment =cavity_grid.getNumCellsInSegment();

		int cellindex_1D;
		/*------- Find a cavity grid cell for inserting the displaced particle -------*/
		while(found_cavityindex == false ){
			// select a segment
			int cavity_segmentindex = int(rand_cavity_segment()*num_cavity_segments);
			// get the number of cavities available for the current particle type in the segment.
			int num_cavities_in_segment = particle_types.getNumCavitiesInSegment(ptcltype_index,
					cavity_segmentindex);
			// Select a cavity out of the list of all cavities in the segment for the current particle type.
			int seg_cavitycount_index = int(rand_cavity*num_cavities_in_segment);

			// Get the first cell index of the segment
			int seg_first_cellindex = cavity_grid.getFirstCellIndexOfSegment(cavity_segmentindex);

			// Get the last cell index of the segment
			int seg_last_cellindex = seg_first_cellindex + num_cells_in_segment-1;

			if(num_cavities_in_segment > 0){
				int it1 = seg_first_cellindex;

				int count_index = 0;
				while((it1 <= seg_last_cellindex) && (count_index != seg_cavitycount_index+1) ){
					cellindex_1D = particle_state.getCellIndexByCavityCountIndex(
							ptcltype_index, it1);
					if(cellindex_1D != -1){
						++count_index;
					}
					++it1;
				}
				found_cavityindex = true;
			}else{
				found_cavityindex = false;
			}
		}
		/*-------Set the new grid cell indices and the position for the displaced particle in
		 * the selected grid cell cavity.-------*/

		/*Get 3D grid cell indices */
		struct_cavitygrid_cellindex_coords cellindex3D_coords =
				cavity_grid.getCellIndex3D(cellindex_1D);
		//	MTRand_open randpos;

		/*Create a copy of the selected particle */
		struct_particle ptclcopy = particle_state.getParticle(ptclindex);

		/*Create a copy of the displaced particle, and change its
		 * position and cavity grid cell location*/

		struct_particle ptcl = ptclcopy;

		ptcl.x = cavity_grid.getX(cellindex3D_coords.icx) + randpos_x * cavity_grid.getHx();
		ptcl.y = cavity_grid.getY(cellindex3D_coords.icy)+ randpos_y * cavity_grid.getHy();
		ptcl.z = cavity_grid.getZ(cellindex3D_coords.icz)+ randpos_z * cavity_grid.getHz();

		ptcl.gridcell_ix = cellindex3D_coords.icx;
		ptcl.gridcell_iy = cellindex3D_coords.icy;
		ptcl.gridcell_iz = cellindex3D_coords.icz;
		ptcl.gridcell_1Dindex = cellindex_1D;

		/*		if(ptcl.label=="Water"){
					std::cout << "[gcmc_calc1.cpp(displaceParticle_1)]: attempting to displace water" << std::endl;

				}*/


		/* Check if there are any other particles that are within the HS cutoff distance.
		 * If yes, then do not proceed to reject this insertion attempt.
		 */
		/*int ptcl_present = particle_state.searchParticleWithinHSCutoff(ptcl,
				particle_types, particlepair_types, cavity_grid);*/

		/*If there is no particle present within  cut-off distance,
		 *  then attempt to displace the particle.
		 */
		//if (ptcl_present == 0) {



		/*------- Calculate the energy change due to removing the particle from its
		 * old position -------*/

		//std::cout << "[gcmc_calc1.cpp(displace)]: particle to be displaced " << ptcl.label << std::endl;
		double delE = calcParticlePairEnergy(1, parameters.use_spm,parameters.use_mwater,
				parameters.solvent_ion_attraction,
				parameters.solvent_dielectric, particle_state,
				particlepair_types, particle_types, mwater,cavity_grid, box, ptclindex,
				ptcltype_index + 1);

		if (parameters.solute_model == all_atom && ptcl.charge !=0) {
			delE += calcSoluteParticleEnergy(1,parameters.use_spm, particle_state,particle_types, solute, box,
					ptclindex);
		}

		/*Set value of current particle's cavity grid cell index to -1*/
		particle_state.updateParticleIndicesByCellindex(ptcltype_index,
				ptclcopy.gridcell_1Dindex, -1);

		/* Remove the cavity grid cell index of the particle's new position
		 * from the cavity list associated with LJ grid cell*/
		if(particle_types.getLJSwitch(ptcltype_index)==true){
			particle_types.removeCavityCellIndexFromLJCell(ptcltype_index,ptclcopy.lj_gridcell_1Dindex,ptclcopy.gridcell_1Dindex);

			//	std::cout << "[gcmc_calc1.cpp(disp)]: displacing ion, " << ptcl.label
			//									<< ", FROM LJ cell index, " << ptclcopy.lj_gridcell_1Dindex << ", cavity cell index, " << ptclcopy.gridcell_1Dindex << std::endl;

			/* Add the cavity grid cell index of the particle's old position
			 * to the cavity list associated with LJ grid cell*/
			ptcl.lj_gridcell_1Dindex = cavity_grid.getLJCellindexByCavityCellIndex(ptcl.gridcell_1Dindex);
			particle_types.addLJCellCavityIndices(ptcltype_index,ptcl.lj_gridcell_1Dindex,ptcl.gridcell_1Dindex);

			//std::cout << "[gcmc_calc1.cpp(disp)]: displacing ion, " << ptcl.label
			//											<< ", TO LJ cell index, " << ptcl.lj_gridcell_1Dindex << ", cavity cell index, " << ptcl.gridcell_1Dindex << std::endl;

		}

		/* Replace the old position of the particle's new position, and
		 * update its cavity cell index grid location */

		particle_state.copyParticlePosition(ptclindex, ptcl);

		particle_state.updateParticleIndicesByCellindex(ptcltype_index,
				cellindex_1D, ptclindex);

		/*------- Calculate the change in energy due to inserting the particle at its new position -------*/
		//	std::cout << "[gcmc_calc1.cpp(displace)]: calculate energy due to displaced location of particle " << ptcl.label << std::endl;

		double delE2 = calcParticlePairEnergy(1, parameters.use_spm,parameters.use_mwater,
				parameters.solvent_ion_attraction,
				parameters.solvent_dielectric, particle_state,
				particlepair_types, particle_types, mwater,cavity_grid, box, ptclindex,
				ptcltype_index + 1);

		if (parameters.solute_model == all_atom && ptcl.charge !=0) {
			delE2 += calcSoluteParticleEnergy(1,parameters.use_spm, particle_state,particle_types, solute, box,
					ptclindex);
		}
		delE = -delE + delE2;

		/*Get the temperature */

		double temp = parameters.temperature;

		/*------- Calculate the acceptance probability for displacing the particle -------*/

		double probdisp = dispProb(delE, temp);

		double energy_new = energy_old + delE;
		/*------- Apply the acceptance criterion for displacing the particle ------*/
		// generate a random number
		//MTRand_closed rand_disp;

		//double randisp = rand_disp();
		//		if (randisp <= probdisp) {
		if (rand_disp <= probdisp) {
			energy_old = energy_new;

			/* Update cavity list of all particle types due to removing the
			 * from its old position */

			particle_state.updateCavityList(ptclcopy, particle_types,
					particlepair_types, cavity_grid, deletion);

			/* De-associate the particle's old cavity grid cell index from its
			 *   index in the particles state vector */
			particle_state.updateParticleIndicesByCellindex(ptcltype_index,
					ptclcopy.gridcell_1Dindex, -1);

			/*Update the cavity list of all particle types due to inserting the particle
			 * at its new position  */
			particle_state.updateCavityList(
					particle_state.getParticle(ptclindex), particle_types,
					particlepair_types, cavity_grid, insertion);

			/* Associate the particle's new cavity grid cell index to its
			 *   index in the particles state vector */
			particle_state.updateParticleIndicesByCellindex(ptcltype_index,
					particle_state.getParticle(ptclindex).gridcell_1Dindex,
					ptclindex);

			/*Increment the number of times a displacement attempt is accepted */
			++num_accept_disp;
			/*	if(ptcl.label=="Water"){
						std::cout << "[gcmc_calc1.cpp(displaceParticle_1)]: displaced water, delE = "
								<< delE << ", probdisp = "<< probdisp << std::endl;

					}*/

			/*------- Get statistics for calculating average quantites and RDFs -------*/

			if (calc_statistics == true) {
				/* Increment the number of sampled states */
				++state_count;

				/* Count the number of each particle type in the state and update the count
				 * from all states sampled so far */
				for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
					double num_avg = particle_types.getParticleType(ind1).num_avg;

					num_avg += particle_types.getParticleType(ind1).num;
					particle_types.updateNumAvg(ind1, num_avg);

				}
				/* Calculate RDF of each particle pair type */
				if (parameters.calculate_particle_pair_rdf == true) {

					if(state_count==1){
						rdf.calcParticlePairBinCounts(particle_types, particle_state);

					}else{
						rdf.calcParticlePairBinCountsPerStateChange(deletion,particle_types,
								particlepair_types, particle_state,ptclcopy,ptclindex, ptcltype_index);

						rdf.calcParticlePairBinCountsPerStateChange(insertion,particle_types,
								particlepair_types, particle_state,particle_state.getParticle(ptclindex),ptclindex, ptcltype_index);
					}
					rdf.incrementParticlePairBinAvgCounts();

				}

				/* Calculate RDF of each particle type around the solute*/
				if (parameters.calculate_solute_particle_rdf == true) {

					if(state_count==1){

						rdf.calcSoluteParticleBinCounts(particle_types, particle_state);
					}else{
						rdf.calcSoluteParticleBinCountsPerStateChange(deletion,particle_types,particle_state,ptclcopy, ptclindex,ptcltype_index);

						rdf.calcSoluteParticleBinCountsPerStateChange(insertion,particle_types,particle_state,particle_state.getParticle(ptclindex), ptclindex,ptcltype_index);

					}
					rdf.incrementSoluteParticleBinAvgCounts();

				}

			}

		} else {

			/* Restore the cavity list of the LJ grid cell*/

			if(particle_types.getLJSwitch(ptcltype_index)==true){
				particle_types.removeCavityCellIndexFromLJCell(ptcltype_index,ptcl.lj_gridcell_1Dindex,ptcl.gridcell_1Dindex);
				//std::cout << "[gcmc_calc1.cpp(disp)]: rejecting displacement ion, " << ptcl.label
				//		<< ", FROM LJ cell index, " << ptcl.lj_gridcell_1Dindex << ", cavity cell index, " << ptcl.gridcell_1Dindex << std::endl;

				particle_types.addLJCellCavityIndices(ptcltype_index,ptclcopy.lj_gridcell_1Dindex,ptclcopy.gridcell_1Dindex);
				//	std::cout << "[gcmc_calc1.cpp(disp)]: rejecting displacement ion, " << ptclcopy.label
				//<< ", moved TO ORIGINAL LJ cell index, " << ptclcopy.lj_gridcell_1Dindex << ", cavity cell index, " << ptclcopy.gridcell_1Dindex << std::endl;

			}

			/* Revert the position of the particle to its old position */
			//std::cout << "[gcmc_calc1.cpp(dispParticle1)]: displacement not accepted and cavity update not required. " << std::endl;
			particle_state.revertToParticleToPreviousPosition(ptclcopy,
					ptclindex, ptcltype_index, cellindex_1D);

		}

		//} // end of IF statement (ptcl_present == 0)
	}
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
