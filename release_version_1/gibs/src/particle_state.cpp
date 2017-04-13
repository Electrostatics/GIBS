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
 * particle_state.cpp
 *
 *  Created on: May 11, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file particle_state.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief CParticleState class constructor and member function definitions
 */


#include "particle_state.hpp"
/**
 * @brief	Constructor for the CParticleState class
 *
 * @param[in] nparticle_types(int) Number of particle types
 * @param[in] cavity_grid(CCavityGrid)	Hard sphere cavity grid data
 * @param[in] particlepair_types(CParticlePairType_t)	Parameters of all particle pair typea
 * @param[in]	max_npart(int)	Maximum number of particles allowed in the simulation box (specified by user)
 *
 * @return void
 */
CParticleState::CParticleState(int nparticle_types,const CCavityGrid_t &cavity_grid,
		const CParticlePairType_t &particlepair_types,int max_npart){

	this->nparticle_types = nparticle_types;
	this->nparticle_pairs = particlepair_types.getNumPairs();
	particles.resize(max_npart);
	npart = 0;

	int num_cells = cavity_grid.getNcx()*cavity_grid.getNcy()*cavity_grid.getNcz();
	this->num_cells = num_cells;

	particle_indices_by_count = new int *[nparticle_types];
	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
		particle_indices_by_count[ind1] = new int[max_npart];
	}

	particle_indices_by_cellindex = new int*[nparticle_types];
	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
		particle_indices_by_cellindex[ind1] = new int[num_cells];
	}
	particle_neighbor_count_hs = new int*[nparticle_types];

	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {

		particle_neighbor_count_hs[ind1] = new int[num_cells];
	}


	for (int i=0;i<nparticle_types;i++){

		for (int j=0;j<num_cells;j++){

			particle_indices_by_cellindex[i][j]=-1;  // initially no particles are present in any of the grid cells
			// assuming a grid cell contain only the center of one particle.
			particle_neighbor_count_hs[i][j] = 0;
		}

	}

	cellindex_by_cavity_countindex = new int *[nparticle_types];


	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
		cellindex_by_cavity_countindex[ind1] = new int[num_cells];
	}

	for (int i=0;i<nparticle_types;i++){
		for (int j=0;j<num_cells;j++){
			cellindex_by_cavity_countindex[i][j]=j;
		}
	}

	// allocate memory to hold cavity count indices for each particle pair type

	this->cavitycount = new struct_cavitycount[nparticle_pairs];

	for (int i=0;i<nparticle_pairs;i++){
		std::string pair_label = particlepair_types.getLabelFromIndex(i);
		this->cavitycount[i].indices.resize(particlepair_types.getHardSphereCellIndexOffsetMatrix(pair_label).nindices);
	}
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Initializes the positions of the particles (ions only) in the simulation box without
 * 			the grid insertion method (solvent is modeled as a dielectric continuum)
 *
 * @param[in]	box(CBox_t)		Simulation box parameters
 * @param[in]	min_ptcl_spacing(double)	Minimum spacing between the particles in Angstrom units
 * @param[in]	niontype(int)	Number of ion types in the simulation box
 *
 * @return void
 */
void CParticleState::initializeParticlePositions_WithoutGridBasedInsertion(const CBox_t &box,double min_ptcl_spacing,
		int niontype){

	int nx = int(box.getBoxXLen()/min_ptcl_spacing);
	int ny = int(box.getBoxYLen() / min_ptcl_spacing);
	int nz = int(box.getBoxZLen() / min_ptcl_spacing);

	double xmax = box.getBoxXMax();
	double	ymax = box.getBoxYMax();
	double	zmax =box.getBoxZMax();

	double xmin = box.getBoxXMin();
	double ymin = box.getBoxYMin();
	double zmin = box.getBoxZMin();

	double *x;
	double *y;
	double *z;
	x = new double[nx];
	y = new double[ny];
	z = new double[nz];



	for (int ind = 0; ind < nx; ind++) {
		x[ind] = xmin + double(ind + 1) * min_ptcl_spacing;

	}
	for (int ind = 0; ind < ny; ind++) {
		y[ind] = ymin + double(ind + 1) * min_ptcl_spacing;

	}
	for (int ind = 0; ind < nz; ind++) {
		z[ind] = zmin + double(ind + 1) * min_ptcl_spacing;
	}

	struct_particle ptcli;

	npart = 0;
	for (int index1 = 0; index1 < nx; index1++) {
		for (int index2 = 0; index2 < ny; index2++) {
			for (int index3 = 0; index3 < nz; index3++) {


				if ((x[index1] < xmax) && (y[index2] < ymax)
						&& (z[index3] < zmax)) {
					++npart;
					ptcli.x = x[index1];
					ptcli.y = y[index2];
					ptcli.z = z[index3];

					particles[npart-1] = ptcli;
				}

			}
		}

	}

	delete[] x;
	delete[] y;
	delete[] z;

	x = NULL;
	y = NULL;
	z = NULL;

}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
/**
 * @brief Assigns ion type parameters to particles placed initially in the simulation box
 *
 * @param[in]	iontypes(CIonType_t)	Parameters of all ion types
 * @param[in]	solute_charge(double)	Total charge on the solute in e units.
 * @param[in]	niontype(int)		Number of ion types
 * @param[in]	particle_types(CParticleType_t) 	Parameters of all particle types
 * @param[in]	flag(int)			Number to indicate the different ways to assign ion types
 *
 * @return void
 */
void CParticleState::initializeParticleTypeWithIonParameters(CIonType_t &iontypes,
		double solute_charge, int niontype,CParticleType_t &particle_types, int flag) {

	double charge_sum;

	MTRand rand1;
	int ind, ind1, ind2, ind3;

	struct_iontype iontype;
	switch (flag) {

	case 0: // assign ion types without trying to neutralize the system
		for (ind = 0; ind < niontype; ind++) {

			iontypes.updateIonTypeNum(ind,0);  // initialize number of ions of each type to zero

		}

		for (ind = 0; ind < this->npart; ind++) {
			ind1 = int(rand1() * niontype);

			iontype = iontypes.getIonType(ind1);

			++iontype.num;

			particle_types.assignIonTypeParametersToParticle(particles[ind],iontype);
			iontypes.updateIonTypeNum(ind1,iontype.num);
		}

		break;
	case 1: // assign ion types such that the system is charge neutral
		charge_sum = solute_charge;
		for (ind = 0; ind < niontype; ind++) {
			iontypes.updateIonTypeNum(ind,0);  // initialize number of ions of each type to zero

		}

		for (ind = 0; ind < this->npart; ind++) {
			ind1 = int(rand1() * niontype);
			iontype = iontypes.getIonType(ind1);
			++iontype.num;
			particle_types.assignIonTypeParametersToParticle(particles[ind],iontype);
			iontypes.updateIonTypeNum(ind1,iontype.num);

			charge_sum = charge_sum + particles[ind].charge;
		}

		std::cout
		<< "[particle_state.cpp (initializeParticleTypeWithIonParameters)]: Total sum of charges in the system (in e units) = "
		<< charge_sum << std::endl;

		ind3 = -1;
		ind = -1;
		// remove particles until total charge goes to zero.
		while (charge_sum!=0 && ind < this->npart){
			ind = ind3 + 1;
			if (charge_sum < 0) {

				if (particles[ind].charge < 0) {

					double chg1 = charge_sum - particles[ind].charge;
					if(chg1 <=0){
						charge_sum -= particles[ind].charge;
						int ion_index = particles[ind].ptype - 1;
						int num = iontypes.getIonType(ion_index).num;
						--num;
						iontypes.updateIonTypeNum(ion_index,num);

						for (int i=ind;i<this->npart-1;i++){
							particles[i]= particles[i+1];
						}
						--npart;

						ind3 = ind-1;
					}
				}else{
					ind3 = ind;
				}
			} else if (charge_sum > 0) {
				if (particles[ind].charge > 0) {

					double chg1 = charge_sum - particles[ind].charge;
					if(chg1 >= 0){
						charge_sum = charge_sum - particles[ind].charge;
						int ion_index = particles[ind].ptype - 1;
						int num = iontypes.getIonType(ion_index).num;
						--num;
						iontypes.updateIonTypeNum(ion_index,num);

						for (int i=ind;i<npart-1;i++){
							particles[i]= particles[i+1];
						}
						--npart;

						ind3 = ind-1;
					} else{
						ind3 = ind;
					}
				}
			}
		}
		std::cout
		<< "[particle_state.cpp (initializeParticleTypeWithIonParameters)]: charge_sum after attempting to neutralize the system (in e units) = "
		<< charge_sum << std::endl;

		std::cout
		<< "[particle_state.cpp (initializeParticleTypeWithIonParameters)]: number of particles after attempting to neutralize the system (in e units) = "
		<< npart << std::endl;

		break;
	default:
		std::cout << "[particle_state.cpp (initializeParticleTypeWithIonParameters)]: "
		<< "flag value is not set correctly." << std::endl;
		exit(0);
		break;

	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Assigns the particle type (ions and/or solvent) parameters to particles
 * 			placed initially in the simulation box
 *
 * @param[in]	particle_types(CParticleType_t)	Parameters of all particle types
 * @param[in]	solute_charge(double)	Total charge on solute in e units
 * @param[in]	flag(int)
 *
 * @return void
 */
void CParticleState::initializeParticleType(CParticleType_t &particle_types,
		double solute_charge, int flag) {

	int nparticle_type = particle_types.getNumParticleTypes();

	double charge_sum = 0.0;

	MTRand rand1;
	int ind, ind1, ind2, ind3;

	switch (flag) {

	case 0: // assign particle types without trying to neutralize the system
		for (ind = 0; ind < nparticle_type; ind++) {

			particle_types.updateNumTypes(ind,0); // initialize number of particles of each type to zero

		}

		for (ind = 0; ind < npart; ind++) {
			ind1 = int(rand1() * nparticle_type);
			int num = particle_types.getParticleType(ind1).num;
			++num;
			particle_types.updateNumTypes(ind1,num);
			particle_types.assignParticleTypeParametersToParticle(particles[ind],ind1);
		}

		break;
	case 1:	// assign particle types such that the system is charge neutral
		charge_sum = solute_charge;
		for (ind = 0; ind < nparticle_type; ind++) {

			particle_types.updateNumTypes(ind,0); // initialize number of particles of each type to zero

		}

		for (ind = 0; ind < npart; ind++) {
			ind1 = int(rand1() * nparticle_type);

			int num = particle_types.getParticleType(ind1).num;
			++num;
			particle_types.updateNumTypes(ind1,num);

			particle_types.assignParticleTypeParametersToParticle(particles[ind],ind1);
			charge_sum += particles[ind].charge;
		}



		std::cout
		<< "[particle_state.cpp (initializeParticleType)]: Total sum of charges in the system (in e units) = "
		<< charge_sum << std::endl;

		ind3 = -1;
		ind = -1;
		// remove particles until total charge goes to zero.
		while (charge_sum!=0 && ind <npart){
			ind = ind3 + 1;
			if (charge_sum < 0) {

				if (particles[ind].charge < 0) {

					double chg1 = charge_sum - particles[ind].charge;
					if(chg1 <=0){
						charge_sum = charge_sum - particles[ind].charge;
						int ion_index = particles[ind].ptype - 1;

						int num = particle_types.getParticleType(ion_index).num;
						--num;
						particle_types.updateNumTypes(ion_index,num);

						for (int i=ind;i<npart-1;i++){
							particles[i]= particles[i+1];
						}
						--npart ;

						ind3 = ind-1;
					}
				}else{
					ind3 = ind;
				}
			} else if (charge_sum > 0) {
				if (particles[ind].charge > 0) {

					double chg1 = charge_sum - particles[ind].charge;
					if(chg1 >= 0){
						charge_sum -= particles[ind].charge;
						int ion_index = particles[ind].ptype - 1;
						int num = particle_types.getParticleType(ion_index).num;
						--num;
						particle_types.updateNumTypes(ion_index,num);

						for (int i=ind;i<npart-1;i++){
							particles[i]= particles[i+1];
						}
						--npart;


						ind3 = ind-1;
					} else{
						ind3 = ind;
					}
				}
			}
		}
		std::cout
		<< "[particle_state.cpp (initializeParticleType)]: charge_sum after attempting to neutralize the system (in e units) = "
		<< charge_sum << std::endl;

		std::cout
		<< "[particle_state.cpp (initializeParticleType)]: number of particles after attempting to neutralize the system (in e units) = "
		<< npart << std::endl;

		break;
	default:
		std::cout << "[particle_state.cpp (initializeParticleType)]: "
		<< "flag value is not set correctly." << std::endl;
		exit(0);
		break;

	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Initializes the particle state of the system with grid-based insertion algorithm
 *
 * @param[in]	particle_types(CParticleType_t)	Parameters of all particle types
 * @param[in]	cavity_grid(CCavityGrid_t)	Cavity grid data
 * @param[in]	particlepair_types(CParticlePairType_t)	Parameters of all particle pair type
 * @param[in]	target_num(int[])	Target number of all particle types
 *
 * @return void
 *
 */
void CParticleState::initializeParticleState_WithGridBasedInsertion(CParticleType_t &particle_types,
		const CCavityGrid_t &cavity_grid,CParticlePairType_t &particlepair_types,int target_num[]){

	// insert particles on particle grid
	MTRand randtype;
	MTRand_open randpos;
	MTRand randcavity;
	struct_cavitygrid_cellindex_coords cellindex3D_coords;

	int targ_total = 0;

	int nparticle_type = particle_types.getNumParticleTypes();
	// randomly select the type of particle
	for (int i=0;i<nparticle_type;i++){
		targ_total += target_num[i];
	}

	//int failed_attempts = 0;

	/*int counter = 0;

	while (counter < targ_total){
		int i = int(randtype()*nparticle_type);*/
	int targ_num = 0;
	for (int i=0;i<nparticle_type;i++){
		if (particle_types.getLabel(i) != "Water"){
			targ_num = target_num[i];
		}else if(particle_types.getLabel(i) == "Water"){
			targ_num = 50;
		}
		int counter = 0;

		//if(particle_types[i].num < target_num[i] && particle_types[i].num_cavities!=0 ){



		while(counter<targ_num){
			int num_cavities = particle_types.getNumCavities(i);
			if(num_cavities > 0 ){

				// randomly select a cavity site for placing the particle

				/*std::cout << "[particle_state.cpp(initializeParticleState_WithGridBasedInsertion)]: "
						"current number of cavities for particle " << particle_types[i].label << " = " << particle_types[i].num_cavities
						<< std::endl;
				 */
				int ind1 = int(randcavity()*num_cavities);
				//int cell_index = cellindex_by_cavity_countindex[i][ind1];

				// Cavity count index need not correspond to the index of the cavity in the cell index array.
				int cell_index;
				int it1=0;

				int count_num = 0;

				/*	std::cout << "[particle_state.cpp] number of particles added: " << counter << ", type = " << i <<
						", number of cavities = " << num_cavities << std::endl;*/


				while ((it1 < num_cells) && (count_num != ind1+1)){
					cell_index = cellindex_by_cavity_countindex[i][it1];
					if(cell_index != -1){
						++count_num;
					}
					++it1;
				}

				// get 3D grid cell indices
				cellindex3D_coords = cavity_grid.getCellIndex3D(cell_index);

				struct_particle ptcl;

				particle_types.assignParticleTypeParametersToParticle(ptcl,i);

				/*randomly choose a position for the particle in the selected grid cell cavity.
				*ptcl.x = cavity_grid.getX(cellindex3D_coords.icx)+ randpos()*cavity_grid.getHx() ;
				*	ptcl.y = cavity_grid.getY(cellindex3D_coords.icy)+ randpos()*cavity_grid.getHy();
				*	ptcl.z = cavity_grid.getZ(cellindex3D_coords.icz)+ randpos()*cavity_grid.getHz();
				*/

				ptcl.x = cavity_grid.getX(cellindex3D_coords.icx)+ 0.5*cavity_grid.getHx() ;
				ptcl.y = cavity_grid.getY(cellindex3D_coords.icy)+ 0.5*cavity_grid.getHy();
				ptcl.z = cavity_grid.getZ(cellindex3D_coords.icz)+ 0.5*cavity_grid.getHz();

				/*
			ptcl.x = cavity_grid.getX(cellindex3D_coords.icx)+ 0.5*cavity_grid.getHx()+(randpos()-0.5)*cavity_grid.getHx() ;
					ptcl.y = cavity_grid.getY(cellindex3D_coords.icy)+ 0.5*cavity_grid.getHy() + (randpos()-0.5)*cavity_grid.getHy();
					ptcl.z = cavity_grid.getZ(cellindex3D_coords.icz)+ 0.5*cavity_grid.getHz() + (randpos()-0.5)*cavity_grid.getHz();
				 */


				/*ptcl.x = cavity_grid->x[cellindex3D_coords.icx];
										ptcl.y = cavity_grid->y[cellindex3D_coords.icy];
										ptcl.z = cavity_grid->z[cellindex3D_coords.icz];*/
				ptcl.gridcell_ix = cellindex3D_coords.icx;
				ptcl.gridcell_iy = cellindex3D_coords.icy;
				ptcl.gridcell_iz = cellindex3D_coords.icz;
				ptcl.gridcell_1Dindex = cell_index;

				if(particle_types.getLJSwitch(i)==true){
					ptcl.lj_gridcell_1Dindex = cavity_grid.getLJCellindexByCavityCellIndex(cell_index);
					particle_types.addLJCellCavityIndices(i,ptcl.lj_gridcell_1Dindex,cell_index);
				}

				/*			std::cout << "[particle_state.cpp(initializeParticleState_WithGridBasedInsertion)] " << std::endl;
						std::cout <<"inserting particle of type: " << ptcl.label << std::endl;
						std::cout <<"particle positions (x,y,z): " << ptcl.x << ", " << ptcl.y << ", " << ptcl.z << std::endl;
						std::cout <<"particle grid cell 3D coordinates: " << cavity_grid.getX(cellindex3D_coords.icx)
								<< ", " << cavity_grid.getY(cellindex3D_coords.icy) << ", " <<
								cavity_grid.getZ(cellindex3D_coords.icz) << std::endl;

						std::cout <<" 1D cell index: " << ptcl.gridcell_1Dindex << std::endl;*/

				/* Check if there are any other particles that are within the HS cutoff distance.
				 * If yes, then do not proceed to reject this insertion attempt.
				 */
				/*int ptcl_present = searchParticleWithinHSCutoff(ptcl,
								particle_types, particlepair_types, cavity_grid);*/
				//	if(ptcl_present == 0){
				// add particle to particle state vector
				particle_types.incrementParticleTypeCountByOne(i);
				int num = particle_types.getNum(i);

				++npart;

				//std::cout <<" number of particles added: " << npart << std::endl;
				//std::cout << " number of current particle type: " << particle_types.getLabel(i) << " = " << num << std::endl;

				particles[npart-1] = ptcl;
				particle_indices_by_count[i][num -1] = npart - 1;
				particle_indices_by_cellindex[i][cell_index] = npart - 1;

				/*std::cout << "[intiial_particle_state.cpp(initializeParticleState_WithGridBasedInsertion)]:particle type added,"
					<< ptcl.label << ",  grid cell index = " << ptcl.gridcell_1Dindex << ",  x = " << ptcl.x <<
					", y = " << ptcl.y << ", z = " << ptcl.z << ", " <<
					", particle index = " << particle_indices_by_count[i][particle_types.getNum(i) -1] << std::endl;*/

				// update cavity list for all particles
				// this routine will also update the particle neighbor count


				this->updateCavityList(ptcl,particle_types, particlepair_types,cavity_grid,insertion);

				++counter;
				/*}else{
						if(failed_attempts == 0){
							std::cout << "[particle_state.cpp(initializeParticleState_WithGridBasedInsertion)]:"
									" first failed attempt at particle type: " << ptcl.label << ", counter= " << counter << std::endl;
						}
						++failed_attempts;
					}*/
			}else{
				std::cout << "[particle_state.cpp(initializeParticleState_WithGridBasedInsertion)]: number of cavities = " <<
						num_cavities << " for particle type, " << particle_types.getLabel(i) << std::endl;

				exit(0);
			}

		}
	}

	//std::cout << "[particle_state.cpp:]: Total number of failed attempts: " << failed_attempts << std::endl;
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Read  and initialize particle state of the system from file and set up the cavity lists
 * 			and particle neighbor counts for each particle type
 *
 * @param[in]	particle_types		Class (CParticleType_t) instance containing parameters of each particle type
 * @param[in]	cavity_grid		Class (CCavityGrid_t) instance containing parameters of cavity grid data
 * @param[in]	particlepair_types	Class(CParticlePairType_t) instance containing parameters of each particle pair type
 * @param[in]	file_num		Suffix number (integer) of state file
 * @param[in]	dir_name		Directory name (string) containing the state files
 * @param[in]	target_num		Target number (integer) of each particle type
 *
 * @return void
 */
void CParticleState::readInitializeParticleState_WithGridBasedInsertion(CParticleType_t &particle_types,
		const CCavityGrid_t &cavity_grid,CParticlePairType_t &particlepair_types,int file_num, int num_lines_skip,
		std::string dir_name,int target_num[]){

	int nparticle_type = particle_types.getNumParticleTypes();

	for (int i=0;i<nparticle_type;i++){
			if(particle_types.getLabel(i) == "Water"){
				target_num[i] = 50;
			}
	}


	double energy;
	std::ifstream statefile;


	int i;

	int ind1;
	int col_index;

	int string_count;
	int line_num;

	std::stringstream ss_suffix;
	std::string statefile_name;
	std::string suffix;

	std::string line;
	std::string* sub_strings;
	std::string sub_string;
	std::stringstream ss;

	struct_particle ptcl;
	const char* x = "";
	const char* y = "";
	const char* z = "";
	const char* radius = "";
	const char* chgcoeff = "";
	const char* ptype = "";

	sub_strings = NULL;

	ss_suffix << file_num;
	suffix = ss_suffix.str();

	statefile_name = dir_name;
	statefile_name.append("particle_state_");
	statefile_name.append(suffix);

	line_num = 0;
	i = -1;

	statefile.open(statefile_name.c_str());

	if (statefile.is_open()) {
		while (getline(statefile, line)) {
			ss.clear();
			ss << "";
			ss << line;
			string_count = 0;

			while (ss >> sub_string) {
				++string_count;
			}
			// store the sub strings of a line into an array of strings
			sub_strings = new std::string[string_count];
			col_index = 0;
			ss.clear();
			ss << "";
			ss << line;
			while (ss >> sub_strings[col_index]) {
				++col_index;
			}

			++line_num ;

			if (line_num == 3){
				energy = atof(sub_strings[1].c_str());
			}
			if (line_num > num_lines_skip) {



				if (string_count == 8) {
					x = sub_strings[1].c_str();
					y = sub_strings[2].c_str();
					z = sub_strings[3].c_str();
					radius = sub_strings[4].c_str();
					chgcoeff = sub_strings[5].c_str();
					ptype = sub_strings[6].c_str();

					ptcl.x = atof(x);
					ptcl.y = atof(y);
					ptcl.z = atof(z);
					ptcl.radius = atof(radius);
					ptcl.charge = atof(chgcoeff);
					ptcl.ptype = atoi(ptype);
					ptcl.label = sub_strings[7];

					// find if the cavity is occupied.
					int lix = (int)(ptcl.x-cavity_grid.getXMin())/cavity_grid.getHx();
					int liy  = (int)(ptcl.y-cavity_grid.getYMin())/cavity_grid.getHy();
					int liz  = (int)(ptcl.z-cavity_grid.getZMin())/cavity_grid.getHz();

					int cellindex_1D = cavity_grid.getCellIndex1D(lix,liy,liz);
					int itype = ptcl.ptype-1;
					if(particle_neighbor_count_hs[itype][cellindex_1D]==0 &&
							particle_types.getNum(itype) < target_num[itype]){
						ptcl.gridcell_1Dindex = cellindex_1D;
						ptcl.gridcell_ix = lix;
						ptcl.gridcell_iy = liy;
						ptcl.gridcell_iz = liz;

						if(particle_types.getLJSwitch(itype)==true){
							ptcl.lj_gridcell_1Dindex = cavity_grid.getLJCellindexByCavityCellIndex(cellindex_1D);
							particle_types.addLJCellCavityIndices(itype,ptcl.lj_gridcell_1Dindex,cellindex_1D);
						}
						++i;
						this->particles[i]=ptcl;  // update particle grid cell coordinate and index.
						++npart;/* total number of particles  in the state file */

						particle_types.incrementParticleTypeCountByOne(itype);
						int num = particle_types.getNum(itype);
						particle_indices_by_count[itype][num-1] =  i;
						particle_indices_by_cellindex[itype][cellindex_1D] = i;
						// update cavity list for all particles

						this->updateCavityList(ptcl,particle_types,particlepair_types,cavity_grid,insertion);

					}


				} else {
					std::cout
					<< "[particle_state.cpp (readInitializeParticleState_WithGridBasedInsertion)] Error reading file: "
					<< statefile_name << ", at line number " << line_num
					<< std::endl;
					std::cout << "No. of columns is not equal to 8." << std::endl;
					exit(0);
				}
			}

			delete[] sub_strings;
			sub_strings = NULL;
		}

	} else {
		std::cout << "[particle_state.cpp (readInitializeParticleState_WithGridBasedInsertion)] Error opening file: "
				<< statefile_name << std::endl;
		exit(0);
	}

	statefile.close();

	std::cout << "[particle_state.cpp (readInitializeParticleState_WithGridBasedInsertion)] finished reading file: "
			<< statefile_name << std::endl;
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Count number of each particle type in the system
 *
 * @param[in,out]	particle_types	Class (instance) of particle type
 *
 * @return void
 */
void CParticleState::countParticleTypeFromState(CParticleType_t &particle_types){


	for (int ind1 = 0; ind1 < npart; ind1++) {

		int itype = particles[ind1].ptype-1;
		// Note: particle type index (and labels) in the input state file and parameter file should be the same.
		particle_types.incrementParticleTypeCountByOne(itype);

	}
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Adds particle to particle state vector
 *
 * Maps the particle's vector element index to its 1D cavity grid cell index and its rank
 * in the number of its particle type
 *
 * @param[in]	ptcl(struct_particle)	Particle data structure added to particle state vector (this->particles)
 * @param[in]	ptcltype_index(int)	Array/Vector index of the added particle's type
 * @param[in]	num(int)		Index of the particle's element in particle state vector
 * @param[in]	cellindex_1D(int)	Index of the cavity grid cell containing the particle's center
 *
 * @return void
 */

void CParticleState::addParticle(struct_particle ptcl,int ptcltype_index,int num,int cellindex_1D){
	++npart;
	particles[npart-1]=ptcl;
	particle_indices_by_count[ptcltype_index][num -1] = npart - 1;
	particle_indices_by_cellindex[ptcltype_index][cellindex_1D] = npart - 1;
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Deletes a particle from the particle state vector
 *
 *  Updates the mappings between the particle's vector element index and its
 *  1D cavity grid cell index and its rank in the number of its particle type.
 *
 * @param[in]	cellindex_1D(int)	Array index of the cavity grid cell containing the particle's center
 * @param[in]	ptclcount_index(int)	Particle count index  referring to its rank among other particles of the same type
 * @param[in]	ptclindex(int)		Array/Vector index of the deleted particle in particle state vector
 * @param[in]	particle_types(CParticleType_t)	Parameters of all particle type
 *
 * @return void
 */

void CParticleState::deleteParticle(int cellindex_1D,int ptclcount_index,int ptcltype_index,
		int ptclindex,const CParticleType_t &particle_types){

	/* Decrease the indices of all particles by 1, if their index before
		* deletion of the particle is greater than the index of the deleted particle*/

	for (int ind1 = 0 ; ind1<nparticle_types; ind1++){
		if(ind1==ptcltype_index){
			/*For the deleted particle type, not only the indices of the remaining particles (greater than ptclindex before deletion)
			*have to be decreased by 1, but also their position in the list of particle indices have to be moved up
			* by 1 element.
			*/
			particle_indices_by_cellindex[ptcltype_index][cellindex_1D] = -1;

			for (int ind2 = ptclcount_index;ind2<particle_types.getNum(ptcltype_index)-1;ind2++){

				particle_indices_by_count[ptcltype_index][ind2] = particle_indices_by_count[ptcltype_index][ind2+1];
				if(particle_indices_by_count[ind1][ind2] > ptclindex){
					--particle_indices_by_count[ind1][ind2];
				}
			}
			particle_indices_by_count[ptcltype_index][particle_types.getNum(ptcltype_index)-1] = -1;

		}else {
			for (int ind2 = 0; ind2 < particle_types.getNum(ind1);ind2++){
				if(particle_indices_by_count[ind1][ind2] > ptclindex){
					--particle_indices_by_count[ind1][ind2];
				}
			}
		}
	}
	for (int ind1 = ptclindex;ind1<npart-1;ind1++){
		particles[ind1] = particles[ind1+1];
		int cellindex = particles[ind1].gridcell_1Dindex;
		particle_indices_by_cellindex[particles[ind1].ptype-1][cellindex] = ind1;

	}

	--npart;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief	Copy particle position data to particle state vector
 *
 * @param[in]	i	The vector element index (integer) of the particle in the particle state vector
 * @param[in]	ptcl	The data structure (struct_particle) of the particle whose position is copied
 *
 * @return void
 * s
 */
void CParticleState::copyParticlePosition(int i,struct_particle ptcl){
	particles[i].x = ptcl.x;
	particles[i].y = ptcl.y;
	particles[i].z = ptcl.z;


	particles[i].gridcell_ix = ptcl.gridcell_ix;
	particles[i].gridcell_iy = ptcl.gridcell_iy;
	particles[i].gridcell_iz = ptcl.gridcell_iz;
	particles[i].gridcell_1Dindex = ptcl.gridcell_1Dindex;
	particles[i].lj_gridcell_1Dindex = ptcl.lj_gridcell_1Dindex;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Searches for a particle with hard sphere cutoff radius
 *
 * @param[in]	ptcl (struct_particle)	Particle data structure
 * @param[in]	particle_types(CParticleType_t)	Parameters of all particle types
 * @param[in]	particlepair_types(CParticlePairType_t)	Parameters of all particle pair types
 * @param[in]	cavity_grid(CCavityGrid_t)		Cavity grid data
 *
 * @return 	Number (int) to indicate whether a particle is present (one) or not (zero)
 */
int CParticleState::searchParticleWithinHSCutoff(struct_particle ptcl,CParticleType_t &particle_types,
		const CParticlePairType_t &particlepair_types,const CCavityGrid_t &cavity_grid){

	int ptcl_present = 0;


	double xlen = cavity_grid.getXLen();
	double ylen = cavity_grid.getYLen();
	double zlen = cavity_grid.getZLen();

	int lix = ptcl.gridcell_ix;
	int liy = ptcl.gridcell_iy;
	int liz = ptcl.gridcell_iz;

	int ljx,ljy,ljz;

	int cellindexi_1D = ptcl.gridcell_1Dindex;

	int j = 0;
	while(j<nparticle_types && ptcl_present == 0){

		std::string label12 = ptcl.label;
		label12.append("_");
		label12.append(particle_types.getLabel(j));


		double hs_cutoff = particlepair_types.getHardSphereCutoff(label12);

		struct_cellindex_offset cellindex_offset =
				particlepair_types.getHardSphereCellIndexOffsetMatrix(label12);
		int ncx = cellindex_offset.ncx;
		int ncy = cellindex_offset.ncy;
		int ncz = cellindex_offset.ncz;

		int ind1 = 0;

		while (ind1<2*ncx+1 && ptcl_present == 0){

			ljx = lix + cellindex_offset.delxj[ind1] + cellindex_offset.bcmap_x[lix][ind1];
			int ind2 = 0;
			while (ind2<2*ncy+1 && ptcl_present ==0){
				ljy = liy + cellindex_offset.delyj[ind2] + cellindex_offset.bcmap_y[liy][ind2];
				int ind3 = 0;
				while (ind3 < 2*ncz+1 && ptcl_present ==0){

					ljz = liz + cellindex_offset.delzj[ind3] + cellindex_offset.bcmap_z[liz][ind3];

					int cellindexj_1D = cavity_grid.getCellIndex1D(ljx,ljy,ljz) ;
					int ptclindex = this->getParticleIndicesByCellindex(j,cellindexj_1D);

					if (ptclindex!=-1 && cellindexj_1D!=cellindexi_1D){

						// there is a particle in that cell



						double xj= this->getParticlePositionX(ptclindex);
						double yj = this->getParticlePositionY(ptclindex);
						double zj = this->getParticlePositionZ(ptclindex);

						// Find the nearest image
						if (xj < ptcl.x - xlen * 0.5) {
							xj += xlen;
						} else if (xj > ptcl.x + xlen * 0.5) {
							xj -= xlen;
						}
						if (yj < ptcl.y - ylen * 0.5) {
							yj += ylen;
						} else if (yj > ptcl.y + ylen * 0.5) {
							yj -= ylen;
						}
						if (zj < ptcl.z - zlen * 0.5) {
							zj += zlen;
						} else if (zj > ptcl.z + zlen * 0.5) {
							zj -= zlen;
						}


						double rij = sqrt((ptcl.x-xj)*(ptcl.x-xj)+(ptcl.y-yj)*(ptcl.y-yj) + (ptcl.z-zj)*(ptcl.z-zj));
						if(rij < hs_cutoff){

							ptcl_present = 1;
						}
					}
					++ind3;
				}
				++ind2;
			}

			++ind1;
		}
		++j;
	}


	return ptcl_present;
}
