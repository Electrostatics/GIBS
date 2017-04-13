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
 * particle_type.cpp
 *
 *  Created on: Apr 25, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file particle_type.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief CParticleType class constructor and member function definitions
 */

#include "particle_type.hpp"


/**
 * @brief	Constructor for the CParticleType class
 *
 * Sets up the parameters and cavity grid data for each particle type
 *
 * @param[in] nparticle_types(int) Number of particle types
 * @param[in] cavity_grid(CCavityGrid)	Hard sphere cavity grid data
 *
 */

CParticleType::CParticleType(int nparticle_types,const CCavityGrid &cavity_grid){

	particle_types.resize(nparticle_types);

	int num_cells = cavity_grid.getNcx()*cavity_grid.getNcy()*cavity_grid.getNcz();

	int num_cells_in_segment = cavity_grid.getNumCellsInSegment();
	int num_cavity_segments = cavity_grid.getNumCavitySegments();


	for (int i=0;i<nparticle_types;i++){
		particle_types[i].num_cavities = num_cells; // initial number of cavities available for each particle type
		particle_types[i].num_cavities_in_segment = new int[num_cavity_segments];

		particle_types[i].lj_switch = false;

		for (int j=0; j < num_cavity_segments; j++){
			// initially all the cells are cavities

			particle_types[i].num_cavities_in_segment[j]=num_cells_in_segment;
		}
	}
	this->nparticle_types = nparticle_types;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief	Initialize the LJ grid-based cavity indices
 *
 *
 * @param[in] nparticle_types(int) Number of particle types
 * @param[in] cavity_grid(CCavityGrid)	Hard sphere cavity grid data
 *
 * @return void
 */

void CParticleType::initializeLJCellCavityIndices(int nparticle_types,const CCavityGrid &cavity_grid){

	for (int i=0;i<nparticle_types;i++){
		// If Lennard Jones grid cells are present, then map it to cavity
		// cell indices that are occupied by particles
		int num_lj_cells = cavity_grid.getLJGridNumCells();
		if(num_lj_cells!=0){
			// initialize LJ cell cavity indices
			struct_ljcell_cavityindices ljcell_cavityindices;

			ljcell_cavityindices.num =0;
			int num_cavitycells_perLJcell = cavity_grid.getCavityGridCellsPerLJCell();
			ljcell_cavityindices.cavity_1Dcellindices.resize(num_cavitycells_perLJcell);

			for (int j=0; j< num_lj_cells; j++){

				particle_types[i].ljcell_cavityindices.insert(std::make_pair(j,ljcell_cavityindices));
			}
		}
	}
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Associates the unoccupied cavity grid cell indices to the LJ grid cell for a given
 * 		particle type
 *
 * @param[in]	ptcltype_index(int)	Array/vector index of the particle type
 * @param[in]	lj_cellindex1D(int)	LJ grid cell index number
 * @param[in]	cavity_cellindex1D(int)	Cavity grid cell index
 *
 * @return void
 */
void CParticleType::addLJCellCavityIndices(int ptcltype_index,int lj_cellindex1D,int cavity_cellindex1D){
	struct_ljcell_cavityindices ljcell_cavityindices = particle_types[ptcltype_index].ljcell_cavityindices[lj_cellindex1D];
	ljcell_cavityindices.cavity_1Dcellindices[ljcell_cavityindices.num]= cavity_cellindex1D;
	++ljcell_cavityindices.num;
	particle_types[ptcltype_index].ljcell_cavityindices[lj_cellindex1D] = ljcell_cavityindices;
}
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief  Remove the cavity grid cell index from the unoccupied cavity list of each LJ grid cell
 *
 * @param[in]	ptcltype_index(int)	Array/vector index of the particle type
 * @param[in]	lj_cellindex1D(int)	LJ grid cell index number
 * @param[in]	cavity_cellindex1D(int)	Cavity grid cell index
 *
 * @return void
 */
void CParticleType::removeCavityCellIndexFromLJCell(int ptcltype_index,int lj_cellindex1D,int cavity_cellindex1D){

	struct_ljcell_cavityindices ljcell_cavityindices = particle_types[ptcltype_index].ljcell_cavityindices[lj_cellindex1D];
	bool index_found = false;
	int i = 0;
	while (i < ljcell_cavityindices.num && index_found==false){
		if(ljcell_cavityindices.cavity_1Dcellindices[i]==cavity_cellindex1D){
			index_found = true;
		}
		++i;
	}
	--i;
	if(index_found == true){
		if (i!=ljcell_cavityindices.num-1){
			for (int j  = i; j < ljcell_cavityindices.num-1; j++){
				ljcell_cavityindices.cavity_1Dcellindices[j]=ljcell_cavityindices.cavity_1Dcellindices[j+1];
			}
		}
		--ljcell_cavityindices.num;
		particle_types[ptcltype_index].ljcell_cavityindices[lj_cellindex1D] = ljcell_cavityindices;
	}else{
		std::cerr << "[CParticleType::removeCavityCellIndexFromLJCell]: cavity cell index " << cavity_cellindex1D
				<< " not found in LJ cell " << lj_cellindex1D << std::endl;
	}

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Reads and sets the molar mass values of each particle type
 *
 * @param[in]	fname(char*)	Name of the file  containing the molar mass values (g/mol)
 * @param[in]	nparticle_types (int)	The number of particle types
 *
 * @return void
 */
void CParticleType::readMolarMassAndSetParameters(const char* fname,int nparticle_types){
	std::ifstream paramfile;
	paramfile.open(fname);

	if (paramfile.is_open()) {
		std::cout << "[particle_type.cpp(readMolarMassAndSetParameters)] Reading file " << fname
				<< std::endl;

		int line_num = 0;
		while (!paramfile.eof()) {
			line_num = line_num + 1;
			std::string label1;
			double value;

			if(line_num==1){
				paramfile >> label1 >> label1;
			} else	if(line_num>1){
				paramfile >> label1 >> value;


				for (int i = 0;i<nparticle_types;i++){

					if(label1==particle_types[i].label){
						particle_types[i].molmass = value;
						particle_types[i].targ_dens = concToDensity(particle_types[i].conc,value);
						particle_types[i].epsim1 = - particle_types[i].targ_dens;
					}
				}  // end of for i loop



			} // end of if(line_num) loop

		} // end of while loop

	} // end of IF loop

	else{
		std::cout
		<< "[particle_type.cpp(readMolarMassAndSetParameters)]Error opening file: "
		<< fname << std::endl;
		exit(1);
	}
	paramfile.close();

	std::cout << "[particle_type.cpp(readMolarMassAndSetParameters)]Finished reading file: "
			<< fname << std::endl;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Reads and sets the Lennard Jones parameters of each particle type
 *
 * @param[in]	fname(char*)	Name of the file  containing the LJ parameters
 * @param[in]	nparticle_types (int)	The number of particle types
 *
 * @return void
 */
void CParticleType::readAndSetLennardJonesParameters(const char* fname,int nparticle_types){
	std::ifstream paramfile;
	paramfile.open(fname);

	if (paramfile.is_open()) {
		std::cout << "[particle_type.cpp(readAndSetLennardJonesParameters)] Reading file " << fname
				<< std::endl;

		int line_num = 0;
		while (!paramfile.eof()) {
			line_num = line_num + 1;
			std::string label1;
			double ljsig,ljeps;

			if(line_num==1){
				paramfile >> label1 >> label1 >> label1;
			} else	if(line_num>1){
				paramfile >> label1 >> ljsig >> ljeps;


				for (int i = 0;i<nparticle_types;i++){

					if(label1==particle_types[i].label){
						particle_types[i].ljeps = ljeps;
						particle_types[i].ljsig = ljsig;

					}
				}  // end of for i loop



			} // end of if(line_num) loop

		} // end of while loop

	} // end of IF loop

	else{
		std::cout
		<< "[particle_type.cpp(readAndSetLennardJonesParameters)]Error opening file: "
		<< fname << std::endl;
		exit(1);
	}
	paramfile.close();

	std::cout << "[particle_type.cpp(readAndSetLennardJonesParameters)]Finished reading file: "
			<< fname << std::endl;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Sets the parameters of each particle type
 *
 * @param[in]	ion(CIonType_t)	Ion type parameters
 * @param[in]	solvent(CSolvent_t)	Solvent parameters
 * @param[in]	niontype(int)	Number of ion types
 * @param[in]	nsolvent(int)	Number of solvent types : can only be 0 or 1
 * @param[in]	temp(double)	Temperature in Kelvin
 * @param[in]	box_vol(double)	Simulation box volume in Angstrom cube
 *
 * @return void
 */
void CParticleType::setParameters(const CIonType_t &ion,const CSolvent_t &solvent,
		int niontype,int nsolvent,double temp,double box_vol){


	for(int i=0;i<niontype;i++){
		struct_iontype iont = ion.getIonType(i);

		particle_types[i].charge = iont.charge;// ion.iontypes[i].charge;
		particle_types[i].label = iont.label;
		particle_types[i].radius = iont.radius;
		particle_types[i].ptype = i+1;
		particle_types[i].ljeps = iont.lj_epsilon;
		particle_types[i].ljsig = iont.lj_sigma;
		particle_types[i].conc = iont.conc;
		particle_types[i].num = 0;//iontypes[i].num;
		particle_types[i].mu_ex = iont.mu_ex;
		particle_types[i].mu_old = iont.mu_old;
		particle_types[i].mu_new = iont.mu_new;
		particle_types[i].conc_avg = 0;//iont.cion_avg;
		particle_types[i].num_avg =0;// iont.nion_avg;

		//	particle_types[i].targ_dens = concToDensity(particle_types[i].conc,particle_types);
		particle_types[i].targ_num = concMolarityToNum(particle_types[i].conc,box_vol);
		//	particle_types[i].epsim1 = - particle_types[i].targ_dens;

		particle_types[i].Bi = 1;

		particle_types[i].targ_mu = particle_types[i].mu_new;
		particle_types[i].mu_i = particle_types[i].mu_new;
		/*std::cout <<"[particletype_methods.cpp (getParticleTypeParameters]: label, charge, ptype, conc =" <<
			particle_types[i].label << ", " << particle_types[i].charge <<
			", "<< particle_types[i].ptype << ", " << particle_types[i].conc << std::endl;

		 */
	}
	if (nsolvent==1){
		int nparticletype = niontype + nsolvent;


		//double packing_fraction = solvent.getPackingFraction();


		particle_types[nparticletype-1].charge = solvent.getCharge();
		particle_types[nparticletype-1].label = solvent.getLabel();
		particle_types[nparticletype-1].radius = solvent.getRadius();
		particle_types[nparticletype-1].ptype = nparticletype;

		particle_types[nparticletype-1].mu_ex = solvent.getExcessChemicalPotential();
		particle_types[nparticletype-1].conc = solvent.getConcentration();
		particle_types[nparticletype-1].mu_old = kTLnConc(solvent.getConcentration(), temp) + solvent.getExcessChemicalPotential();
		particle_types[nparticletype-1].mu_new = 	particle_types[nparticletype-1].mu_old;
		particle_types[nparticletype-1].conc_avg = 0;//solvent.getConcentration();
		particle_types[nparticletype-1].num_avg = 0;//concMolarityToNum(solvent.getConcentration(), box_vol);;
		particle_types[nparticletype-1].num = 0;
		/*	std::cout <<"[particletype_methods.cpp (getParticleTypeParameters]: label, charge, ptype, conc =" <<
					particle_types[nparticletype-1].label << ", " << particle_types[nparticletype-1].charge <<
					", "<< particle_types[nparticletype-1].ptype << ", " << particle_types[nparticletype-1].conc << std::endl;
		 */


		//particle_types[nparticletype-1].targ_dens = concToDensity(particle_types[nparticletype-1].conc);
		particle_types[nparticletype-1].targ_num = concMolarityToNum(particle_types[nparticletype-1].conc,box_vol);

		//	particle_types[nparticletype-1].epsim1 = - particle_types[nparticletype-1].targ_dens;

		particle_types[nparticletype -1].Bi = 1;

		particle_types[nparticletype-1].targ_mu = particle_types[nparticletype-1].mu_new;
		particle_types[nparticletype-1].mu_i= particle_types[nparticletype-1].mu_new;
	}

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Reads and sets the parameters of the PID method for calculating the excess chemical potential
 * 		of each particle type.
 *
 * @param[in]	fname	Name of the file (char*) containing the parameters of the PID algorithm
 *
 * @return void
 */
void CParticleType::readAndSetPIDInitValues(const char* fname){

	std::ifstream paramfile;
	paramfile.open(fname);
	std::string label1;
	double bi,epsi,mu_ex;

	if (paramfile.is_open()) {
		std::cout << "[CParticleType::readAndSetPIDInitValues] Reading file " << fname
				<< std::endl;

		int line_num = 0;
		while (!paramfile.eof()) {
			line_num = line_num + 1;


			if(line_num==1){
				paramfile >> label1 >> label1 >> label1 >> label1;
			} else	if(line_num>1){
				paramfile >> label1 >> mu_ex >> bi >> epsi;


				for (int i = 0;i<this->nparticle_types;i++){

					if(label1==particle_types[i].label){
						particle_types[i].mu_ex = mu_ex;
						//particle_types[i].targ_mu = mu_ex + kTLnConc(particle_types[i].conc,);
						particle_types[i].Bi = bi;
						particle_types[i].epsim1 = epsi;
						particle_types[i].mu_ex_m1 = mu_ex;

					}
				}  // end of for i loop



			} // end of if(line_num) loop

		} // end of while loop

	} // end of IF loop

	else
		std::cout
		<< "[CParticleType::readAndSetPIDInitValues]Error opening file: "
		<< fname << std::endl;

	paramfile.close();

	std::cout << "[CParticleType::readAndSetPIDInitValues]Finished reading file: "
			<< fname << std::endl;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Assigns particle type parameters to a particle in the simulation box
 *
 * @param[in,out]	ptcl	Data structure of particle (struct_particle)
 * @param[in]		i		 Index (integer) of the particle type in particle_types vector
 *
 * @return void
 */
void CParticleType::assignParticleTypeParametersToParticle(struct_particle &ptcl, int i) {

	ptcl.ptype =particle_types[i].ptype;
	ptcl.charge = particle_types[i].charge;
	//ptcl.ljeps =particle_types[i].ljeps;
	//ptcl.ljsig = particle_types[i].ljsig;
	ptcl.radius = particle_types[i].radius;
	ptcl.label = particle_types[i].label;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Assign Ion type parameters to a particle in the simulation box
 *
 * @param[in,out]	ptcl(struct_particle)	Data structure of particle
 * @param[in]		iontype(struct_iontype)	Data structure of ion type
 *
 * @return void
 */
void CParticleType::assignIonTypeParametersToParticle(struct_particle &ptcl, struct_iontype iontype) {

	ptcl.ptype = iontype.itype;
	ptcl.charge = iontype.charge;
	//ptcl.ljeps = iontype.lj_epsilon;
	//ptcl.ljsig = iontype.lj_sigma;
	ptcl.radius = iontype.radius;
	ptcl.label = iontype.label;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Writes the chemical potential (excess, total) of each particle type to a file after each iteration
 * 			of a chemical potential calculation
 *
 * @param[in]	iter(int)	The iteration number at which the chemical potential is updated
 *
 * @return voidl
 */
void CParticleType::writeChemicalPotentialToFile(int iter) {

	std::ofstream writefile[nparticle_types];
	int ind;

	std::string filename;
	const char* fname = "";

	for (ind = 0; ind < nparticle_types; ind++) {

		filename = "./outputfiles/chempot_particletype_";

		filename.append(particle_types[ind].label);
		filename.append(".out");
		fname = filename.c_str();

		if (iter == -1) {
			writefile[ind].open(fname);
			writefile[ind] << "iteration#" << '\t' << "avg_conc(M)" << '\t'
					<< "avg_num" << '\t' << "mu_ex(kcal/mol)" << '\t'
					<< "mu_old(kcal/mol)" << '\t' << "mu_new(kcal/mol)" << '\t'
					<< "target_concentration(M)" << std::endl;

		}else {
			writefile[ind].open(fname, std::ios::app);

		}
		/* else if (iter >= 0) {
			writefile[ind] << iter + 1 << '\t' << particle_types[ind].conc_avg << '\t'
					<< particle_types[ind].num_avg << '\t' << particle_types[ind].mu_ex << '\t'
					<< particle_types[ind].mu_old << '\t' << particle_types[ind].mu_new << '\t'
					<< particle_types[ind].conc << std::endl;
		}*/
		writefile[ind] << iter + 1 << '\t' << particle_types[ind].conc_avg << '\t'
				<< particle_types[ind].num_avg << '\t' << particle_types[ind].mu_ex << '\t'
				<< particle_types[ind].mu_old << '\t' << particle_types[ind].mu_new << '\t'
				<< particle_types[ind].conc << std::endl;
	}

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Updates the chemical potential during each iteration of a chemical potential calculation
 * 			based on the adaptive charge-corrected GCMC method
 *
 * @param[in]	iter(int) The iteration number at which the chemical potential is updated
 * @param[in]	parameters(CInputParameters_t)	The input parameters of the simulation
 * @param[in]	box(CBox_t)	Simulation box parameters
 *
 * @return void
 */
void CParticleType::updateChemicalPotential(int iter,const CInputParameters_t &parameters,const CBox_t box){

	double temp = parameters.temperature;
	double solvdiel = parameters.solvent_dielectric;
	int icorr = parameters.icorr;
	double vol = box.getBoxVol();

	double corrTerm = 0.0;
	double boxlen = pow(vol,1.0/3.0);
	double qavg = 0.0;
	if(icorr==1){
		// double KTerm = 12.0*log(2.0+sqrt(3))-2.0*pi;
		//KTerm = 9.520309;
		for(int ind1=0;ind1<nparticle_types;ind1++)
		{
			qavg = qavg + particle_types[ind1].num_avg*particle_types[ind1].charge;
			// charge in e units
		}
	}

	for(int ind1=0;ind1<nparticle_types;ind1++)
	{
		//	double ptcl_num = particle_types[ind1].num_avg;
		//particle_types[ind1].conc_avg=numToConcMolarity(ptcl_num,vol);
		/* MGB iteration process */
		particle_types[ind1].mu_ex = particle_types[ind1].mu_old-kTLnConc(particle_types[ind1].conc_avg,temp);
		if(icorr==1){
			//corrTerm = -avognum*convJtoCal*0.001*ioni.valency*fundchg*fundchg*qavg*KTerm/
			//(32.0*pi*vacpermit*solvdiel*boxlen*pow(10,-10));
			//corrTerm = -393.8258973*ioni.valency*qavg/(solvdiel*boxlen); (until version 3_17_2013)
			corrTerm = -395.1623633*particle_types[ind1].charge*qavg/(solvdiel*boxlen);
		}
		particle_types[ind1].mu_new = kTLnConc(particle_types[ind1].conc,temp)+particle_types[ind1].mu_ex + corrTerm;


		std::cout << "[particle_type.cpp(CParticleType::updateChemicalPotential)] Iteration# = " << iter+1 << std::endl;
		std::cout << "[update_chemical_potential.cpp: chemPotUpdate] average number of particles of type, " << particle_types[ind1].label << " = "
				<< particle_types[ind1].num_avg << ", <mu_ex> = " << particle_types[ind1].mu_ex << " mu_new = " << particle_types[ind1].mu_new << ", mu_old = " << particle_types[ind1].mu_old << std::endl;
		std::cout << "[particle_type.cpp(CParticleType::updateChemicalPotential)] conc_avg = " << particle_types[ind1].conc_avg << std::endl;
		std::cout << "-----" << std::endl;

		/* --- store the new updated value into mu_old --- */
		//        ion1[ind1].mu_old = ion1[ind1].mu_new;

	}

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb                                                                                                                                      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Updates the chemical potential of each particle type during a chemical potential calculation
 * 			using PID method, based on the target mass concentration of the particle type
 *
 * @param[in]	ind1(int)	Index of the particle type in the the particle_types vector
 * @param[in]	vol(double)		Volume of the simulation box in Angstrom cube
 * @param[in]	temp(double)	Temperature  in Kelvin
 *
 * @return void
 */
void CParticleType::updateChemPotentialPID(int ind1,double vol,double temp){


		particle_types[ind1].epsi = numToDensity(particle_types[ind1].num,vol,particle_types[ind1].molmass)
				- particle_types[ind1].targ_dens;

	//particle_types[ind1].epsi = numToConcMolarity(particle_types[ind1].num,vol) - particle_types[ind1].conc;


	particle_types[ind1].Bip1 = particle_types[ind1].Bi - 112.0*(particle_types[ind1].epsi -
			particle_types[ind1].epsim1) -112*particle_types[ind1].epsi /13000;

	particle_types[ind1].Bi = particle_types[ind1].Bip1;
	particle_types[ind1].epsim1 = particle_types[ind1].epsi;

	particle_types[ind1].mu_ex = 0.00198719137 * temp * (particle_types[ind1].Bip1 - log(particle_types[ind1].num));
	particle_types[ind1].mu_i = particle_types[ind1].mu_old;
	particle_types[ind1].mu_old = kTLnConc(particle_types[ind1].conc,temp)+particle_types[ind1].mu_ex;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Updates the chemical potential of each particle type in the PID simulations based
 * 			on target concentration of each particle type, with Sloth Sorensen charge correction
 *
 * @param[in]	ind1(int)	Index of the particle type in the the particle_types vector
 * @param[in]	vol(double)		Volume of the simulation box in Angstrom cube
 * @param[in]	temp(double)	Temperature (double) in Kelvin
 * @param[in]	solvdiel(double)	Solvent dielectic coefficient
 * @param[in]	boxlen(double)		Box length, assumed to be the same along all axes in chemical potential calculations
 *
 * @return void
 */
void CParticleType::updateChemPotentialPIDwSlothCorr(int ind1,double vol,double temp,double solvdiel,double boxlen){

	particle_types[ind1].epsi = numToConcMolarity(particle_types[ind1].num,vol) - particle_types[ind1].conc;


	particle_types[ind1].Bip1 = particle_types[ind1].Bi - 112.0*(particle_types[ind1].epsi -
			particle_types[ind1].epsim1) -112*particle_types[ind1].epsi /13000;

	particle_types[ind1].Bi = particle_types[ind1].Bip1;
	particle_types[ind1].epsim1 = particle_types[ind1].epsi;

	particle_types[ind1].mu_ex = 0.00198719137 * temp * (particle_types[ind1].Bip1 - log(particle_types[ind1].num));

	particle_types[ind1].mu_i = particle_types[ind1].mu_old;
	double qavg;
	for(int ind2=0;ind2<nparticle_types;ind2++)
	{
		qavg = qavg + particle_types[ind2].num*particle_types[ind2].charge;
		// charge in e units
	}

	for (int ind2=0;ind2<nparticle_types;ind2++){
		double corrTerm = -395.1623633*particle_types[ind2].charge*qavg/(solvdiel*boxlen);
		particle_types[ind2].mu_old = kTLnConc(particle_types[ind2].conc,temp)+particle_types[ind2].mu_ex + corrTerm;

	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Updates chemical potential of each particle type in the PID simulations based on
 * 			target number of each particle type
 *
 * @param[in]	ind1(int	Index of the particle type in the the particle_types vector
 * @param[in]	vol(double)		Simulation box volume in Angstrom units
 * @param[in]	temp(double)	Temperature in Kelvin
 *
 * @return void
 */
void CParticleType::updateChemPotentialPIDtargetNum(int ind1,double vol,double temp){

	particle_types[ind1].epsi = particle_types[ind1].num - particle_types[ind1].targ_num;


	particle_types[ind1].Bip1 = particle_types[ind1].Bi - 112.0*(particle_types[ind1].epsi -
			particle_types[ind1].epsim1) -112*particle_types[ind1].epsi /13000;

	particle_types[ind1].Bi = particle_types[ind1].Bip1;
	particle_types[ind1].epsim1 = particle_types[ind1].epsi;

	particle_types[ind1].mu_ex = 0.00198719137 * temp * (particle_types[ind1].Bip1 - log(particle_types[ind1].num));
	particle_types[ind1].mu_i = particle_types[ind1].mu_old;
	particle_types[ind1].mu_old = kTLnConc(particle_types[ind1].conc,temp)+particle_types[ind1].mu_ex;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Updates chemical potential of each particle type in the PID simulations based on
 * 			target chemical potential of each particle type
 *
 * @param[in]	ind1(int)	Index (integer) of the particle type in the the particle_types vector
 * @param[in]	vol(double)		Simulation box volume  in Angstrom units
 * @param[in]	temp(double)	Temperature in Kelvin
 *
 * @return void
 */
void CParticleType::targetChemPotentialPID(int ind1,double vol,double temp){

	particle_types[ind1].epsi = particle_types[ind1].mu_old - particle_types[ind1].targ_mu;
	particle_types[ind1].mu_ex = particle_types[ind1].mu_ex_m1 - 112.0*(particle_types[ind1].epsi -
			particle_types[ind1].epsim1) -112*particle_types[ind1].epsi /13000;

	particle_types[ind1].mu_ex_m1 = particle_types[ind1].mu_ex;

	double conc_avg = numToConcMolarity(particle_types[ind1].num, vol);

	particle_types[ind1].mu_old = particle_types[ind1].mu_ex + kTLnConc(conc_avg,temp);
	particle_types[ind1].mu_i = particle_types[ind1].mu_old;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Replaces the old chemical potential value with new one
 *
 * @return void
 */
void CParticleType::replaceOldChemicalPotentialWithNew(){
	for (int ind2 = 0; ind2 < nparticle_types; ind2++) {
		particle_types[ind2].mu_old = particle_types[ind2].mu_new; /* store the new updated value of mu into mu_old */
	}
}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxs
/**
 * @brief Writes out the parameters and chemical potential values at any given step of a PID simulation
 *
 * @param[in]	istep1(int)	Simulation step number
 *
 * @return void
 */
void CParticleType::writePIDStateParameters(int istep1){

	std::ofstream writefile;

	int ind;
	std::string filename = "./outputfiles/pid_state_params.out";
	const char* fname = filename.c_str();

	for (ind = 0; ind < nparticle_types; ind++) {


		if(istep1>-1){
			writefile.open(fname, std::ios::app);
		}else{
			writefile.open(fname);
			writefile << "PARTICLE_TYPE" << '\t' << "STATE_NUM" << '\t' << "PARTICLE_COUNT"
					<< '\t' << "BI" << '\t' << "EPSI" << '\t' << "MUEX" << '\t' << "TARGET_MASS_DENSITY" << std::endl;
		}

		writefile << istep1 + 1 << '\t' << particle_types[ind].label << '\t' << istep1 + 1 <<
				'\t' << particle_types[ind].num << '\t' << particle_types[ind].Bi << '\t' <<
				particle_types[ind].epsi << '\t' << particle_types[ind].mu_ex << '\t' <<
				particle_types[ind].targ_dens << std::endl;
		writefile.close();

	}
}
