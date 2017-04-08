/*
 * particlepair_type.cpp
 *
 *  Created on: Apr 25, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file particlepair_type.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief CParticlePairType class constructor and member function definitions
 *
 */

#include "particlepair_type.hpp"
#include "input_parameters.hpp"

/**
 * @brief Reads and sets the square well potential parameters for each particle pair type
 *
 * @param[in]	fname(char*)	Name of the file from which the square well potential parameters are
 *
 * @return void
 */

void CParticlePairType::readAndSetSquareWellParameters(const char* fname){
	std::ifstream paramfile;

	paramfile.open(fname);

	if (paramfile.is_open()) {
		std::cout << "[particlepair_type.cpp(CParticlePairType::readAndSetSquareWellParameters)] Reading file " << fname
				<< std::endl;

		int line_num = 0;
		while (!paramfile.eof()) {
			line_num = line_num + 1;
			std::string label1,label2;
			double well_depth;
			double wdth2contact_ratio;  // well with to contact ratio

			if(line_num==1){

				paramfile >> label1 >> label1 >> label1 >> label1;  // first line is header line
			}else if(line_num>1){
				paramfile >> label1 >> label2 >> well_depth >> wdth2contact_ratio;
				label1.append("_");
				label1.append(label2);

				for (int i = 0;i<nparticle_pairs;i++){

					if(label1==particlepair_types[i].label12 || label1==particlepair_types[i].label21){
						particlepair_types[i].sqrwell_depth = well_depth;
						particlepair_types[i].sqrwell_wdth2contact_ratio = wdth2contact_ratio;

						std::cout << "[CParticlePairType::readSquareWellParameters]: well distance:"
								<< (particlepair_types[i].sqrwell_wdth2contact_ratio+1.0)*(
										this->getHardSphereCutoff(label1))  << std::endl;
					}

				}  // end of for i loop


			} // end of if(line_num) loop

		} // end of while loop

	} // end of IF loop

	else
		std::cout
		<< "[particlepair_type.cpp(CParticlePairType::readAndSetSquareWellParameters)]Error opening file: "
		<< fname << std::endl;

	paramfile.close();

	std::cout << "[particlepair_type.cpp(CParticlePairType::readAndSetSquareWellParameters)]Finished reading file: "
			<< fname << std::endl;
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Sets the Lennard Jones potential parameters for each particle pair type
 *
 * @param[in] particle_types(CParticleType_t) Parameters of each particle type
 *
 * @return void
 */
void CParticlePairType::setLennardJonesParameters(const CParticleType_t &particle_types){
	for(int i=0;i<nparticle_pairs;i++){
		if (this->particlepair_types[i].lennard_jones_potential==true){
			int itype1 = particlepair_types[i].ptype1-1;
			int itype2 = particlepair_types[i].ptype2-1;
			particlepair_types[i].ljsig = 0.5*(particle_types.getParticleType(itype1).ljsig + particle_types.getParticleType(itype2).ljsig);
			particlepair_types[i].ljeps = sqrt(particle_types.getParticleType(itype1).ljeps * particle_types.getParticleType(itype2).ljeps);

			particlepair_types[i].Aij = 4.0*particlepair_types[i].ljeps*
					pow(particlepair_types[i].ljsig,12.0);

			particlepair_types[i].Bij = 4.0*particlepair_types[i].ljeps*
					pow(particlepair_types[i].ljsig,6.0);

			std::cout << "[particlepair_type.cpp(setLennardJonesParameters)]: particle pair " << particlepair_types[i].label12 << std::endl;
			std::cout << "LJ eps = " << particlepair_types[i].ljeps << std::endl;
			std::cout << " LJ sig = " << particlepair_types[i].ljsig << std::endl;
			std::cout << " LJ Aij = " << particlepair_types[i].Aij << std::endl;
			std::cout << " LJ Bij = " << particlepair_types[i].Bij << std::endl;
		}

	}
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Reads and sets the types of interactions specified by the user
 * 			for each particle pair type
 *
 *
 * @param[in]	fname(char*)	Name of the file where the user specifies interactions
 *
 * @return void
 */
void CParticlePairType::readAndSetInteractionTypes(const char* fname){

	std::ifstream paramfile;

	paramfile.open(fname);

	if (paramfile.is_open()) {
		std::cout << "[particlepair_type.cpp(readAndSetInteractionTypes)] Reading file " << fname
				<< std::endl;

		int line_num = 0;
		while (!paramfile.eof()) {
			line_num = line_num + 1;
			std::string label1,label2;
			bool hs,coul,lj,sqrp,lookup;

			if(line_num==1){
				paramfile >> label1 >> label1 >> label1 >> label1 >> label1 >> label1 >> label1;  // for sake of reading the headers.
			} else if(line_num>1){
				paramfile >> label1 >> label2 >> hs >> coul >> lj >> sqrp >> lookup;
				label1.append("_");
				label1.append(label2);

				for (int i = 0;i<nparticle_pairs;i++){

					if(label1==particlepair_types[i].label12 || label1==particlepair_types[i].label21){
						particlepair_types[i].hard_sphere_repulsion = hs;
						particlepair_types[i].coulomb_potential = coul;
						particlepair_types[i].lennard_jones_potential = lj;
						particlepair_types[i].square_well_potential = sqrp;
						particlepair_types[i].pmf_lookup_table = lookup;


					}

				}  // end of for i loop


			} // end of if(line_num) loop

		} // end of while loop

	} // end of IF loop

	else
		std::cout
		<< "[particlepair_type.cpp(readAndSetInteractionTypes)]Error opening file: "
		<< fname << std::endl;

	paramfile.close();

	std::cout
	<< "[particlepair_type.cpp(readAndSetInteractionTypes)] Finished reading file: "
	<< fname << std::endl;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Sets the id and labels of each particle pair type
 *
 * @param[in]	nparticle_types(int)	Number of particle types
 * @param[in]	particle_types(CParticleType_t)	Instance of particle type class
 *
 * @return void
 */

void CParticlePairType::setTypeAndLabels(int nparticle_types,const CParticleType_t &particle_types){
	std::string label1, label2, label12, label21;
	for(int i=0;i<nparticle_types;i++){

		label1= particle_types.getParticleType(i).label;
		label2 = label1;
		label12 = label1;
		label12.append("_");
		label12.append(label2);

		label21 = label2;
		label21.append("_");
		label21.append(label1);

		particlepair_types[i].index = i;
		particlepair_types[i].label1 = label1;
		particlepair_types[i].label2 = label2;
		particlepair_types[i].label12 = label12;
		particlepair_types[i].label21 = label21;
		particlepair_types[i].ptype1 = particle_types.getParticleType(i).ptype;
		particlepair_types[i].ptype2 = particle_types.getParticleType(i).ptype;
		particlepair_types[i].radius1 = particle_types.getParticleType(i).radius;
		particlepair_types[i].radius2 = particle_types.getParticleType(i).radius;

		label2index.insert(std::make_pair(label12,i));
		index2label.insert(std::make_pair(i,label12));

	}

	int index = nparticle_types-1;
	for(int i=0;i<nparticle_types;i++){
		for (int j=i+1;j<nparticle_types;j++){
			label1= particle_types.getParticleType(i).label;
			label2 = particle_types.getParticleType(j).label;
			label12 = label1;
			label12.append("_");
			label12.append(label2);



			label21 = label2;
			label21.append("_");
			label21.append(label1);

			++index;
			particlepair_types[index].index = index;
			particlepair_types[index].label1 = label1;
			particlepair_types[index].label2 = label2;
			particlepair_types[index].label12 = label12;
			particlepair_types[index].label21 = label21;
			particlepair_types[index].ptype1 = particle_types.getParticleType(i).ptype;
			particlepair_types[index].ptype2 = particle_types.getParticleType(j).ptype;
			particlepair_types[index].radius1 = particle_types.getParticleType(i).radius;
			particlepair_types[index].radius2 = particle_types.getParticleType(j).radius;

			label2index.insert(std::make_pair(label12,index));
			label2index.insert(std::make_pair(label21,index));  // unique index; multiple labels of the same pair containing different particle types
			index2label.insert(std::make_pair(index,label12));

		}
	}

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Reads and sets the hard sphere potential cutoff distance between the two particle
 * 			types of each particle pair type
 *
 * @param[in]	fname(char*)	Name of the file where the user specifies the cut-off distance
 *
 * @return void
 */
void CParticlePairType::readAndSetHardSphereCutoffs(const char* fname){

	std::ifstream paramfile;

	paramfile.open(fname);

	if (paramfile.is_open()) {
		std::cout << "[particlepair_type.cpp(readAndSetHardSphereCutoffs)] Reading file " << fname
				<< std::endl;

		int line_num = 0;
		while (!paramfile.eof()) {
			line_num = line_num + 1;
			std::string label1,label2, hs_dist;


			if(line_num==1){
				paramfile >> label1 >> label1 >> label1;  // for sake of reading the headers.
			} else if(line_num>1){
				paramfile >> label1 >> label2 >> hs_dist;
				label1.append("_");
				label1.append(label2);

				//	std::cout << label1 << "," << hs << ", " << coul  << ", " << lj << ", " << sqrp
				//		<< std::endl;

				for (int i = 0;i<nparticle_pairs;i++){

					if(label1==particlepair_types[i].label12 || label1==particlepair_types[i].label21){
						this->hs_cutoff.insert(std::make_pair(particlepair_types[i].label12,atof(hs_dist.c_str())));
						this->hs_cutoff.insert(std::make_pair(particlepair_types[i].label21,atof(hs_dist.c_str())));
					}

				}  // end of for i loop


			} // end of if(line_num) loop

		} // end of while loop

	} // end of IF loop

	else
		std::cout
		<< "[particlepair_type.cpp(readAndSetHardSphereCutoffs)]Error opening file: "
		<< fname << std::endl;

	paramfile.close();

	std::cout
	<< "[particlepair_type.cpp(readAndSetHardSphereCutoffs)]Finished reading file: "
	<< fname << std::endl;

}
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Sets the Lennard-Jones potential cutoff distance between the two particle types
 * 			of each particle pair type
 *
 *
 * @param[in]	lj_cutoff(double)	L-J cutoff distance in Angstrom
 *
 * @return void
 */
void CParticlePairType::setLJCutoffs(double lj_cutoff){

	std::string label1,label2;

	for (int i = 0;i<nparticle_pairs;i++){

		this->lj_cutoff.insert(std::make_pair(particlepair_types[i].label12,lj_cutoff));
		this->lj_cutoff.insert(std::make_pair(particlepair_types[i].label21,lj_cutoff));

	}  // end of for i loop
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Reads and sets the Lennard-Jones potential cutoff distance between the two
 * 		particle types of each particle pair
 *
 * @param[in]	fname(char*)	Name  of the file where the user specifies the cutoff distances
 *
 * @return void
 */
void CParticlePairType::readAndSetLJCutoffs(const char* fname){

	std::ifstream paramfile;

	paramfile.open(fname);

	if (paramfile.is_open()) {
		std::cout << "[particlepair_type.cpp(readAndSetLJCutoffs)] Reading file " << fname
				<< std::endl;

		int line_num = 0;
		while (!paramfile.eof()) {
			line_num = line_num + 1;
			std::string label1,label2, lj_dist;


			if(line_num==1){
				paramfile >> label1 >> label1 >> label1;  // for sake of reading the headers.
			} else if(line_num>1){
				paramfile >> label1 >> label2 >> lj_dist;
				label1.append("_");
				label1.append(label2);

				for (int i = 0;i<nparticle_pairs;i++){

					if(label1==particlepair_types[i].label12 || label1==particlepair_types[i].label21){
						this->lj_cutoff.insert(std::make_pair(particlepair_types[i].label12,atof(lj_dist.c_str())));
						this->lj_cutoff.insert(std::make_pair(particlepair_types[i].label21,atof(lj_dist.c_str())));
					}

				}  // end of for i loop


			} // end of if(line_num) loop

		} // end of while loop

	} // end of IF loop

	else
		std::cout
		<< "[particlepair_type.cpp(readAndSetLJCutoffs)] Error opening file: "
		<< fname << std::endl;

	paramfile.close();

	std::cout
	<< "[particlepair_type.cpp(readAndSetLJCutoffs)] Finished reading file: "
	<< fname << std::endl;

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/*void CParticlePairType::setHardSphereCutOff(std::string label12,double hs_dist){
	hs_cutoff.insert(std::make_pair(label12,hs_dist));

	double hs_dist;
	for (int i=0; i < nparticle_pairs; i++){
		double hs_dist = particlepair_types[i].radius1 +particlepair_types[i].radius2;


		hs_cutoff.insert(std::make_pair(particlepair_types[i].label12,hs_dist));
		if(particlepair_types[i].label12!=particlepair_types[i].label21){
			hs_cutoff.insert(std::make_pair(particlepair_types[i].label21,hs_dist));  // mapping multiple label pairs to the same hard sphere cut-off
		}
	}

}*/
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Sets the monoatomic water potential cutoff distance between two solvent (water) hard spheres
 *
 * @param[in]	rc_cutoff(double)	MWater potential cutoff distance in Angstrom
 *
 * @return void
 */
void CParticlePairType::setMWaterPotentialCutoff(const double rc_cutoff){

	mw_cutoff.insert(std::make_pair("Water_Water",rc_cutoff));
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Sets the offset values for determining the cavity grid cell indices
 * that fall within the hard sphere cut off distance of each particle pair type
 *
 * @param[in]	cavity_grid(CCavityGrid_t)	Cavity grid data
 *
 * @return void
 */

void CParticlePairType::setHardSphereCellIndexOffsetMatrix(const CCavityGrid_t &cavity_grid){


	for (int ind1 = 0; ind1 < nparticle_pairs; ind1++){
		struct_cellindex_offset cellindex_offset;

		int Ncx = cavity_grid.getNcx();
		int Ncy = cavity_grid.getNcy();
		int Ncz = cavity_grid.getNcz();


		cellindex_offset.Ncx = Ncx;
		cellindex_offset.Ncy = Ncy;
		cellindex_offset.Ncz = Ncz;
		cellindex_offset.hcx = cavity_grid.getHx();
		cellindex_offset.hcy = cavity_grid.getHy();
		cellindex_offset.hcz = cavity_grid.getHz();


		cellindex_offset.particlepair_index = ind1;
		cellindex_offset.particlepair_label = particlepair_types[ind1].label12;

		double cutoff = getHardSphereCutoff(particlepair_types[ind1].label12);

		int ncx = int (ceil(cutoff/cellindex_offset.hcx))-1;
		int ncy = int (ceil(cutoff/cellindex_offset.hcy))-1;
		int ncz = int (ceil(cutoff/cellindex_offset.hcz))-1;


		std::cout << "[particlepair_type.cpp(setHardSphereCellIndexOffsetMatrix)]: particle pair label,  "
				<< cellindex_offset.particlepair_label << std::endl;
		std::cout << "[particlepair_type.cpp(setHardSphereCellIndexOffsetMatrix)]: ncx =  "
				<< ncx << std::endl;
		std::cout << "[particlepair_type.cpp(setHardSphereCellIndexOffsetMatrix)]: ncy =  "
				<< ncy << std::endl;
		std::cout << "[particlepair_type.cpp(setHardSphereCellIndexOffsetMatrix)]: ncz =  "
				<< ncz << std::endl;

		cellindex_offset.ncx = ncx;
		cellindex_offset.ncy = ncy;
		cellindex_offset.ncz = ncz;
		cellindex_offset.rcutoff = cutoff;
		cellindex_offset.delxj = new int[2*ncx+1];
		cellindex_offset.delyj = new int[2*ncy+1];
		cellindex_offset.delzj = new int[2*ncz+1];


		for(int i=0;i<2*ncx+1;i++){
			cellindex_offset.delxj[i]=-ncx+i;
		}
		for(int i=0;i<2*ncy+1;i++){
			cellindex_offset.delyj[i]=-ncy+i;
		}
		for(int i=0;i<2*ncz+1;i++){
			cellindex_offset.delzj[i]=-ncz+i;
		}

		// map periodic boundary conditions on cell neighbor indices.
		cellindex_offset.bcmap_x = new int*[Ncx];
		for (int i=0;i<Ncx;i++){
			cellindex_offset.bcmap_x[i] = new int[2*ncx+1];
		}
		cellindex_offset.bcmap_y = new int*[Ncy];
		for (int i=0;i<Ncy;i++){
			cellindex_offset.bcmap_y[i] = new int[2*ncy+1];
		}
		cellindex_offset.bcmap_z = new int*[Ncz];
		for (int i=0;i<Ncz;i++){
			cellindex_offset.bcmap_z[i] = new int[2*ncz+1];
		}

		for(int i=0; i < Ncx; i++){
			for (int j=0;j<2*ncx+1;j++){
				cellindex_offset.bcmap_x[i][j] = 0;
				if(i+cellindex_offset.delxj[j] > Ncx-1){
					cellindex_offset.bcmap_x[i][j] = -Ncx;
				}else if(i+cellindex_offset.delxj[j] < 0){
					cellindex_offset.bcmap_x[i][j] = Ncx;
				}

			}

		}


		for(int i=0; i < Ncy; i++){
			for (int j=0;j<2*ncy+1;j++){
				cellindex_offset.bcmap_y[i][j] = 0;
				if(i+cellindex_offset.delyj[j] > Ncy-1){
					cellindex_offset.bcmap_y[i][j] = -Ncy;
				}else if(i+cellindex_offset.delyj[j] < 0){
					cellindex_offset.bcmap_y[i][j] = Ncy;
				}

			}

		}


		for(int i=0; i < Ncz; i++){
			for (int j=0;j<2*ncz+1;j++){
				cellindex_offset.bcmap_z[i][j] = 0;
				if(i+cellindex_offset.delzj[j] > Ncz-1){
					cellindex_offset.bcmap_z[i][j] = -Ncz;
				}else if(i+cellindex_offset.delzj[j] < 0){
					cellindex_offset.bcmap_z[i][j] = Ncz;
				}

			}

		}


		cellindex_offset.nindices = (2*ncx+1)*(2*ncy+1)*(2*ncz+1);

		hs_cellindex_neighbor_offset.insert(std::make_pair(particlepair_types[ind1].label12, cellindex_offset));

		if(particlepair_types[ind1].label12!=particlepair_types[ind1].label21){
			hs_cellindex_neighbor_offset.insert(std::make_pair(particlepair_types[ind1].label21, cellindex_offset));  //adding the same pair with labels arranged in reverse order.
		}


		//particlepair_types[ind1].cavitycount_indices.resize(cellindex_offset.nindices);
	}

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Sets the offset values for determining the cavity grid cell indices
 * that fall within the square well potential cut off distance of each particle pair type
 *
 * @param[in]	cavity_grid(CCavityGrid_t)		Cavity grid data
 *
 * @return void
 */
void CParticlePairType::setSqrWellWidthCellIndexOffsetMatrix(const CCavityGrid_t &cavity_grid){

	for (int ind1 = 0; ind1 < nparticle_pairs; ind1++){
		struct_cellindex_offset cellindex_offset;

		int Ncx = cavity_grid.getNcx();
		int Ncy = cavity_grid.getNcy();
		int Ncz = cavity_grid.getNcz();


		cellindex_offset.Ncx = Ncx;
		cellindex_offset.Ncy = Ncy;
		cellindex_offset.Ncz = Ncz;
		cellindex_offset.hcx = cavity_grid.getHx();
		cellindex_offset.hcy = cavity_grid.getHy();
		cellindex_offset.hcz = cavity_grid.getHz();

		cellindex_offset.particlepair_index = ind1;
		cellindex_offset.particlepair_label = particlepair_types[ind1].label12;

		std::string label12=this->index2label.find(ind1)->second;

		double cutoff = (particlepair_types[ind1].sqrwell_wdth2contact_ratio + 1.0)*(
				this->getHardSphereCutoff(label12));

		int ncx = int (ceil(cutoff/cellindex_offset.hcx))-1;
		int ncy = int (ceil(cutoff/cellindex_offset.hcy))-1;
		int ncz = int (ceil(cutoff/cellindex_offset.hcz))-1;

		std::cout << "[CParticlePairType::setSqrWellWidthCellIndexOffsetMatrix]: particle pair label,  "
				<< cellindex_offset.particlepair_label << std::endl;
		std::cout << "[CParticlePairType::setSqrWellWidthCellIndexOffsetMatrix]: ncx =  "
				<< ncx << std::endl;

		std::cout << "[CParticlePairType::setSqrWellWidthCellIndexOffsetMatrix]: ncy =  "
				<< ncy << std::endl;
		std::cout << "[CParticlePairType::setSqrWellWidthCellIndexOffsetMatrix]: ncz =  "
				<< ncz << std::endl;

		cellindex_offset.ncx = ncx;
		cellindex_offset.ncy = ncy;
		cellindex_offset.ncz = ncz;
		cellindex_offset.rcutoff = cutoff;
		cellindex_offset.delxj = new int[2*ncx+1];
		cellindex_offset.delyj = new int[2*ncy+1];
		cellindex_offset.delzj = new int[2*ncz+1];



		for(int i=0;i<2*ncx+1;i++){
			cellindex_offset.delxj[i]=-ncx+i;
		}
		for(int i=0;i<2*ncy+1;i++){
			cellindex_offset.delyj[i]=-ncy+i;
		}
		for(int i=0;i<2*ncz+1;i++){
			cellindex_offset.delzj[i]=-ncz+i;
		}

		// map periodic boundary conditions on cell neighbor indices.
		cellindex_offset.bcmap_x = new int*[Ncx];
		for (int i=0;i<Ncx;i++){
			cellindex_offset.bcmap_x[i] = new int[2*ncx+1];
		}
		cellindex_offset.bcmap_y = new int*[Ncy];
		for (int i=0;i<Ncy;i++){
			cellindex_offset.bcmap_y[i] = new int[2*ncy+1];
		}
		cellindex_offset.bcmap_z = new int*[Ncz];
		for (int i=0;i<Ncz;i++){
			cellindex_offset.bcmap_z[i] = new int[2*ncz+1];
		}
		for(int i=0; i < Ncx; i++){
			for (int j=0;j<2*ncx+1;j++){
				cellindex_offset.bcmap_x[i][j] = 0;
				if(i+cellindex_offset.delxj[j] > Ncx-1){
					cellindex_offset.bcmap_x[i][j] = -Ncx;
				}else if(i+cellindex_offset.delxj[j] < 0){
					cellindex_offset.bcmap_x[i][j] = Ncx;
				}

			}

		}


		for(int i=0; i < Ncy; i++){
			for (int j=0;j<2*ncy+1;j++){
				cellindex_offset.bcmap_y[i][j] = 0;
				if(i+cellindex_offset.delyj[j] > Ncy-1){
					cellindex_offset.bcmap_y[i][j] = -Ncy;
				}else if(i+cellindex_offset.delyj[j] < 0){
					cellindex_offset.bcmap_y[i][j] = Ncy;
				}

			}

		}


		for(int i=0; i < Ncz; i++){
			for (int j=0;j<2*ncz+1;j++){
				cellindex_offset.bcmap_z[i][j] = 0;
				if(i+cellindex_offset.delzj[j] > Ncz-1){
					cellindex_offset.bcmap_z[i][j] = -Ncz;
				}else if(i+cellindex_offset.delzj[j] < 0){
					cellindex_offset.bcmap_z[i][j] = Ncz;
				}

			}

		}


		cellindex_offset.nindices = (2*ncx+1)*(2*ncy+1)*(2*ncz+1);

		sqrwell_cellindex_neighbor_offset.insert(std::make_pair(particlepair_types[ind1].label12, cellindex_offset));

		if(particlepair_types[ind1].label12!=particlepair_types[ind1].label21){
			sqrwell_cellindex_neighbor_offset.insert(std::make_pair(particlepair_types[ind1].label21, cellindex_offset));//adding the same pair with labels arranged in reverse order.
		}

	}


}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Sets the offset values for determining the cavity grid cell indices
 * that fall within the monoatomic water potential cut off distance of solvent (water) pair type
 *
 * @param[in]	cavity_grid(CCavityGrid_t)		Cavity grid data
 *
 * @return void
 */

void CParticlePairType::setMWaterCellIndexOffsetMatrix(const CCavityGrid_t &cavity_grid){
	int ind1 = this->label2index.find("Water_Water")->second;


	struct_cellindex_offset cellindex_offset;

	int Ncx = cavity_grid.getNcx();
	int Ncy = cavity_grid.getNcy();
	int Ncz = cavity_grid.getNcz();


	cellindex_offset.Ncx = Ncx;
	cellindex_offset.Ncy = Ncy;
	cellindex_offset.Ncz = Ncz;
	cellindex_offset.hcx = cavity_grid.getHx();
	cellindex_offset.hcy = cavity_grid.getHy();
	cellindex_offset.hcz = cavity_grid.getHz();

	cellindex_offset.particlepair_index = ind1;
	cellindex_offset.particlepair_label = particlepair_types[ind1].label12;



	double cutoff1 = mw_cutoff[particlepair_types[ind1].label12];
	double cutoff = ceil(cutoff1); // extra length added

	int ncx = int (ceil(cutoff/cellindex_offset.hcx))-1;
	int ncy = int (ceil(cutoff/cellindex_offset.hcy))-1;
	int ncz = int (ceil(cutoff/cellindex_offset.hcz))-1;

	std::cout << "[CParticlePairType::setMWaterCellIndexOffsetMatrix]: particle pair label,  "
			<< cellindex_offset.particlepair_label << std::endl;
	std::cout << "[CParticlePairType::setMWaterCellIndexOffsetMatrix]: ncx =  "
			<< ncx << std::endl;

	std::cout << "[CParticlePairType::setMWaterCellIndexOffsetMatrix]: ncy =  "
			<< ncy << std::endl;
	std::cout << "[CParticlePairType::setMWaterCellIndexOffsetMatrix]: ncz =  "
			<< ncz << std::endl;

	cellindex_offset.ncx = ncx;
	cellindex_offset.ncy = ncy;
	cellindex_offset.ncz = ncz;
	cellindex_offset.rcutoff = cutoff;
	cellindex_offset.delxj = new int[2*ncx+1];
	cellindex_offset.delyj = new int[2*ncy+1];
	cellindex_offset.delzj = new int[2*ncz+1];



	for(int i=0;i<2*ncx+1;i++){
		cellindex_offset.delxj[i]=-ncx+i;
	}
	for(int i=0;i<2*ncy+1;i++){
		cellindex_offset.delyj[i]=-ncy+i;
	}
	for(int i=0;i<2*ncz+1;i++){
		cellindex_offset.delzj[i]=-ncz+i;
	}

	// map periodic boundary conditions on cell neighbor indices.
	cellindex_offset.bcmap_x = new int*[Ncx];
	for (int i=0;i<Ncx;i++){
		cellindex_offset.bcmap_x[i] = new int[2*ncx+1];
	}
	cellindex_offset.bcmap_y = new int*[Ncy];
	for (int i=0;i<Ncy;i++){
		cellindex_offset.bcmap_y[i] = new int[2*ncy+1];
	}
	cellindex_offset.bcmap_z = new int*[Ncz];
	for (int i=0;i<Ncz;i++){
		cellindex_offset.bcmap_z[i] = new int[2*ncz+1];
	}
	for(int i=0; i < Ncx; i++){
		for (int j=0;j<2*ncx+1;j++){
			cellindex_offset.bcmap_x[i][j] = 0;
			if(i+cellindex_offset.delxj[j] > Ncx-1){
				cellindex_offset.bcmap_x[i][j] = -Ncx;
			}else if(i+cellindex_offset.delxj[j] < 0){
				cellindex_offset.bcmap_x[i][j] = Ncx;
			}

		}

	}


	for(int i=0; i < Ncy; i++){
		for (int j=0;j<2*ncy+1;j++){
			cellindex_offset.bcmap_y[i][j] = 0;
			if(i+cellindex_offset.delyj[j] > Ncy-1){
				cellindex_offset.bcmap_y[i][j] = -Ncy;
			}else if(i+cellindex_offset.delyj[j] < 0){
				cellindex_offset.bcmap_y[i][j] = Ncy;
			}

		}

	}


	for(int i=0; i < Ncz; i++){
		for (int j=0;j<2*ncz+1;j++){
			cellindex_offset.bcmap_z[i][j] = 0;
			if(i+cellindex_offset.delzj[j] > Ncz-1){
				cellindex_offset.bcmap_z[i][j] = -Ncz;
			}else if(i+cellindex_offset.delzj[j] < 0){
				cellindex_offset.bcmap_z[i][j] = Ncz;
			}

		}

	}


	cellindex_offset.nindices = (2*ncx+1)*(2*ncy+1)*(2*ncz+1);

	mwater_cellindex_neighbor_offset.insert(std::make_pair(particlepair_types[ind1].label12, cellindex_offset));

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Reads PMF data for each particle pair type
 *
 * The name of the file containing the PMF data should be in the format
 * lookup_pmf_{particle pair type label}.in (e.g., lookup_pmf_Na_Cl.in,
 * lookup_pmf_Na_Na.in)
 *
 * @param[in]	parameters(CInputParameters_t)	Simulation run parameters
 *
 * @return void
 */
void CParticlePairType::readPMFLookupFile(const CInputParameters_t &parameters){

	double delr_pmf = parameters.lookup_pmf_delr;
	int numr_pmf = parameters.lookup_pmf_numr;


	std::ifstream pmflookup_file;
	std::string line;

	int index;
	std::string sub_strings[2];
	std::string sub_string;

	std::stringstream ss;
	std::vector<double> r_lookup;
	std::vector<double> pmf_lookup;
	int n_lookup = 0;
	int inc_size = 0;
	int nsize = 100;
	r_lookup.resize(nsize);
	pmf_lookup.resize(nsize);

	for (int i=0;i<nparticle_pairs;i++){
		if(particlepair_types[i].pmf_lookup_table==true){


			std::string filename = "./inputfiles/lookup_pmf_";
			filename.append(particlepair_types[i].label12);
			filename.append(".in");
			const char* fname = filename.c_str();

			std::cout << "[particlepair_types.cpp(CParticlePairType::readPMFLookupFile)]: reading file, " <<
					filename << std::endl;


			pmflookup_file.open(fname);



			// count the number of lines to read
			int num_lines = 0;

			if(pmflookup_file.is_open()){
				while(getline(pmflookup_file,line)){
					ss.clear();
					ss << "";
					ss << line;
					int index = 0;

					while (ss >> sub_strings[index]) {
						++index;
					}
					if(num_lines>=r_lookup.size()){
						++inc_size;
						r_lookup.resize((inc_size+1)*nsize);
						pmf_lookup.resize((inc_size+1)*nsize);
					}
					r_lookup[num_lines] = strtod(sub_strings[0].c_str(),NULL);
					pmf_lookup[num_lines] = strtod(sub_strings[1].c_str(),NULL);

					/*std::cout << "[particlepair_types.cpp(CParticlePairType::readPMFLookupFile)]: num_lines = " << num_lines
							<< ", r = " << r_lookup[num_lines] << ", pmf = " << pmf_lookup[num_lines] << std::endl;*/
					++num_lines;

				}

				n_lookup = num_lines;

				std::cout << "[particlepair_types.cpp(CParticlePairType::readPMFLookupFile)]: number of lines read (n_lookup) = "
						<< n_lookup << std::endl;

				std::vector<struct_spline_coeff> spline_coeff;
				spline_coeff.resize(n_lookup);

				create_splineCoeff(pmf_lookup,spline_coeff,n_lookup);

				particlepair_types[i].n_lookup = n_lookup;
				particlepair_types[i].r_lookup.resize(n_lookup);
				particlepair_types[i].pmf_lookup.resize(n_lookup);
				particlepair_types[i].spline_coeff.resize(n_lookup);

				particlepair_types[i].r_lookup = r_lookup;
				particlepair_types[i].pmf_lookup = pmf_lookup;
				particlepair_types[i].spline_coeff = spline_coeff;


				// Get the PMF values for r values in the range set by the parameters, delr and numr
				/*particlepair_types[i].n_lookup = numr_pmf;
				particlepair_types[i].r_lookup.resize(numr_pmf);
				particlepair_types[i].pmf_lookup.resize(numr_pmf);*/

				/*std::cout << "[particlepair_types.cpp(CParticlePairType::readPMFLookupFile)]: interpolating values for user-specified "
						"range and precision of distance values." << std::endl;


				for (int indk = 0; indk < numr_pmf; indk++) {
					double	rij = particlepair_types[i].radius1 + particlepair_types[i].radius2 + indk * delr_pmf;

					std::cout << "[particlepair_types.cpp(CParticlePairType::readPMFLookupFile)]: indk = " << indk <<
							", rij = " << rij << std::endl;

					particlepair_types[i].r_lookup[indk] = rij;

					 interpolate the pmf value from the first lookup table \


					if (rij < r_lookup[0]) {
						// give a very large value, if the radial distance between the two ions is
						// less than the lowest value found in the lookup table.

						particlepair_types[i].pmf_lookup[indk] = 100000.0;

					} else if (rij > r_lookup[n_lookup - 1]) {

						particlepair_types[i].pmf_lookup[indk] = 0;

					} else {

						int	indm = binarySearch(r_lookup, n_lookup,	rij);

						if (indm != n_lookup - 1) {
							double	x = (rij - r_lookup[indm])
																			/ (r_lookup[indm+1]- r_lookup[indm]);

							particlepair_types[i].pmf_lookup[indk] =
									cspline_interp(spline_coeff[indm].cubic,x);

						} else if (indm == n_lookup - 1) {

							particlepair_types[i].pmf_lookup[indk] =
									pmf_lookup[n_lookup- 1];
						}

					}

				}

				// re-calculate the spline coefficients.


				create_splineCoeff(particlepair_types[i].pmf_lookup,
						particlepair_types[i].spline_coeff,particlepair_types[i].n_lookup);*/

			}else {
				std::cerr
				<< "[particlepair_type.cpp(CParticlePairType::readPMFLookupFile)]Error opening file: "
				<< fname << std::endl;
				exit(0);
			}

			// close file
			pmflookup_file.close();
		}
	}


}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Writes PMF lookup table and spline coefficients
 *  of each particle pair type to a file
 *
 * The file name has the format pmf_lookup_{particle pair type label}.out (e.g.,
 * pmf_lookup_Na_Cl.out)
 *
 * @return void
 */
void CParticlePairType::writePMFLookupTableToFile() {


	std::ofstream writefile;
	const char* fname = "";

	for (int i = 0; i < nparticle_pairs; i++) {
		if(particlepair_types[i].pmf_lookup_table == true){
			std::string filename = "./outputfiles/pmf_lookup_";
			filename.append(particlepair_types[i].label12);

			filename.append(".out");
			fname = filename.c_str();

			writefile.open(fname, std::ios::trunc);
			writefile << "r" << '\t' << "PMF(kcal/mol)" << '\t' << "a0" << '\t'
					<< '\t' << "a1" << '\t' << "a2" << '\t' << "a3" << std::endl;

			writefile.close();
			writefile.open(fname, std::ios::app);

			for (int j = 0; j < particlepair_types[i].n_lookup; j++) {
				if (j < particlepair_types[i].n_lookup - 1) {
					writefile << particlepair_types[i].r_lookup[j] << '\t'
							<< particlepair_types[i].pmf_lookup[j] << '\t'
							<< particlepair_types[i].spline_coeff[j].cubic[0] << '\t'
							<< particlepair_types[i].spline_coeff[j].cubic[1] << '\t'
							<< particlepair_types[i].spline_coeff[j].cubic[2] << '\t'
							<< particlepair_types[i].spline_coeff[j].cubic[3] << std::endl;
				} else {
					writefile << particlepair_types[i].r_lookup[j] << '\t'
							<< particlepair_types[i].pmf_lookup[j] << '\t' << "" << '\t' << ""
							<< '\t' << "" << '\t' << "" << std::endl;
				}
			}
			writefile.close();

		}
	}
}
