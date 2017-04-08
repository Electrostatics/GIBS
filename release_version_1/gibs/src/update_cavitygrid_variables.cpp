/*
 * update_cavitygrid_variables.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file update_cavitygrid_variables.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief CParticleState member function definitions for updating the cavity grid
 *  		variables
 */

#include "particle_state.hpp"

/**
 * @brief Updates the cavity grid variables based on the position and number of all the particles
 * 			in the simulation box
 *
 * @param[in]	particle_types(CParticleType_t)		Selected particle type data
 * @param[in]	cavity_grid(CCavityGrid_t)			Cavity grid data
 * @param[in]	particlepair_types(CParticlePairType_t)	Selected particle pair type data
 *
 * @return void
 */
void CParticleState::updateCavityGridVariables_AllParticles(CParticleType_t &particle_types,
		const CCavityGrid_t &cavity_grid,const CParticlePairType_t &particlepair_types){


	int nparticle_type = particle_types.getNumParticleTypes();

	struct_cavitygrid_cellindex_coords cellindex3D_coords;
	int ptcl_num[nparticle_type];
	for (int itype=0;itype<nparticle_type;itype++){
		ptcl_num[itype] = 0;
	}

	for (int ipart=0;ipart<npart; ipart++){
		struct_particle ptcl = particles[ipart];

		// get the cell index coordinates of the particle

		cellindex3D_coords.icx = (int)(ptcl.x-cavity_grid.getXMin())/cavity_grid.getHx();
		cellindex3D_coords.icy = (int)(ptcl.y-cavity_grid.getYMin())/cavity_grid.getHy();
		cellindex3D_coords.icz = (int)(ptcl.z-cavity_grid.getZMin())/cavity_grid.getHz();

		int lix = cellindex3D_coords.icx;
		int liy = cellindex3D_coords.icy;
		int liz = cellindex3D_coords.icz;

		int cellindex_1D = cavity_grid.getCellIndex1D(lix,liy,liz);

		ptcl.gridcell_1Dindex = cellindex_1D;
		ptcl.gridcell_ix = lix;
		ptcl.gridcell_iy = liy;
		ptcl.gridcell_iz = liz;

		int itype = ptcl.ptype-1;
		if(particle_types.getLJSwitch(itype)==true){
			ptcl.lj_gridcell_1Dindex = cavity_grid.getLJCellindexByCavityCellIndex(cellindex_1D);
			particle_types.addLJCellCavityIndices(itype,ptcl.lj_gridcell_1Dindex,cellindex_1D);
		}

		particles[ipart] = ptcl;  // update particle grid cell coordinate and index.

		++ptcl_num[itype];

		particle_indices_by_count[itype][ptcl_num[itype]-1] =  ipart;
		particle_indices_by_cellindex[itype][cellindex_1D] = ipart;
		// update cavity list for all particles
		// this routine will also update the particle neighbor count

		this->updateCavityList(ptcl,particle_types,particlepair_types,cavity_grid,insertion);

	}

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Updates the cavity grid variables based on the position and number of all the atoms
 * 			in the all-atom solute model
 *
 * @param[in]	particle_type(CParticleType_t)		Selected particle type data
 * @param[in]	cavity_grid(CCavityGrid_t)		Cavity grid data
 * @param[in]	solute_model(CSoluteModel_t)		Solute model data
 *
 * @return void
 */
void CParticleState::updateCavityGridVariables_SoluteAllAtom(CParticleType_t &particle_types,
		const CCavityGrid_t &cavity_grid,const CSoluteModel_t &solute_model){

	std::vector<int> cavitycount_indices;
	cavitycount_indices.resize(8000);

	int natoms = solute_model.getSoluteAllAtom_NumAtoms();
	int nparticle_types = particle_types.getNumParticleTypes();

	int num_cells_in_segment = cavity_grid.getNumCellsInSegment();
	int Ncx = cavity_grid.getNcx();
	int Ncy = cavity_grid.getNcy();
	int Ncz = cavity_grid.getNcz();

	double hcx = cavity_grid.getHx();
	double hcy = cavity_grid.getHy();
	double hcz = cavity_grid.getHz();


	double xlen = cavity_grid.getXLen();
	double ylen = cavity_grid.getYLen();
	double zlen = cavity_grid.getZLen();


	for (int i=0;i<natoms;i++){
		double xatom = solute_model.getSoluteAllAtom_AtomCoordX(i);
		double yatom = solute_model.getSoluteAllAtom_AtomCoordY(i);
		double zatom = solute_model.getSoluteAllAtom_AtomCoordZ(i);

		int lix = (int)(xatom-cavity_grid.getXMin())/cavity_grid.getHx();
		int liy = (int) (yatom-cavity_grid.getYMin())/cavity_grid.getHy();
		int liz = (int) (zatom-cavity_grid.getZMin())/cavity_grid.getHz();


		int atom_cellindex_1D = cavity_grid.getCellIndex1D(lix,liy,liz);

		for (int j=0;j<nparticle_types;j++){

			double cutoff = solute_model.getSoluteAllAtom_AtomRadius(i)+particle_types.getRadius(j);

			int ncx = int (ceil(cutoff/hcx))-1;
			int ncy = int (ceil(cutoff/hcy))-1;
			int ncz = int (ceil(cutoff/hcz))-1;

			int ljx,ljy,ljz;

			int ncavity = 0;

			for (int ind1=0;ind1<2*ncx+1;ind1++){


				ljx = lix -ncx +ind1;// cellindex_offset.delxj[ind1];


				if (ljx > Ncx-1){
					ljx = ljx - Ncx;
				}else if(ljx < 0){
					ljx = ljx + Ncx;
				}


				for (int ind2 = 0; ind2<2*ncy+1; ind2++){
					ljy = liy -ncy + ind2; //cellindex_offset.delyj[ind2];

					if(ljy > Ncy-1){
						ljy = ljy - Ncy;
					}else if(ljy < 0){
						ljy = ljy + Ncy;
					}

					for (int ind3 = 0; ind3 < 2*ncz+1; ind3++){

						ljz = liz -ncz + ind3;
						if (ljz > Ncz - 1){
							ljz = ljz - Ncz;
						}else if(ljz < 0){
							ljz = ljz + Ncz;
						}

						int cellindex_1D = cavity_grid.getCellIndex1D(ljx,ljy,ljz) ;

						double	xj = cavity_grid.getXCenter(ljx);
						double	yj = cavity_grid.getYCenter(ljy);
						double	zj = cavity_grid.getZCenter(ljz);

						// Find the nearest image
						if (xj < xatom - xlen * 0.5) {
							xj += xlen;
						} else if (xj > xatom + xlen * 0.5) {
							xj -= xlen;
						}
						if (yj < yatom - ylen * 0.5) {
							yj += ylen;
						} else if (yj > yatom + ylen * 0.5) {
							yj -= ylen;
						}
						if (zj < zatom - zlen * 0.5) {
							zj += zlen;
						} else if (zj > zatom + zlen * 0.5) {
							zj -= zlen;
						}


						double rc = sqrt((xatom-xj)*(xatom-xj)+(yatom-yj)*(yatom-yj)+(zatom-zj)*(zatom-zj));

						if(cutoff > rc || cellindex_1D==atom_cellindex_1D) { //cubic cell is within the cut-off radius

							if(particle_neighbor_count_hs[j][cellindex_1D]==0 ){

								cavitycount_indices[ncavity]=cellindex_1D;
								++ncavity;

								int iseg = int(cellindex_1D/num_cells_in_segment);
								particle_types.updateNumCavitiesInSegments(j,iseg,-1);

							}

							++particle_neighbor_count_hs[j][cellindex_1D];

						} // end of IF statement (cutoff > rc)

					}

				}


			}

			if (ncavity > 0){

				for (int ind4 = 0; ind4 < ncavity; ind4++){
					cellindex_by_cavity_countindex[j][cavitycount_indices[ind4]] = -1;
				}

				int num = particle_types.getNumCavities(j);
				num -= ncavity;
				particle_types.updateNumCavities(j,num);

			} // end of IF statement (ncavity > 0)
		}
	}

	cavitycount_indices.clear();
}
