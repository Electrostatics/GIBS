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
 * update_cavity_list.cpp
 *
 *  Created on: Jun 8, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file update_cavity_list.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief CParticleState member function definitions for updating the cavity list of each particle type in
 *  		the simulation box
 *
 */

#include "particle_state.hpp"

/**
 * @brief Updates the cavity list of each particle type after a state change due to a GCMC move
 *
 * @param[in]	ptcl(struct_particle)	 Data structure of the particle selected for insertion, deletion, or single-particle displacement
 * @param[in]	particle_types(CParticleType_t)		Class  instance of the selected particle type
 * @param[in]	particlepair_types(CParticlePairType_t)	Class  instance of the selected particle pair type
 * @param[in]	cavity_grid(CCavityGrid_t)		Class  instance containing the cavity grid data
 * @param[in]	state_change(enum_state_change)	Type of state change
 *
 * @return void
 *
 */
void CParticleState::updateCavityList(struct_particle ptcl,CParticleType_t &particle_types,
		const CParticlePairType_t &particlepair_types,const CCavityGrid_t &cavity_grid,enum_state_change state_change){

	/*If particles are inserted (state_change: insertion), then some sites for each
	 * particle type will be deleted from its  cavity list.
	 * If particles are (state_change: deletion), then some sites for each particle
	 * type may be added.
	 */

	// current cell index

	int nparticle_type = particle_types.getNumParticleTypes();

	int lix = ptcl.gridcell_ix;
	int liy = ptcl.gridcell_iy;
	int liz = ptcl.gridcell_iz;

	double x = ptcl.x;
	double y = ptcl.y;
	double z = ptcl.z;

	double xlen = cavity_grid.getXLen();
	double ylen = cavity_grid.getYLen();
	double zlen = cavity_grid.getZLen();

	int ptcl_cellindex_1D = ptcl.gridcell_1Dindex;

	int num_cells_in_segment = cavity_grid.getNumCellsInSegment();
	for (int j=0; j<nparticle_type;j++){

		std::string label12 = ptcl.label;
		label12.append("_");
		label12.append(particle_types.getLabel(j));

		struct_cellindex_offset cellindex_offset = particlepair_types.getHardSphereCellIndexOffsetMatrix(label12);
		int pair_index = particlepair_types.getIndexFromLabel(label12);

		double hs_cutoff = particlepair_types.getHardSphereCutoff(label12);

		int ncx = cellindex_offset.ncx;
		int ncy = cellindex_offset.ncy;
		int ncz = cellindex_offset.ncz;

		int ljx,ljy,ljz;

		int ncavity = 0;

		for (int ind1=0;ind1<2*ncx+1;ind1++){

			ljx = lix + cellindex_offset.delxj[ind1] + cellindex_offset.bcmap_x[lix][ind1];

			for (int ind2 = 0; ind2<2*ncy+1; ind2++){
				ljy = liy + cellindex_offset.delyj[ind2] + cellindex_offset.bcmap_y[liy][ind2];

				for (int ind3 = 0; ind3 < 2*ncz+1; ind3++){

					ljz = liz + cellindex_offset.delzj[ind3] + cellindex_offset.bcmap_z[liz][ind3];

					int cellindex_1D = cavity_grid.getCellIndex1D(ljx,ljy,ljz) ;

					// check if the cell is enclosed within the hard sphere cutoff radius.

					// get the center of the cavity grid cell
					double	xj = cavity_grid.getXCenter(ljx);
					double	yj = cavity_grid.getYCenter(ljy);
					double	zj = cavity_grid.getZCenter(ljz);

					// Find the nearest image
					if (xj < x - xlen * 0.5) {
						xj += xlen;
					} else if (xj > x + xlen * 0.5) {
						xj -= xlen;
					}
					if (yj < y - ylen * 0.5) {
						yj += ylen;
					} else if (yj > y + ylen * 0.5) {
						yj -= ylen;
					}
					if (zj < z - zlen * 0.5) {
						zj += zlen;
					} else if (zj > z + zlen * 0.5) {
						zj -= zlen;
					}

					double rc = sqrt((x-xj)*(x-xj)+(y-yj)*(y-yj)+(z-zj)*(z-zj));

					if (hs_cutoff > rc || cellindex_1D ==ptcl_cellindex_1D) { //square cell is within the cut-off radius

						if(particle_neighbor_count_hs[j][cellindex_1D]==0 && state_change==insertion){

							cavitycount[pair_index].indices[ncavity]=cellindex_1D;
							++ncavity;

							int iseg = int(cellindex_1D/num_cells_in_segment);
							particle_types.updateNumCavitiesInSegments(j,iseg,-1);

						}
						if(state_change==insertion){
							++particle_neighbor_count_hs[j][cellindex_1D];

						}else if (state_change == deletion){

							--particle_neighbor_count_hs[j][cellindex_1D];

							if(particle_neighbor_count_hs[j][cellindex_1D]==0){
								int num = particle_types.getNumCavities(j);
								++num;
								particle_types.updateNumCavities(j,num);

								cellindex_by_cavity_countindex[j][cellindex_1D] = cellindex_1D;

								int iseg = int(cellindex_1D/num_cells_in_segment);
								particle_types.updateNumCavitiesInSegments(j,iseg,1);

							}
						}

					}  // end of if statement ( rc > hs_cutoff)

				}
			}
		}


		if (state_change==insertion && ncavity > 0){

			// The number of indices in the cavity list
			int num = particle_types.getNumCavities(j);

			for (int ind4 = 0; ind4 < ncavity; ind4++){
				cellindex_by_cavity_countindex[j][cavitycount[pair_index].indices[ind4]] = -1;
			}

			num -= ncavity;
			particle_types.updateNumCavities(j,num);

		}

	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
