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
 * cavity_grid.cpp
 *
 *  Created on: Apr 13, 2015
 *      Author:  Dennis G. Thomas
 *
 * @file cavity_grid.cpp
 * @author Dennis G. Thomas
 *
 *  @brief	CCavityGrid class member function definitions
 *
 */


#include "cavity_grid.hpp"

/**
 * @brief	CCavityGrid class constructor
 *
 * This function sets up the cavity grid dimensions, calls setCellIndexMaps to set the
 * arrays for mapping 1D to 3D cell indices and vice-versa, and the number of
 * cavity grid segments. The total number of cells should be divisible by the number of cavity segments, so that
 * they can be equally divided among the segments.
 *
 * @param[in,out] box(CBox_t)				Simulation box parameters
 * @param[in] grid_spacing(double)			Cavity grid spacing in Angstrom units
 * @param[in] num_cavity_segments(int) 		Number of cavity grid segments
 *
 */

CCavityGrid::CCavityGrid(const CBox_t &box,double grid_spacing,int num_cavity_segments){

	xmin = box.getBoxXMin();
	ymin = box.getBoxYMin();
	zmin = box.getBoxZMin();

	xlen = box.getBoxXLen();
	ylen = box.getBoxYLen();
	zlen = box.getBoxZLen();


	if(fmod(xlen,grid_spacing) ==0){
		nx = (int)xlen/grid_spacing + 1;
		hx = grid_spacing;
	}else {
		nx = int(xlen/grid_spacing)+2;
		hx = xlen/(nx-1);
	}


	if(fmod(ylen,grid_spacing) ==0){
		ny = (int)ylen/grid_spacing + 1;
		hy = grid_spacing;
	}else {
		ny = int(ylen/grid_spacing)+2;
		hy = ylen/(ny-1);
	}


	if(fmod(zlen,grid_spacing) ==0){
		nz = (int)zlen/grid_spacing + 1;
		hz = grid_spacing;
	}else {
		nz = int(zlen/grid_spacing)+2;
		hz = zlen/(nz-1);
	}


	std::cout << "[cavity_grid.cpp(CCavityGrid::CCavityGrid)]: cavity grid spacing in x = " << hx << std::endl;
	std::cout << "[cavity_grid.cpp(CCavityGrid::CCavityGrid)]: cavity grid spacing in y = " << hy << std::endl;
	std::cout << "[cavity_grid.cpp(CCavityGrid::CCavityGrid)]: cavity grid spacing in z = " << hz << std::endl;

	// set up cavity grid points
	x = new double[nx];
	y = new double[ny];
	z = new double[nz];

	for (int i=0; i<nx; i++){
		x[i] = xmin+i*hx;
		//std::cout << "[CCavityGrid::CCavityGrid] i = " << i << ", x = " << x[i] << std::endl;
	}

	for (int i=0; i<ny; i++){
		y[i] = ymin+i*hy;
		//	std::cout << "[CCavityGrid::CCavityGrid] i = " << i << ", y = " << y[i] << std::endl;
	}

	for (int i=0; i<nz; i++){
		z[i] = zmin+i*hz;
		//	std::cout << "[CCavityGrid::CCavityGrid] i = " << i << ", z = " << z[i] << std::endl;
	}

	Ncx = nx - 1;
	Ncy = ny - 1;
	Ncz = nz - 1;

	std::cout << "[cavity_grid.cpp(CCavityGrid::CCavityGrid)]: number of cavity grid cells along x = " << Ncx << std::endl;
	std::cout << "[cavity_grid.cpp(CCavityGrid::CCavityGrid)]: number of cavity grid cells along y = " << Ncy << std::endl;
	std::cout << "[cavity_grid.cpp(CCavityGrid::CCavityGrid)]: number of cavity grid cells along z = " << Ncz << std::endl;

	// set up the centers of cavity grid cells

	xcenter = new double[Ncx];
	ycenter = new double[Ncy];
	zcenter = new double[Ncz];

	for (int i=0; i<Ncx; i++){
		xcenter[i] =x[i] + 0.5*hx;

	}

	for (int i=0; i<Ncy; i++){
		ycenter[i] = y[i]+ 0.5*hy;

	}

	for (int i=0; i<Ncz; i++){
		zcenter[i] = z[i] + 0.5*hz;
	}

	// Set cell index maps
	setCellIndexMaps();

	this->num_cells = Ncx*Ncy*Ncz;

	/*get the number of segments */
	/*Total number of cells should be divisible by the number of cavity segements, so that
		they can be equally divided among the segments.*/
	this->num_cavity_segments = num_cavity_segments;


	this->num_cells_in_segment = this->num_cells/num_cavity_segments;
	this->first_cellindex_of_segment = new int[num_cavity_segments];

	for (int i=0; i < num_cavity_segments; i++){
		this->first_cellindex_of_segment[i] = i*this->num_cells_in_segment;

	}

	this->num_lj_cells = 0;// number of LJ cells = 0 if LJ is not used.

}



/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/**
 * @brief Sets the array for mapping the 1D and 3D cell indices of the cavity grid
 *
 * @return void
 */
void CCavityGrid::setCellIndexMaps(){

	int count = 0;
	struct_cavitygrid_cellindex_coords coords;

	cellindex_3Dto1Dmap = new int**[Ncx];
	for (int ind1 = 0; ind1 < Ncx; ind1++) {
		cellindex_3Dto1Dmap[ind1] = new int*[Ncy];
		for (int ind2 = 0; ind2 < Ncy; ind2++) {
			cellindex_3Dto1Dmap[ind1][ind2] = new int[Ncz];
		}
	}

	for (int i=0; i< Ncx;i++){
		for(int j=0; j < Ncy; j++){
			for (int k=0; k<Ncz; k++){

				coords.icx = i;
				coords.icy=j;
				coords.icz =k;


				cellindex_1Dto3Dmap.insert(std::make_pair(count,coords));

				cellindex_3Dto1Dmap[i][j][k] = count;

				++count;

			}
		}
	}


}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Sets up the cell index maps for Lennard Jones cells
 *
 * This functions sets the arrays for mapping the cavity grid cell indices
 * to Lennard-Jones (L-J) grid cell indices.
 *
 * @param[in] lj_cutoff(double)	Lennard-Jones potential cut-off distance
 *
 * @return void
 */
void CCavityGrid::setLJCellIndexMaps(double lj_cutoff){

	double xlen = this->xlen;
	double ylen = this->ylen;
	double zlen =this->zlen;

	double hcx = xlen/floor(xlen/(lj_cutoff/sqrt(3)));
	double hcy = ylen/floor(ylen/(lj_cutoff/sqrt(3)));
	double hcz = zlen/floor(zlen/(lj_cutoff/sqrt(3)));

	int Ncx, Ncy, Ncz;

	if(fmod(xlen,hcx) ==0){
		Ncx = (int)xlen/hcx;
		std::cout << "[CCavityGrid::setLJCellIndexMaps]:  grid spacing in x = " << hcx << std::endl;

	}else {
		Ncx = int(xlen/hcx);
		hcx = (double) xlen/Ncx;
		//std::cerr << "[CCavityGrid::setLJCellIndexMaps]: Ncx is not an integer." << std::endl;
		std::cout << "[CCavityGrid::setLJCellIndexMaps]:  grid spacing in x = " << hcx << std::endl;

	}


	if(fmod(ylen,hcy) ==0){
		Ncy = (int)ylen/hcy;
		std::cout << "[CCavityGrid::setLJCellIndexMaps]:  grid spacing in y = " << hcy<< std::endl;

	}else {
		Ncy = int(ylen/hcy);
		hcy = (double) ylen/Ncy;
		//	std::cerr << "[CCavityGrid::setLJCellIndexMapsx]: Ncy is not an integer." << std::endl;
		std::cout << "[CCavityGrid::setLJCellIndexMaps]:  grid spacing in y = " << hcy<< std::endl;

	}
	if(fmod(zlen,hcz) ==0){
		Ncz = (int)zlen/hcz;
		std::cout << "[CCavityGrid::setLJCellIndexMaps]:  grid spacing in z = " << hcz << std::endl;

	}else {
		Ncz = int(zlen/hcz);
		hcz = (double) zlen/Ncz;
		//	std::cerr << "[CCavityGrid::setLJCellIndexMaps]: Ncz is not an integer." << std::endl;
		std::cout << "[CCavityGrid::setLJCellIndexMaps]:  grid spacing in z = " << hcz << std::endl;

	}

	int ncx = int (floor(lj_cutoff/hcx));
	int ncy = int (floor(lj_cutoff/hcy));
	int ncz = int (floor(lj_cutoff/hcz));


	this->lj_cellindex_offset.Ncx = Ncx;
	this->lj_cellindex_offset.Ncy = Ncy;
	this->lj_cellindex_offset.Ncz = Ncz;
	this->lj_cellindex_offset.hcx = hcx;
	this->lj_cellindex_offset.hcy = hcy;
	this->lj_cellindex_offset.hcz = hcz;
	this->lj_cellindex_offset.ncx = ncx;
	this->lj_cellindex_offset.ncy = ncy;
	this->lj_cellindex_offset.ncz = ncz;

	this->lj_cellindex_offset.rcutoff = lj_cutoff;

	// map LJ grid cell indices to cavity grid cell indices

	// set up LJ grid points

	this->xcenter_ljgrid = new double[Ncx];
	this->ycenter_ljgrid = new double[Ncy];
	this->zcenter_ljgrid = new double[Ncz];

	for (int i=0; i<Ncx; i++){
		this->xcenter_ljgrid[i] =this->xmin + 0.5*hcx;

	}

	for (int i=0; i<Ncy; i++){
		this->ycenter_ljgrid[i] =this->ymin + 0.5*hcy;

	}

	for (int i=0; i<Ncz; i++){
		this->zcenter_ljgrid[i] =this->zmin + 0.5*hcz;

	}

	// set  1D  to 3D cell index maps for LJ grid and vice-versa.
	int count = 0;
	struct_cavitygrid_cellindex_coords coords;

	this->ljcellindex_3Dto1Dmap = new int**[Ncx];
	for (int ind1 = 0; ind1 < Ncx; ind1++) {
		this->ljcellindex_3Dto1Dmap[ind1] = new int*[Ncy];
		for (int ind2 = 0; ind2 < Ncy; ind2++) {
			this->ljcellindex_3Dto1Dmap[ind1][ind2] = new int[Ncz];
		}
	}

	for (int i=0; i< Ncx;i++){
		for(int j=0; j < Ncy; j++){
			for (int k=0; k<Ncz; k++){

				coords.icx = i;
				coords.icy=j;
				coords.icz =k;

				this->ljcellindex_1Dto3Dmap.insert(std::make_pair(count,coords));
				this->ljcellindex_3Dto1Dmap[i][j][k] = count;
				++count;

			}
		}
	}


	this->num_cavitycells_perLJcell = lj_cutoff/this->hx; // this is a rough estimate, assume hx = hy = hz

	// done setting cell index maps for LJ grid

	// mapping 1D cellindices of cavity grid and LJ cutoff grid
	for (int i=0; i < this->Ncx; i++){

		double xc = this->xcenter[i];

		int icx = (int) floor((xc-this->xmin)/hcx);
		for (int j = 0; j < this->Ncy; j++){

			double yc = this->ycenter[j];
			int icy = (int) floor((yc - this->ymin)/hcy);

			for(int k = 0; k < this->Ncz; k++){
				double zc = this->zcenter[k];
				int icz = (int) floor((zc - this->zmin)/hcz);

				int cellindex_1D = this->cellindex_3Dto1Dmap[i][j][k];
				int lj_cellindex_1D = this->ljcellindex_3Dto1Dmap[icx][icy][icz];
				this->cavity_to_LJ_cellindex_1Dmap.insert(std::make_pair(cellindex_1D, lj_cellindex_1D));

			}
		}
	}

	// set cellindex offset distances and matrices

	this->lj_cellindex_offset.delxj = new int[2*ncx+1];
	this->lj_cellindex_offset.delyj = new int[2*ncy+1];
	this->lj_cellindex_offset.delzj = new int[2*ncz+1];

	for(int i=0;i<2*ncx+1;i++){
		this->lj_cellindex_offset.delxj[i]=-ncx+i;
	}
	for(int i=0;i<2*ncy+1;i++){
		this->lj_cellindex_offset.delyj[i]=-ncy+i;
	}
	for(int i=0;i<2*ncz+1;i++){
		this->lj_cellindex_offset.delzj[i]=-ncz+i;
	}

	// map periodic boundary conditions on cell neighbor indices.
	this->lj_cellindex_offset.bcmap_x = new int*[Ncx];
	for (int i=0;i<Ncx;i++){
		this->lj_cellindex_offset.bcmap_x[i] = new int[2*ncx+1];
	}
	this->lj_cellindex_offset.bcmap_y = new int*[Ncy];
	for (int i=0;i<Ncy;i++){
		this->lj_cellindex_offset.bcmap_y[i] = new int[2*ncy+1];
	}
	this->lj_cellindex_offset.bcmap_z = new int*[Ncz];
	for (int i=0;i<Ncz;i++){
		this->lj_cellindex_offset.bcmap_z[i] = new int[2*ncz+1];
	}

	for(int i=0; i < Ncx; i++){
		for (int j=0;j<2*ncx+1;j++){
			this->lj_cellindex_offset.bcmap_x[i][j] = 0;
			if(i+this->lj_cellindex_offset.delxj[j] > Ncx-1){
				this->lj_cellindex_offset.bcmap_x[i][j] = -Ncx;
			}else if(i+this->lj_cellindex_offset.delxj[j] < 0){
				this->lj_cellindex_offset.bcmap_x[i][j] = Ncx;
			}

		}

	}


	for(int i=0; i < Ncy; i++){
		for (int j=0;j<2*ncy+1;j++){
			this->lj_cellindex_offset.bcmap_y[i][j] = 0;
			if(i+this->lj_cellindex_offset.delyj[j] > Ncy-1){
				this->lj_cellindex_offset.bcmap_y[i][j] = -Ncy;
			}else if(i+this->lj_cellindex_offset.delyj[j] < 0){
				this->lj_cellindex_offset.bcmap_y[i][j] = Ncy;
			}

		}

	}

	for(int i=0; i < Ncz; i++){
		for (int j=0;j<2*ncz+1;j++){
			this->lj_cellindex_offset.bcmap_z[i][j] = 0;
			if(i+this->lj_cellindex_offset.delzj[j] > Ncz-1){
				this->lj_cellindex_offset.bcmap_z[i][j] = -Ncz;
			}else if(i+this->lj_cellindex_offset.delzj[j] < 0){
				this->lj_cellindex_offset.bcmap_z[i][j] = Ncz;
			}

		}

	}

	this->lj_cellindex_offset.nindices = (2*ncx+1)*(2*ncy+1)*(2*ncz+1);


	std::vector<int> neighbor_indices;
	neighbor_indices.resize(this->lj_cellindex_offset.nindices);
	this->num_lj_cells = Ncx*Ncy*Ncz;
	// creating 1D neighbor list of LJ cellindices for each LJ cell index.
	for (int i=0;i<this->num_lj_cells;i++){

		struct_cavitygrid_cellindex_coords cellindex3D_coords = this->ljcellindex_1Dto3Dmap.find(i)->second;
		int ljx,ljy,ljz;
		int count = 0;
		for (int ind3=0;ind3<2*ncx+1;ind3++){
			ljx = cellindex3D_coords.icx + this->lj_cellindex_offset.delxj[ind3] + this->lj_cellindex_offset.bcmap_x[cellindex3D_coords.icx][ind3];


			for (int ind4 = 0; ind4<2*ncy+1; ind4++){
				ljy = cellindex3D_coords.icy + this->lj_cellindex_offset.delyj[ind4] + this->lj_cellindex_offset.bcmap_y[cellindex3D_coords.icy][ind4];

				for (int ind5 = 0; ind5 < 2*ncz+1; ind5++){
					ljz = cellindex3D_coords.icz + this->lj_cellindex_offset.delzj[ind5] + this->lj_cellindex_offset.bcmap_z[cellindex3D_coords.icz][ind5];

					int cellindex_1D = this->ljcellindex_3Dto1Dmap[ljx][ljy][ljz];
					neighbor_indices[count] = cellindex_1D;
					++count;
				}
			}
		}

		this->ljcellindex_1D_neighbor.insert(std::make_pair(i,neighbor_indices));

	}

}
