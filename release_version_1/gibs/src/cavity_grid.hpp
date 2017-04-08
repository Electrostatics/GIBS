/*
 * cavity_grid.hpp
 *
 *  Created on: Apr 13, 2015
 *      Author:  Dennis G. Thomas
 *
 * @file cavity_grid.hpp
 * @author Dennis G. Thomas
 *
 * @brief Class declaration, struct definitions, and function prototypes for cavity grid
 */

#ifndef CAVITY_GRID_HPP_
#define CAVITY_GRID_HPP_

#include "generic.hpp"
#include "box.hpp"

/**
 * @brief Struct data type for storing the cavity grid cell indices along x,y, and z dimensions
 *
 * The cavity grid is a 3-dimensional grid with uniform spacing along x, y, and z-axes of the simulation box.
 * The number of grid points along x, y, and z-axes are denoted as nx, ny, and nz
 * The x-axis grid points are stored in an array x of dimension nx.
 * The y-axis grid points are stored in an array y of dimension ny.
 * The z-axis grid points are stored in an array z of dimension nz.
 * The grid cell indices are stored as a 3-dimensional and as 1-dimensional array.
 * The 3D array is of size nx, ny and nz. The 1-D array is of size nx*ny*nz.
 * The 1D cell indices are grouped into segments
 *
 *	Use defined type struct_cavitygrid_cellindex_coords for struc declarations
 */
struct cavitygrid_cellindex_coords_struct{

	int icx; // Integer variable to store the index of a cavity grid cell along x-axis
	int icy; // Integer variable to store the index of a cavity grid cell along y-axis
	int icz;  //Integer variable to store the  index of a cavity grid cell along z-axis

	// Values of icx: 0,1,...,Ncx-1
	// Values of icy: 0,1,...,Ncy-1
	// Values of icz: 0,1,...,Ncz-1

};
typedef struct cavitygrid_cellindex_coords_struct struct_cavitygrid_cellindex_coords;

/**
 * @brief Struct data type for neighboring cells along x,y, and z dimensions
 *
 * Use defined type struct_cellindex_offset for struct declarations
 */

struct cellindex_offset_struct{

	/* String variable to store the particle pair labels
	(e.g., Na_Cl for sodium chloride ion pair, Cl_Cl for chloride ion pair)*/
	std::string particlepair_label;
	// Floating point variable to store any potential's cut-off distance
	double rcutoff;
	// Floating point variables to store the grid cell edge distances along each axis dimension
	double hcx,hcy,hcz;
	// Integer variable to store the array index of the particle pair
	int particlepair_index;

	// Integer variables to store the number of grid cells along each axis dimension
	int Ncx,Ncy,Ncz;
	/* Integer variables to store the number of neigboring cells from a given cell index
	  within the cut-off distance "rcutoff" */
	int ncx,ncy,ncz;
	int nindices;
	int *delxj;
	int *delyj;
	int *delzj;
	int **bcmap_x;
	int **bcmap_y;
	int **bcmap_z;

	/* bcmap_x, bcmap_y, and bcmap_z take values -1, 0, or 1:
		-1 if central cell index is > Ncx-1-ncx (similarly for y and z)
		1 if central cell index < ncx (similarly for y and z)
		 0 , otherwise.
	 */
};
typedef struct cellindex_offset_struct struct_cellindex_offset;

/**
 * @brief Class defining the parameters, arrays, and matrices associated with cavity grid
 *
 * Use defined type CCavityGrid_t for class declarations
 */
class CCavityGrid{

	/* Integer variables to store the number of cavity grid points along x, y, and z axes */
	int nx,ny,nz;
	/* Integer variable store the number of grid cells along each dimension */
	int Ncx,Ncy,Ncz;
	/*Floating point variable to store the box length (in Angstrom units) along x, y, and z axes */
	double xlen,ylen,zlen;
	/* Floating point variables to store the uniform grid spacing values along x, y, and z axes*/
	double hx,hy,hz;
	/*Floating point value to store the the lowest value of x, y, and z */
	double xmin,ymin,zmin;

	/*Allocatable array of floating point values to store the x-coordinate values of the grid points along x-axis */
	double *x;
	/*Allocatable array of floating point values to store the y-coordinate values of the grid points along y-axis */
	double *y;
	/*Allocatable array of floating point values to store the z-coordinate values of the grid points along z-axis */
	double *z;

	/*Allocatable array of floating point values to store the x-coordinates of the center of each grid cell along x-axis */
	double *xcenter;
	/*Allocatable array of floating point values to store the y-coordinates of the center of each grid cell along y-axis */
	double *ycenter;
	/*Allocatable array of floating point values to store the z-coordinates of the center of each grid cell along z-axis */
	double *zcenter;


	/*Integer variable to store the number of cavity grid cells*/
	int num_cells;
	/*Integer variable to store the number of cavity grid cell segments */
	int num_cavity_segments;
	/*Integer variable to store the number of cavity grid cells in each segment*/
	int num_cells_in_segment;
	/*Integer variable to store the index of the first cell of each segment */
	int *first_cellindex_of_segment;

	/*3-dimensional array to map the 3D cell indices to 1D cell indices of the cavity grid*/
	int ***cellindex_3Dto1Dmap;
	/*A map variable to store the x,y,z indices of each cell, given the cell's 1D index */
	std::tr1::unordered_map<int, struct_cavitygrid_cellindex_coords> cellindex_1Dto3Dmap;

	/*A map variable to associate the LJ grid 1D cell index to each cavity grid 1D cell index*/
	std::tr1::unordered_map<int, int> cavity_to_LJ_cellindex_1Dmap;
	/*Allocatable array of floating point values to store the x-coordinate value of the center of
	 * of each LJ grid cell */
	double *xcenter_ljgrid;
	/*Allocatable array of floating point values to store the y-coordinate value of the center of
		 * of each LJ grid cell */
	double *ycenter_ljgrid;
	/*Allocatable array of floating point values to store the z-coordinate value of the center of
		 * of each LJ grid cell */
	double *zcenter_ljgrid;

	struct_cellindex_offset lj_cellindex_offset;

	/*3-dimensional array to map the 3D cell indices to 1D cell indices of LJ grid */
	int ***ljcellindex_3Dto1Dmap;

	/*A map variable to store the x,y,z indices of each LJ grid cell, given the cell's 1D index */
	std::tr1::unordered_map<int, struct_cavitygrid_cellindex_coords> ljcellindex_1Dto3Dmap;

	std::tr1::unordered_map<int, std::vector<int> > ljcellindex_1D_neighbor;

	/*Number of LJ grid cells */
	int num_lj_cells;
	/*Number of cavity grid cells in each LJ cell */
	int num_cavitycells_perLJcell;

	void setCellIndexMaps();
public:

	// constructor
	CCavityGrid(const CBox_t &box,double grid_spacing,int num_cavity_segments);

	// destructor
	~CCavityGrid(){

		delete[] first_cellindex_of_segment;

		delete[] x;
		delete[] y;
		delete[] z;

		x = NULL;
		y = NULL;
		z = NULL;

		delete[] xcenter;
		delete[] ycenter;
		delete[] zcenter;

		xcenter = NULL;
		ycenter = NULL;
		zcenter = NULL;

		for (int ind1 = 0; ind1 < Ncx; ind1++) {
			for (int ind2 = 0; ind2 < Ncy; ind2++) {
				delete[] cellindex_3Dto1Dmap[ind1][ind2];
			}
			delete[] cellindex_3Dto1Dmap[ind1];
		}
		delete[] cellindex_3Dto1Dmap;
		cellindex_3Dto1Dmap = NULL;

		std::cout << "[cavity_grid.hpp(~CCavityGrid]: called CCavityGrid destructor." << std::endl;
	}

	void setLJCellIndexMaps(double lj_cutoff);

	std::vector<int> getLJNeighborCellIndices(int lj_cellindex_1D) const{
		return this->ljcellindex_1D_neighbor.find(lj_cellindex_1D)->second;
	}


	int getLJCellindexByCavityCellIndex(int cellindex_1D) const{
		return this->cavity_to_LJ_cellindex_1Dmap.find(cellindex_1D)->second;
	}

	int getLJCellindexOffsetNindices() const{
		return this->lj_cellindex_offset.nindices;
	}

	int getCavityGridNumCells() const{
		return this->num_cells;
	}

	int getCavityGridCellsPerLJCell() const{
		return this->num_cavitycells_perLJcell;
	}

	int getLJGridNumCells() const{
		return this->num_lj_cells;
	}

	int getNx()const {return nx;}
	int getNy()const{return ny;}
	int getNz()const {return nz;}
	double getXLen()const {return xlen;}
	double getYLen()const {return ylen;}
	double getZLen()const {return zlen;}

	double getHx()const {return hx;}
	double getHy()const {return hy;}
	double getHz()const {return hz;}
	double getXMin()const {return xmin;}
	double getYMin()const {return ymin;}
	double getZMin()const {return zmin;}

	int getNcx()const {return Ncx;}
	int getNcy()const{return Ncy;}
	int getNcz()const {return Ncz;}

	double getX(int i)const {return x[i];}
	double getY(int i)const {return y[i];}
	double getZ(int i)const {return z[i];}

	double getXCenter(int i)const {return xcenter[i];}
	double getYCenter(int i)const {return ycenter[i];}
	double getZCenter(int i)const {return zcenter[i];}

	int getCellIndex1D(int i, int j, int k)const {return cellindex_3Dto1Dmap[i][j][k];}
	struct_cavitygrid_cellindex_coords getCellIndex3D(int i)const {return cellindex_1Dto3Dmap.find(i)->second;}

	int getNumCellsInSegment()const {return num_cells_in_segment;}
	int getNumCavitySegments()const {return num_cavity_segments;}
	int getFirstCellIndexOfSegment(int iseg)const{return first_cellindex_of_segment[iseg];}
};
typedef class CCavityGrid CCavityGrid_t;

#endif /* CAVITY_GRID_HPP_ */
