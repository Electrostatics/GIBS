/*
 * particle_state.hpp
 *
 *  Created on: May 11, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file particle_state.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Class definitions and function prototypes for particle state class
 */

#ifndef PARTICLE_STATE_HPP_
#define PARTICLE_STATE_HPP_

#include "generic.hpp"
#include "enum_types.hpp"
#include "cavity_grid.hpp"
#include "box.hpp"
#include "particle_type.hpp"
#include "particlepair_type.hpp"
#include "ion_type.hpp"
#include "particle_structs.hpp"
#include "solute_type.hpp"

/**
 * @brief Data structure to store the list of cavity indices for each
 * 			particle type in the simulation box
 *
 */
struct cavitycount_struct{
	std::vector<int> indices;
};
typedef struct cavitycount_struct struct_cavitycount;

/**
 * @brief Class for storing data of the particle state of the system
 *
 * Use defined type CParticleState_t for class declarations
 */
class CParticleState{

	std::vector<struct_particle> particles;
	int npart;
	int nparticle_types;
	int nparticle_pairs;
	int num_cells;
	int **particle_neighbor_count_hs;
	int **particle_indices_by_count;
	int **particle_indices_by_cellindex;

	int **cellindex_by_cavity_countindex;

	struct_cavitycount *cavitycount;

public:
	// class constructor
	CParticleState(int nparticle_types,const CCavityGrid_t &cavity_grid,const CParticlePairType_t &particlepair_types,int max_npart);

	// class destructor
	~CParticleState(){
		for (int ind1=0; ind1 < nparticle_types; ind1++){
			delete[] cellindex_by_cavity_countindex[ind1];
			delete[] particle_indices_by_count[ind1];
			delete[] particle_indices_by_cellindex[ind1];
			delete[] particle_neighbor_count_hs[ind1];

		}
		delete[] cellindex_by_cavity_countindex;
		delete[] particle_indices_by_count;
		delete[] particle_indices_by_cellindex;
		delete[] particle_neighbor_count_hs;

		cellindex_by_cavity_countindex = NULL;
		particle_indices_by_count = NULL;
		particle_indices_by_cellindex = NULL;
		particle_neighbor_count_hs = NULL;

		for (int i=0;i<nparticle_pairs;i++){
			cavitycount[i].indices.clear();
		}
		delete[] cavitycount;
		cavitycount = NULL;

	}



	// initialization methods
	void initializeParticlePositions_WithoutGridBasedInsertion(const CBox_t &box,double min_ptcl_spacing,
			int niontype);
	void initializeParticleTypeWithIonParameters(CIonType_t &iontypes,
			double solute_charge, int niontype,CParticleType_t &particle_types, int flag);
	void initializeParticleType(CParticleType_t &particle_types,double solute_charge, int flag);
	void initializeParticleState_WithGridBasedInsertion(CParticleType_t &particle_types,
			const CCavityGrid_t &cavity_grid,CParticlePairType_t &particlepair_types,int target_num[]);
	void readInitializeParticleState_WithGridBasedInsertion(CParticleType_t &particle_types,
			const CCavityGrid_t &cavity_grid,CParticlePairType_t &particlepair_types,int file_num, int num_lines_skip,
			std::string dir_name,int target_num[]);
	void readFromFile(int file_num, int num_lines_skip, double &energy_old,std::string dir_name) ;

	void countParticleTypeFromState(CParticleType_t &particle_types);
	// update methods
	void updateParticleNum(int n){npart=n;}
	void updateParticleIndicesByCount(int itype,int countindex, int val)
	{particle_indices_by_count[itype][countindex] = val; }

	void updateParticleIndicesByCellindex(int itype,int cellindex, int val)
	{particle_indices_by_cellindex[itype][cellindex] = val; }
	void updateParticleNeighborCountHS(int itype,int cellindex, int val)
	{particle_neighbor_count_hs[itype][cellindex] = val; }
	void updateCavityList(struct_particle ptcl,CParticleType_t &particle_types,
			const CParticlePairType_t &particlepair_types,const CCavityGrid_t &cavity_grid,enum_state_change state_change);


	void incrementParticleCountByOne(){ ++npart;}
	void decrementParticleCountByOne(){ --npart;}

	void updateCavityGridVariables_AllParticles(CParticleType_t &particle_types,
			const CCavityGrid_t &cavity_grid,const CParticlePairType_t &particlepair_types);
	void updateCavityGridVariables_SoluteAllAtom(CParticleType_t &particle_types,
			const CCavityGrid_t &cavity_grid,const CSoluteModel_t &solute_model);

	void addParticle(struct_particle ptcl,int ptcltype_index,int num,int cellindex_1D);

	void revertToParticleToPreviousPosition(struct_particle ptcl,int ptclindex, int ptcltype_index,int cellindex_1D){
		particles[ptclindex] = ptcl;
		particle_indices_by_cellindex[ptcltype_index][ptcl.gridcell_1Dindex] = ptclindex;
		particle_indices_by_cellindex[ptcltype_index][cellindex_1D] = -1;
	}

	void deleteParticle(int cellindex_1D,int ptclcount_index,int ptcltype_index,
			int ptclindex,const CParticleType_t &particle_types);

	void copyParticlePosition(int i,struct_particle ptcl);

	// get methods
	int getParticleNum() const {return npart;}


	struct_particle getParticle(int i) const {return particles[i];}
	int getParticleIndicesByCount(int itype,int countindex)const {return particle_indices_by_count[itype][countindex];}
	int getParticleIndicesByCellindex(int itype,int cellindex)const {return particle_indices_by_cellindex[itype][cellindex];}
	int getParticleNeighborCountHS(int itype,int cellindex)const {return particle_neighbor_count_hs[itype][cellindex];}
	int getCellIndexByCavityCountIndex(int itype,int i)const {return cellindex_by_cavity_countindex[itype][i];}
	//int getCavityCountIndexByCellIndex(int itype,int i)const {return cavity_countindex_by_cellindex[itype][i];}
	double getParticleCharge(int i) const { return particles[i].charge;}
	double getParticlePositionX(int i)const {return particles[i].x;}
	double getParticlePositionY(int i)const {return particles[i].y;}
	double getParticlePositionZ(int i) const {return particles[i].z;}
	int getParticleGridCellIndexX (int i) const {return particles[i].gridcell_ix;}
	int getParticleGridCellIndexY (int i) const {return particles[i].gridcell_iy;}
	int getParticleGridCellIndexZ (int i) const {return particles[i].gridcell_iz;}
	int getParticleGridCellIndex1D (int i) const {return particles[i].gridcell_1Dindex;}

	int searchParticleWithinHSCutoff(struct_particle ptcl,CParticleType_t &particle_types,
			const CParticlePairType_t &particlepair_types,const CCavityGrid_t &cavity_grid);

	// write methods
	void writeToFile(int istep, int file_index,const CBox_t &box, double energy, std::string filename);



};

typedef class CParticleState CParticleState_t;



#endif /* PARTICLE_STATE_HPP_ */
