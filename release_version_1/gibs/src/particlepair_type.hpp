/*
 * particlepair_type.hpp
 *
 *  Created on: Apr 25, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file particlepair_type.hpp
 *  @author Dennis G. Thomas
 *
 * @brief Class definitions and function prototypes for particle pair type class
 */

#ifndef PARTICLEPAIR_TYPE_HPP_
#define PARTICLEPAIR_TYPE_HPP_

#include "generic.hpp"
#include "cavity_grid.hpp"
#include "particle_type.hpp"
#include "particlepairtype_structs.hpp"
#include "splines.hpp"

/**
 * @brief Class for storing and accessing the parameters of each particle pair type
 *
 * Use defined type CParticlePairType_t for class declarations
 */

class CParticlePairType{

	/*Integer variable to store the number of particle pair types*/
	int nparticle_pairs;

	/*Vector of struct_particlepairtype data structures to store
	 * parameters and data of each particle pair type*/
	std::vector<struct_particlepairtype> particlepair_types;

	/*Mapping variable to access the vector element index of each particle
	 *  pair type by its label */
	std::tr1::unordered_map<std::string, int> label2index;

	/*Mapping variable to access the label of each particle pair
	 * type by its vector element index */
	std::tr1::unordered_map<int, std::string> index2label;

	/*Mapping variable to access the hard sphere potential cutoff distance of each
	 * particle pair type by its label */
	std::tr1::unordered_map<std::string, double> hs_cutoff;

	/*Mapping variable to access the Lennard-Jones potential cutoff distance of each
	 * particle pair type by its label */
	std::tr1::unordered_map<std::string, double> lj_cutoff;

	/*Mapping variable to access the monoatomic water potential cutoff distance of each
	 * particle pair type by its label */
	std::tr1::unordered_map<std::string, double> mw_cutoff;

	/*Mapping variable to access the offset values for determining the  cell indices
	 * that fall within the hard sphere cut off distance of each particle pair type,by its label */
	std::tr1::unordered_map<std::string, struct_cellindex_offset> hs_cellindex_neighbor_offset;

	/*Mapping variable to access the offset values for determining the cell indices
	 * that fall within the square well potential cut off distance of each particle pair type,by its label */
	std::tr1::unordered_map<std::string, struct_cellindex_offset> sqrwell_cellindex_neighbor_offset;

	/*Mapping variable to access the offset values tfor determining the cell indices
	 * that fall within the monoatomic water potential cut off distance of the solvent pair,by its label */
	std::tr1::unordered_map<std::string, struct_cellindex_offset> mwater_cellindex_neighbor_offset;

public:

	// constructor
	CParticlePairType(int n){
		nparticle_pairs = n * (n + 1) / 2;
		particlepair_types.resize(nparticle_pairs);
	}

	// destructor
	~CParticlePairType(){}

	// read methods
	/**
	 * @brief Function to reads PMF data for each particle pair type
	 *
	 * @return void
	 */
	void readPMFLookupFile(const CInputParameters_t &parameters);

	/**
	 * @brief Function to read and set the types of interactions specified by the user
	 * for each particle pair type
	 *
	 * @return void
	 */
	void readAndSetInteractionTypes(const char* fname);

	/**
	 * @brief Function to read and set the hard sphere potential cutoff distance between the two particle
	 * 		types of each particle pair type
	 *
	 * @return void
	 */
	void readAndSetHardSphereCutoffs(const char* fname);

	/**
	 * @brief  Function to read and set the Lennard-Jones potential cutoff distance
	 * between the two particle types of each particle pair
	 *
	 * @return void
	 */
	void readAndSetLJCutoffs(const char* fname);

	/**
	 * @brief Function to set the Lennard-Jones potential cutoff distance between the
	 * two particle types of each particle pair type
	 *
	 * @return void
	 */
	void setLJCutoffs(double lj_cutoff);

	/**
	 * @brief Function to set the id and labels of each particle pair type
	 *
	 * @return void
	 */

	void setTypeAndLabels(int nparticle_types,const CParticleType_t &particle_types);

	/**
	 * @brief Function to set the Lennard Jones potential parameters for each
	 * particle pair type
	 *
	 * @return void
	 */
	void setLennardJonesParameters(const CParticleType_t &particle_types);

	/**
	 * @brief Function to reads and set the square well potential parameters for
	 * each particle pair type
	 *
	 * @return void
	 */
	void readAndSetSquareWellParameters(const char* fname);


	//void setHardSphereCutOff(std::string label112,double hs_dist);
	/**
	 * @brief Function to set the offset values for determining the cavity grid cell indices
	 * that fall within the hard sphere cut off distance of each particle pair type
	 *
	 * @return void
	 */
	void setHardSphereCellIndexOffsetMatrix(const CCavityGrid_t &cavity_grid);

	/**
	 * @brief Function to sets the offset values for determining the cavity grid cell indices
	 * that fall within the square well potential cut off distance of each particle pair type
	 *
	 * @return void
	 */
	void setSqrWellWidthCellIndexOffsetMatrix(const CCavityGrid_t &cavity_grid);

	/**
	 * @brief Function to sets the monoatomic water potential cutoff distance between
	 * two solvent (water) hard spheres
	 *
	 * @return void
	 */
	void setMWaterPotentialCutoff(const double rc_cutoff);

	/**
	 * @brief Function to set the offset values for determining the cavity grid cell indices
	 * that fall within the monoatomic water potential cut off distance of solvent (water) pair type
	 *
	 * @return void
	 */
	void setMWaterCellIndexOffsetMatrix(const CCavityGrid_t &cavity_grid);

	/**
	 * @brief Function to change the hard sphere cutoff distance for each particle pair type
	 *
	 * @param[in]	pairlabel(string) Label of each particle pair type
	 * @param[in]	cutoff(double)	Hard sphere cutoff distance
	 *
	 * @return void
	 */
	void changeHardSphereCutoff(std::string pairlabel,double cutoff){
		hs_cutoff[pairlabel]=cutoff;
	}

	// get methods

	//std::vector<int> getCavityCountIndices(int i) const {return particlepair_types[i].cavitycount_indices;}

	/**
	 * @brief Function to get the number of particle pairs
	 *
	 * @return Number (int) of particle pair types
	 */
	int getNumPairs() const {return nparticle_pairs;}

	/**
	 * @brief Function to get the the vector element index of a particle pair type, by its label
	 *
	 * @param[in]	label(string)	Label of the particle pair type
	 *
	 * @return 	Index value (std::tr1::unordered_map<std::string, int>) of particle pair type
	 */

	int getIndexFromLabel(const std::string &label) const{	return label2index.find(label)->second;}

	/**
	 * @brief Function to get the label of a particle pair type, by its vector element index
	 *
	 * @param[in]	i(int)	vector element index of the particle pair type
	 *
	 * @return 	Label (std::tr1::unordered_map<int, std::string>)of particle pair type
	 */

	std::string getLabelFromIndex(int i)const {	return index2label.find(i)->second;}

	/**
	 * @brief Function to get the data structure of the particle pair type, by its vector element index
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return  Particle pair type data struct (struct_particlepairtype)
	 */
	struct_particlepairtype getType(int i) const{return particlepair_types[i];}

	/**
	 * @brief Function to get the hard sphere cutoff distance of each particle pair type, by its label
	 *
	 * @param[in]	pair_label(string)	Label  of the particle pair type
	 *
	 * @return 	Hard sphere cutoff distance(std::tr1::unordered_map<std::string, double>) in Angstrom
	 */

	double getHardSphereCutoff(std::string pair_label) const{return hs_cutoff.find(pair_label)->second;}

	/**
	 * @brief Function to get the Lennard-Jones potential cutoff distance of each particle pair type, by its label
	 *
	 * @param[in]	pair_label(string)	Label of the particle pair type
	 *
	 * @return 	Hard sphere cutoff distance (std::tr1::unordered_map<std::string, double>)
	 */
	double getLJCutoff(std::string pair_label) const{return lj_cutoff.find(pair_label)->second;}

	/**
	 * @brief Function to get the offset values for determining the  cell indices
	 * that fall within the hard sphere cut off distance of each particle pair type,by its label
	 *
	 * @param[in]	label (string)	Label of each particle pair type
	 *
	 * @return 	Data structure (std::tr1::unordered_map<std::string, struct_cellindex_offset>) of cell index offset values
	 */

	struct_cellindex_offset getHardSphereCellIndexOffsetMatrix(const std::string &label)
	const {return hs_cellindex_neighbor_offset.find(label)->second;}


	/**
	 * @brief Function to get the offset values for determining the  cell indices
	 * that fall within the square well potential cut off distance of each particle pair type,by its label
	 *
	 * @param[in]	label (string)	Label of each particle pair type
	 *
	 * @return Data structure (std::tr1::unordered_map<std::string, struct_cellindex_offset>) of cell index offset values
	 */

	struct_cellindex_offset getSqrWellWidthCellIndexOffsetMatrix(const std::string label)
	const{return sqrwell_cellindex_neighbor_offset.find(label)->second;}

	/**
	 * @brief Function to get the offset values for determining the  cell indices
	 * that fall within the Monoatomic water potential cut off distance of solvent pair type,by its label
	 *
	 * @param[in]	label (string)	Label of solvent pair type
	 *
	 * @return 	Data structure (std::tr1::unordered_map<std::string, struct_cellindex_offset>) of cell index offset values
	 */

	struct_cellindex_offset getMWaterCellIndexOffsetMatrix(const std::string label)
	const{return mwater_cellindex_neighbor_offset.find(label)->second;}

	/**
	 * @brief Function to determine whether the two particle types of a
	 * 			particle pair type interact via Coulomb potential
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return True or False (bool)
	 */

	bool getCoulombPotential(int i)const {return particlepair_types[i].coulomb_potential;}
	/**
	 * @brief Function to determine whether the two particle types of a
	 * 			particle pair type interact via Lennard Jones potential
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return True or False (bool)
	 */
	bool getLennardJonesPotential(int i)const {return particlepair_types[i].lennard_jones_potential;}
	/**
	 * @brief Function to determine whether the two particle types of a
	 * 			particle pair type interact via a PMF specified as a look up table
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return 	True or False (bool)
	 */
	bool getPMFLookTable(int i) const {return particlepair_types[i].pmf_lookup_table;}

	/**
	 * @brief Function to determine whether the two particle types of a
	 * 			particle pair type interact via a square well potential
	 *
	 * @param[in]	i(int)	Vector element index (integer) of the particle pair type
	 *
	 * @return 	True or False (bool)
	 */
	bool getSqrWellPotential(int i) const {return particlepair_types[i].square_well_potential;}

	/**
	 * @brief Function to get the Lennard-Jones potential potential well depth of each
	 * 			particle pair type
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return 	Well depth (double) of the Lennard-Jones potential in kcal/mol
	 */
	double getLJeps(int i) const {return particlepair_types[i].ljeps;}
	/**
	 * @brief Function to get the Lennard-Jones potential collision diameter of each
	 * 			particle pair type
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return 	L-J collision diameter (double) in Angstrom units
	 */

	double getLJsig (int i) const {return particlepair_types[i].ljsig;}
	/**
	 * @brief Function to get the numerator of the long-range term in the Lennard-Jones potential
	 *
	 * Aij = 4.0*ljeps*ljsig^12
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return 	Numerator (double) of the long-range term in the Lennard-Jones potential
	 */
	double getAij(int i) const {return particlepair_types[i].Aij;}
	/**
	 * @brief Function to get the numerator of the short-range term in the Lennard-Jones potential
	 *
	 * Bij = 4.0*ljeps*ljsig^6
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return 	Numerator (double) of the short-range term in the Lennard-Jones potential
	 */
	double getBij (int i) const {return particlepair_types[i].Bij;}

	/**
	 * @brief Function to get the radius of the first particle type of each particle pair type
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return	Radius (double) of the first particle type in the particle pair type, in Angstom units
	 *
	 */
	double getRadius1(int i) const {return particlepair_types[i].radius1;}
	/**
	 * @brief Function to get the radius of the second particle type of each particle pair type
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return 	Radius (double) of the second particle type in the particle pair type, in Angstom units
	 *
	 */
	double getRadius2(int i) const {return particlepair_types[i].radius2;}

	/**
	 * @brief Function to get the width to contact ratio of the square well potential between
	 * 		the particle types of each particle pair type
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return 	Square well potential width	to contact ratio (double)	 *
	 */

	double getSqrwell_wdth2contact_ratio(int i) const {return particlepair_types[i].sqrwell_wdth2contact_ratio;}

	/**
	 * @brief Function to get the square well potential well depth of each particle pair type
	 *
	 * @param[in]	i(int) Vector element index of the particle pair type
	 *
	 * @return Square well potential well depth (double)	 *
	 */

	double getSqrWellDepth (int i) const {return particlepair_types[i].sqrwell_depth;}
	/**
	 * @brief Function to get the ID of the first particle type in the particle pair type
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return 	ID (integer) of the first particle type in the particle pair type
	 */
	int getTypeIndex1 (int i) const {return particlepair_types[i].ptype1;}
	/**
	 * @brief Function to get the ID of the second particle type in the particle pair type
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 *
	 * @return	ID (integer) of the second particle type in the particle pair type
	 */
	int getTypeIndex2 (int i) const {return particlepair_types[i].ptype2;}

	int getNLookup(int i) const {return particlepair_types[i].n_lookup;}

	/**
	 * @brief Function to select a PMF value from the PMF look up table
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 * @param[in]	j(int)	Vector element index of the PMF look up table
	 *
	 * @return 	PMF value (double) in kcal/mol/e
	 */

	double getPMFLookup(int i,int j) const {return particlepair_types[i].pmf_lookup[j];}

	/**
	 * @brief Function to get the radial distance value from the
	 * 		PMF look up table of a particle pair type
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 * @param[in]	j(int)	Vector element index of the PMF look up table
	 *
	 * @return r_lookup	Radial distance (double) in Angstrom
	 */
	double getRLookup(int i,int j) const {return particlepair_types[i].r_lookup[j];}

	/**
	 * @brief Function to get the coefficients for cubic and bicubic spline interpolation
	 *
	 * @param[in]	i(int)	Vector element index of the particle pair type
	 * @param[in]	j(int)	Vector element index of the PMF look up table
	 *
	 * @return 	Coefficients (struct_spline_coeff) for spline interpolation
	 */
	struct_spline_coeff getSplineCoeff(int i,int j) const {return particlepair_types[i].spline_coeff[j];}
	// write methods

	void writePMFLookupTableToFile();
};
typedef class CParticlePairType CParticlePairType_t;


#endif /* PARTICLEPAIR_TYPE_HPP_ */
