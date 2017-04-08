/*
 * rdf_structs.hpp
 *
 *  Created on: Jul 9, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file rdf_structs.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Data type structures for Radial distribution functions
 *
 */

#ifndef RDF_STRUCTS_HPP_
#define RDF_STRUCTS_HPP_

#include "generic.hpp"

/**
 * @brief  Data type structure for storing radial distribution functions
 *
 * Use defined type struct_rdf_data for RDF struct declarations
 */
struct rdf_data_struct{
	/*Allocatable array of floating point values to store the cumulative sum of particle counts
	 * less than the maximum distance value of each bin*/

	double *ndist;

	/*Allocatable array of floating point values to store the actual
	 * count of particles in each bin */
	double *nbin; //

	// ndist and nbin values are updated after every insertion, deletion, or displacement of a particle

	/* Allocatable array of floating point values to store the 'ndist' values before a state change */
	double *ndist_prev;  //
	/*Allocatable array of floating point values to store the 'nbin' values before a state change*/
	double *nbin_prev;	//

	/*Allocatable array of floating point values to store the 'nbin' values averaged over all states*/
	double *nbin_avg;

	/*Allocatable array of floating point values to store the square of the 'nbin' value averaged
	 * over all states*/
	double *nbin_avg_sqr;

	/*Allocatable array of floating point values to store the radial distribution */
	double *grdf;

	/*Allocatable array of floating point values to store the 'ndist' values averaged over
	 *  all states */
	double *ndist_avg;

	/*Allocatable array of floating point values to store the square of the 'ndist' values
	 * averaged over all states */
	double *ndist_avg_sqr;
};

typedef struct rdf_data_struct struct_rdf_data;

/**
 * @brief Data type structure for storing the RDF data and parameters of each particle pair type
 *
 * Use defined type struct_particlepair_rdf for RDF struct declarations
 */

struct particlepair_rdf_struct{

	/* Data type structure for storing RDF data */
	struct_rdf_data *data;
	/* Integer variable to store the total number of state changes from which the
	 * statistics for RDF data were collected*/
	int ndist_icount;
	/* Integer variable to store the total count of the reference particle type*/
	int ref_ptcl_count;
	/* Integer variable to store the particle type ID of the reference particle type*/
	int ref_ptcl_type;

	/*Integer variable to store the particle type of which the RDF is computed in r
	 * reference to the 'ref_ptcl_type'*/
	int rdf_ptcl_type;

	/*Floating point variable to store the mid-point radial distance value of the first bin */
	double rmid_firstbin;

	/*String variable to store the particle pair label used for identifying the
	 * particle pair RDF*/
	std::string particlepair_label;

	/*Allocatable array of floating point values to store the radial distance
	 * values over which the RDF data was obtained */

	double *rdist;
};

typedef particlepair_rdf_struct struct_particlepair_rdf;


/**
 * @brief Data type structure to store the parameters for setting up the RDF data structures
 * 			of each particle pair type
 *
 * Use defined type struct_rdf_parameters for struct declarations
 */
struct rdf_parameters_struct{
	/*Integer variable to store the number of particle types */
	int nparticle_types;
	/*Integer variable to store the number of particle pair types */
	int nparticle_pairs;
	/*Integer variable to store the number of bins */
	int nrdfbin;
	/*Floating point variable to store the bin size */
	double delr_rdf;	// bin size
	/*Floating opint variable to store the length of the simulation
	 * box along x-axis (Angstrom)*/
	double xlen;
	/*Floating point variable to store the length of the simulation
	 * box along y-axis (Angstrom)*/
	double ylen;
	/*Floating point variable to store the length of the simulation
	 * box along z-axis (Angstrom)*/
	double zlen;
	/*Floating point variable to store the lowest value of x along the x-axis
	 * box dimension (Angstrom)*/
	double xmin;
	/*Floating point variable to store the lowest value of y along the z-axis
	 * box dimension (Angstrom)*/
	double ymin;
	/*Floating point variable to store the lowest value of z along the z-axis
	 * box dimension (Angstrom)*/
	double zmin;
	/*Floating point variable to store the volume of the simulation box (Angstrom cube)*/
	double box_vol;

};
typedef rdf_parameters_struct struct_rdf_parameters;

/**
 * @brief Data type structure to store the parameters and radial distribution function
 * 		of each particle type  from the solute
 *
 * Use defined type struct_soluteallatom_particle_rdf for struct declarations
 */

struct soluteallatom_particle_rdf_struct{

	/* Data type structure for storing RDF data */
	struct_rdf_data *data;
	/*Integer variable to store the index of particle type */
	int particletype_index;
	/*Integer variable to store the number of bins */
	int nrdfbin;
	/*Integer variable to store the bin size (Angstrom) */
	double delr_rdf;
	/*Floating point variable to store the length of the axis used as reference for
	 * computing the RDF  */
	double ref_axis_length;
	/*Floating point variable to store the x-coordinate value at the center of the reference axis
	 * from which the RDF of particle type is computed */
	double ref_xcoord;
	/*Floating point variable to store the y-coordinate value at the center of the reference axis
	 * from which the RDF of particle type is computed */
	double ref_ycoord;
	/*Floating point variable to store the z-coordinate value at the center of the reference axis
	 * from which the RDF of particle type is computed */
	double ref_zcoord;

	/*Floating point variable to store the value 0 if the cylindrical RDF reference axis
	 * is the x-axis, and 1 if the reference axis is y-axis or z-axis*/
	double coeff_x;
	/*Floating point variable to store the value 0 if the cylindrical RDF reference axis
	 * is the y-axis, and 1 if the reference axis is x-axis or z-axis*/
	double coeff_y;
	/*Floating point variable to store the value 0 if the cylindrical RDF reference axis
	 * is the z-axis, and 1 if the reference axis is x-axis or y-axis*/
	double coeff_z;

	/*Floating point variable to store the minimum value of x
	 * among the (x,y,z) atom coordinates*/

	double solute_x_min;
	/*Floating point variable to store the minimum value of y
		 * among the (x,y,z) atom coordinates*/
	double solute_y_min;
	/*Floating point variable to store the minimum value of z
		 * among the (x,y,z) atom coordinates*/
	double solute_z_min;

	/*Radius of the reference atom of solute (or of a macroion) from which
	 * the RDFs are computed  */
	double atom_sphere_radius;

	/* Integer variable to store the total number of state changes from which the
		 * statistics for RDF data were collected*/
	int ndist_icount;

	double *rdist;
	/*Enumeration type variable to store the cylindrical axis used as reference for
	 * computing the RDF*/
	enum_cylindrical_rdf_axis_direction  cylindrical_rdf_axis_dir;	// -  X-, Y-, or Z- axis

	/*String variable to store the particle type label used for identifying the particle type */
	std::string particletype_label;
};
typedef soluteallatom_particle_rdf_struct struct_soluteallatom_particle_rdf;

#endif /* RDF_STRUCTS_HPP_ */
