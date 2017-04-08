/*
 * particlepairtype_structs.hpp
 *
 *  Created on: Aug 7, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file particlepairtype_structs.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains declarations of data structures for particle pair types
 */

#ifndef PARTICLEPAIRTYPE_STRUCTS_H_
#define PARTICLEPAIRTYPE_STRUCTS_H_


#include "splines.hpp"

/**
 * @brief Data structure for each particle pair type
 */

struct particlepairtype_struct{

	/*Integer variable to store the vector/array element index of
	 * each particle pair, starting at 0 */
	int index;

	/*Integer variable to store the id of the first particle type in the particle pair.
	 * The ID is 1 more than the array index of the particle type in the patrticle type
	 * struct vector  */
	int ptype1;

	/*Integer variable to store the id of the second particle type in the particle pair.
	 * The ID is 1 more than the array index of the particle type in the patrticle type
	 * struct vector  */
	int ptype2;

	int n_lookup;

	/*A boolean variable to indicate whether the two particle types in the particle pair type
	 * interact via Coulomb potential  */
	bool coulomb_potential;

	/*A boolean variable to indicate whether the two particle types in the particle pair type
	 * interact via hard sphere repulsion potential  */
	bool hard_sphere_repulsion;

	/*A boolean variable to indicate whether the two particle types in the particle pair type
	 * interact via Lennard-Jones potential  */
	bool lennard_jones_potential;

	/*A boolean variable to indicate whether the two particle types in the particle pair type
	 * interact via attractive square well potential  */
	bool square_well_potential;

	/*A boolean variable to indicate whether the two particle types in the particle pair type
	 * interact via a PMF that is provided by the user as a lookup table*/
	bool pmf_lookup_table;

	/*Radius of the first particle type in the particle pair type */
	double radius1;
	/*Radius of the second particle type in the particle pair type */
	double radius2;

	/*LJ potential well depth of the particle pair type (kcal/mol) */
	double ljeps;
	/*LJ collision diameter of the particle pair type (Angstrom) */
	double ljsig;

	double Aij;  // Aij = 4.0*ljeps*ljsig^12
	double Bij; // Bij = 4.0*ljeps*ljsig^6

	/*Well depth of the square well potential for the particle pair type*/
	double sqrwell_depth;

	/* Width to contact ratio of square well potential */
	double sqrwell_wdth2contact_ratio;

	/*Label of the first particle type in the particle pair type */
	std::string label1;
	/*Label of the second particle type in the particle pair type */
	std::string label2;

	/*Label of the particle pair type, created by concatenating the labels of the
	 * first and second particle types */
	std::string label12;  // label1_label2

	/*Label of the particle pair type, created by concatenating the labels of the
	 * second and first particle types */
	std::string label21; // label2_label1;

	//	std::vector<int> cavitycount_indices;
	/*Vector of floating point values to store the radial distance values read from
	 * PMF look up table file*/
	std::vector<double> r_lookup;
	/*Vector of floating point values to store the PMF values read from
	 * PMF look up table file*/
	std::vector<double> pmf_lookup;
	/*Vector of struct_spline_coeff containing the coefficients for cubic and bicubic
	 * spline interpolation*/
	std::vector<struct_spline_coeff> spline_coeff;
};

typedef struct particlepairtype_struct struct_particlepairtype;


#endif /* PARTICLEPAIRTYPE_STRUCTS_H_ */
