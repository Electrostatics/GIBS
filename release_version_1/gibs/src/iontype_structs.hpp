/*
 * iontype_structs.hpp
 *
 *  Created on: Aug 7, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file iontype_structs.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains data structure definitions for ion types
 */

#ifndef IONTYPE_STRUCTS_HPP_
#define IONTYPE_STRUCTS_HPP_

/**
 * @brief Data structure to store ion type parameters and thermodynamic properties
 *
 */
struct iontype_struct {

	/*
	 * Floating point variable to store the charge in (e units) of the
	 * ion type (e.g, 1 for Na+, -1 for Cl-, 2 for Mg2+)
	 */
	double charge;
	/*
	 * Floating point variable to store the
	 *  radius of the ion type, in Angstrom unit
	 */
	double radius;

	/*
	 * Floating point variable to record the bulk concentration of
	 * the ion type in molar (M) units
	 */
	double conc;
	/*
	 * Integer variable to record the number of ions present in
	 * the simulation box during simulation
	 */
	long int num;

	/* Integer variable to assign an identification number (starting from 1)
	 *  to each ion type in a simulation in consecutive order
	 */
	int itype;

	/*
	 * Floating point variable to record the excess chemical
	 * potential of the ion species (kcal/mol)
	 */
	double mu_ex;
	/*
	 * Floating point variable to record the chemical potential of the ion species,
	 * calculated in the nth iteration of a GCMC simulation
	 */
	double mu_old;
	/*
	 * Floating point variable to record the chemical potential of the ion species,
	 * calculated in the (n+1)th iteration of a GCMC simulation
	 */
	double mu_new;
	/*
	 * Floating point variable to record the average number of the ion species
	 * at equilibrium
	 */
	double nion_avg;
	/*
	 * Floating point variable to record the average concentration (M)
	 * of the ion species at equilibrium
	 *
	 */
	double cion_avg;
	/*
	 * Floating point variable to record the average of the
	 * square of the number of ion species at equilibrium
	 *
	 */
	double nion_sqr_avg; // square of the average concentration
	/*
	 * Floating point variable to record the standard deviation
	 * of the number of ion species at equilibrium.
	 */
	double nion_sd;
	/*
	 * Floating point variable to store the bin size used for calculating
	 * the ion pair radial distribution \
	 * function (RDF) - same value is used for all ion-ion pair RDFs.
	 */
	double delr_rdf;
	/*
	 * Integer variable to store the number of bins over which
	 * the RDFs are computed.
	 */
	int nrdfbin;	// number of bins over which the RDFs are computed.
	/*
	 * Floating point variable to store the well depth of the Lennard-Jones
	 * potential (kcal/mol).
	 */
	double lj_epsilon;
	/*
	 * Floating point variable to store the collision diameter (Angstrom) for
	 * pairs of identical ions, at which the Lennard-Jones potential is zero.
	 */
	double lj_sigma;
	/*
	 * Ion-accessibility map computed based on the parameters given in the
	 * solute's charge/coordinate file.
	 * It is computed only in GCMC_WITH_SOLUTE_ATOM_MODEL simulation.
	 */
	double ***accessibility;

	/*
	 * String variable to store the label used to identify the ion type (e.g., Na, Cl, Mg, etc.)
	 */
	std::string label;
	/*
	 *  Integer variable to store the number of neighboring cells
	 *  occupied by the ion along x dimension
	 */
	int hs_ncx;
	/*
	 *  Integer variable to store the number of neighboring cells
	 *  occupied by the ion along y dimension
	 */
	int hs_ncy;
	/*
	 *  Integer variable to store the number of neighboring cells
	 *  occupied by the ion along z dimension
	 */
	int hs_ncz;
};
typedef iontype_struct struct_iontype;


#endif /* IONTYPE_STRUCTS_HPP_ */
