/**
 * input_parameters.hpp
 *
 *  Created on: Apr 13, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file input_parameters.hpp
 *  @author Dennis G. Thomas
 */

#ifndef INPUT_PARAMETERS_HPP_
#define INPUT_PARAMETERS_HPP_

#include "generic.hpp"
#include "enum_types.hpp"
#include "constants.hpp"
#include "iontype_structs.hpp"

/**
 * @brief Class for storing the simulation set up parameters specified by the user
 *
 * Use defined type CInputParameters_t for class declarations
 */
class CInputParameters{

public:
	/* Type of simulation to run. */
	enum_simulation_type simulation_type;

	/*Enumerated variable to store the type of solute model (allowed values: none,all_atom) */
	enum_solute_models solute_model;
	/*Floating point variable to store the cavity grid spacing in Angstrom */
	double cavity_grid_spacing;
	/*Floating point variable to store the maximum allowed number of particles
	 * (value is used for pre-allocating the size of arrays and vectors that will
	 * store particle information)*/
	int maximum_particle_number;
	/* Integer variable to store a non-zero positive integer value assigned
	 * internally (in gcmc_main.cpp file) for each simulation type
		 based on the value of simulation_type.*/
	int run_type;

	/* Integer variable to store the total number of iterations for updating excess
	 * chemical potential. niter = 1 for all simulation types (simulation_type)
		 except for GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION. If PID method is used
		 for updating the chemical potential, then it is sufficient to set niter to 1*/
	int niter;
	/* Integer variable to store the total number of GCMC steps in the simulation. */
	int nsteps;
	/*Integer variable to store the number of insertion/deletion cycles in a
	 * GCMC simulation step. If the simulation_type is GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION,
	 * its value is set to a small number (e.g., 3). For the other simulation types,
	 * its value is typically set to the average expected number of total particles in the system */
	int gcmc_cycl;
	/* Integer variable to store the number of single-particle displacement cycle/steps
	 * in a GCMC simulation step. If the simulation_type is GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION,
	 * its value is set to a small number, typically greater than gcmc_cycl. For the other simulation
	 * types, its value is set such that at least 60-70% of simulations steps are single-particle
	 * displacements. */
	int disp_cycl;
	/* Integer variable to store the number of GCMC steps for equilibration */
	int neqlb_steps;
	/* Integer variable to store the number of ion types in the simulation */
	int nion_type;
	/* Floating point variable to store the length of the simulation box along x-axis in Angstrom units. */
	double xlen;
	/* Floating point variable to store the length of the simulation box along y-axis in Angstrom units.*/
	double ylen;
	/* Floating point variable to store the length of the simulation box along z-axis in Angstrom units. */
	double zlen;
	/* Floating point variable to store the simulation box grid spacing along x-axis in Angstrom units. */
	double hx;
	/* Floating point variable to store the simulation box grid spacing along y-axis in Angstrom units. */
	double hy;
	/* Floating point variable to store the simulation box grid spacing along z-axis in Angstrom units. */
	double hz;

	/* Floating point variable to store the lowest x-axis value of a right-handed
	 * coordinate system,where the leftmost bottom corner of the simulation box is located. */
	double xmin;
	/* Floating point variable to store the lowest y-axis value of a right-handed
	 * coordinate system,where the leftmost bottom corner of the simulation box is located. */
	double ymin;
	/* Floating point variable to store the lowest z-axis value of a right-handed
	 * coordinate system,where the leftmost bottom corner of the simulation box is located. */
	double zmin;

	/* Array of ion type data structure (struct_iontype). The array dimension is determined by the number
	 * of ion types specified by the user in the input parameter file*/
	std::vector<struct_iontype> iontypes;

	/* Floating point variable to store the temperature in Kelvin (K) units*/
	double temperature;
	/* Floating point variable to store the solvent dielectric constant or relative permittivity */
	double solvent_dielectric;
	/* Floating point variable to store the solute dielectric constant or relative permittivity.
	 * Not used so far in any simulations*/
	double solute_dielectric;
	/*Integer variable to store the value (0 or 1) for selecting the
	 * Sloth-Sorensen correction term in GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION simulation,
	 * with the adaptive GCMC algorithm (use_pid=FALSE).
	 * Use 1 to apply the Sloth-Sorensen Correction Term; use 0 otherwise. */
	int icorr;

	/* Integer variable to store the number of ion pairs (computed internally
	 * based on the value of nion_type) */

	int nionpair;

	/*-------Parameters for the spherical macro-ion; used only for the
	 * GCMC_WITH_SPHERICAL_MACROION (run_type = 3) simulation.The center of the
	 * spherical macro-ion coincides with the center of the simulation box (0,0,0).
	 -------*/

	/* Floating point variable to store the charge of the spherical macro-ion (e units),
	 * located at its center (0,0,0).*/
	double sphere_charge;
	/* Floating point variable to store the radius of the spherical macro-ion (in Angstrom). */
	double sphere_radius;

	/* ------- Parameters for the cylindrical polyion; used only for the
	 * GCMC_WITH_CYLINDRICAL_POLYION (run_type = 4) simulation.The center of the cylinder
	 * coincides with the center of the simulation box.The cylindrical axis is vertically
	 * aligned parallel to the z-axis.------- */

	/* Floating point variable to store the length (in Angstrom units) of each monomer unit
	 * of the cylindrical polyion.*/
	double cylinder_monomer_length;

	/* Integer variable to store the number of monomers in the cylindrical polyion. */
	int cylinder_num_monomers;

	/* Floating point variable to store the charge (in e units) per monomer. */
	double cylinder_charge_per_monomer;

	/* Floating point variable to store the radius of the cylindrical polyion
	 * (in Angstrom units).*/
	double cylinder_radius;

	/* Integer variable to store the value (0 or 1) for determining whether the charges of
	 * the cylindrical polyion are discretely or continuously distributed along the cylindrical
	 * axis. Value = 0 if discrete charges of 1 e unit is placed at the center of each
	 * monomer. Value = 1 for continuous charge distribution, in which case the interaction
	 * energy is computed using the formula for the potential between a cylindrical polyion (continuous charge
	 * distribution) and an ion.
	 */

	int select_cylinder_polyion_potential;

	/* Floating point variable to store the length of the cylindrical polyion (in Angstrom). It is
	 * used only if select_cylinder_polyion_potential = 1.
	 */
	double cylinder_length;

	/* Integer variable to store the number of segments (monomers) over which the RDF of an
	 * ion from the cylindrical axis is computed.*/

	int rdf_cylinder_num_segments;

	/* Floating point variable to store the minimum spacing (Angstrom) between particles
	 * initially placed in the box, without the grid insertion algorithm*/

	double min_ptcl_spacing;

	/* Integer variable to store the value (0 or 1) for specifying whether to create the initial
	 * state configuration of the system or read the state from a file. Value is 0 if the
	 *  initial configuration of the particles in the system is randomly set by
	 *  the program. Value is 1 if the initial configuration of the particles is read from a file.
	 */

	int start_state;

	/* Integer variable to store the initial number of particles in the system.
	 * Its value is set through the keyword INITIAL_PTCL_NUMBER in the inputparameters.in
	 * file, if the state has to be created instead of reading the state from a file.*/

	int npart_init;

	/* Enumerated variable to store the format of solute's charge/atom-coordinate file*/

	enum_molecule_coordinate_file_type solute_allatom_coordinate_filetype;

	/* Enumerated variable to store the cylindrical axis direction, perpendicular
	 * to which the cylindrical radial distribution functions have to be computed.
	 */

	enum_cylindrical_rdf_axis_direction cylindrical_rdf_axis_direction;

	/* Boolean variable to specify whether to calculate the solute-particle RDF*/

	bool calculate_solute_particle_rdf;
	/* Boolean variable to specify whether to calculate the particle pair RDFs */

	bool calculate_particle_pair_rdf;

	/* Enumerated variable to specify whether to compute the RDF from a reference line or
	 * from a reference point. It is used only in GCMC_WITH_SOLUTE_ATOM_MODEL simulation. */

	enum_rdf_ref rdf_ref;

	/* Floating point variable to store the bin size for calculating the radial distribution
	 * functions (RDFs). If the value is zero, it is computed internally based on the number
	 * of bins (num_bins).
	 */

	double delr_rdf;

	/* Integer variable to store the number of bins over which the RDFs are computed.
	 * If a non-zero positive bin size value (delr_rdf) is read from the input file,
	 * then the value of nrdfbin is computed internally. */

	int nrdfbin;
	/* Floating point variable to store the maximum radial distance within which the
	 * RDFs are computed.  */

	double rdf_max_distance;
	/* Floating point variable to store the fraction of the solute's principal axis length
	 * from which the cylindrical RDF is computed. For example, if the solute is B-DNA,
	 * one can use a value1.0 to include the whole length of the B-DNA or a fraction of that
	 * length (e.g, 0.8). */

	double rdf_ref_axis_length_fraction;

	/*Floating point variable (double) to store the radial distance spacing (in Angstrom)
	 * at which the particle pair PMFs are listed in the look up table */
	double lookup_pmf_delr;

	/*Integer number to store the number of rows in the PMF lookup table */
	int lookup_pmf_numr;

	/* Boolean variable to specify whether to record the acceptance rates for GCMC moves */
	bool record_acceptance_rate;

	/* Boolean variable to specify whether to record count of each particle type in the
	 * simulation box after every n GCMC steps (record_every_n_steps)*/
	bool record_particletype_count;

	/* Boolean variable to specify whether to record energy of the system after
	 * every n GCMC steps (record_every_n_steps)
	 */
	bool record_system_energy;

	/* Integer variable to store the step period for recording the particle counts and energy */
	long unsigned int record_every_n_steps;

	/* Integer variable to store the number of states to record after equilibration. */

	int num_record_states;
	/* Floating point variable to store the radius of the solvent hard sphere */
	double solvent_radius;
	/* Floating point variable to store the solvent concentration in Molar (M) units */
	double solvent_conc;
	/* Floating point variable to store the excess chemical potential of the solvent in
	 * kcal/mol units */
	double solvent_mu_ex;
	/* Floating point variable to store the solvent packing fraction*/
	double solvent_packing_fraction;
	/* String variable to store the label for identifying the solvent particle type*/
	std::string solvent_label;

	/* String variable to store the name of the open DX formatted file comprising the
	 * solute electrostatic potential values at every grid point. The potential is read in
	 * units of kT/e but internally converted to kcal/(mol e).The potential is computed in
	 * the absence of ions using an external Poisson solver (APBS).*/

	std::string pot_dxmap_file;

	/* String variable to store the name of the solute's charge/atom-coordinate file */
	std::string solute_allatom_coordinate_filename;

	/* Boolean variable to specify whether to select the solvent primitive model, in which the
	 * solvent molecules are represented as hard spheres. */
	bool use_spm;

	/* Boolean variable to specify whether to include an attractive (short-range)
	 * interaction between solvent and ions*/
	bool solvent_ion_attraction;
	/* Integer variable to store the number of segments for grouping the cavity grid cells*/
	int num_cavitygrid_segments;

	/* Boolean variable to use the monoatomic water model (Molinero and Moore 2008)	 */
	bool use_mwater;

	/* Boolean variable to specify whether to update the chemical potential
	 * after every insertion/deletion attempt using PID controller method.
	 */

	bool use_pid;

	/* Floating point variable to store the Lennard-Jones potential cut-off distance
	 * (typically used value is 12 A) */
	double lj_cutoff;

	// read methods
	void readFromFile();

	// write methods
	void printOut();

};

typedef class CInputParameters CInputParameters_t;


#endif /* INPUT_PARAMETERS_HPP_ */
