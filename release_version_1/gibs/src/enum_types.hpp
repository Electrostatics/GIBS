/**
 * enum_types.hpp
 *
 *  Created on: Apr 14, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file enum_types.hpp
 *  @author Dennis G. Thomas
 *
 *
 *  @brief Contains the enumerated data types
 */

#ifndef ENUM_TYPES_HPP_
#define ENUM_TYPES_HPP_

/* enum declarations */

/**
 * @brief Enumerates the different types of solute models
 *
 * Use defined type enum_solute_models for enum declarations
 */

enum solute_models_enum{
	none = 0,
	all_atom = 1,
	cylindrical_polyion = 2,
	spherical_macroion = 3
};

typedef enum solute_models_enum enum_solute_models;

/**
 * @brief Enumerates the different types of GCMC simulation
 *
 * Use defined type enum_simulation_type for enum declarations
 */
enum simulation_type_enum {

	GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION = 0, // compute excess chemical potential and mean activity coefficients of individual ions and salt.
	GCMC_WITH_SOLUTE_ALLATOM_MODEL = 1, // GCMC simulation for computing ion distributions around a fixed solute represented in atomic detail.
	GCMC_WITH_SPHERICAL_MACROION = 2, // GCMC simulation for computing ion distributions around a fixed spherical macroion.
	GCMC_WITH_CYLINDRICAL_POLYION = 3 // GCMC simulation for computing ion distributions around a fixed cylindrical polyion of finite length.

};

typedef enum simulation_type_enum enum_simulation_type;


/**
 * @brief Enumerates the different ways of changing the number of particle species in the system
 *
 * Use defined type enum_state_change for enum declarations
 */

enum state_change_enum{
	insertion = 0,
	deletion = 1

};
typedef enum state_change_enum enum_state_change;

/**
 * @brief Enumerates the axis directions perpendicular to which the cylindrical radial distribution
 * 	  functions have to be computed.
 *
 * Use defined type	enum_cylindrical_rdf_axis_direction for enum declarations
 *
 * 	- Used only in GCMC_SOLUTE_ATOM_MODEL simulations
 */
enum cylindrical_rdf_axis_direction_enum {

	X_AXIS = 0, // Use this if the longest principal axis of the solute molecule aligns with the x-axis.
	Y_AXIS = 1, // Use this if the longest principal axis of the molecule aligns with the y-axis.
	Z_AXIS = 2 // Use this if the longest principal axis of the molecule aligns with the z-axis.
};

typedef enum cylindrical_rdf_axis_direction_enum enum_cylindrical_rdf_axis_direction;


/**
 * @brief Enumerates the reference position from which the radial distribution function is computed.
 *
 * 	Use defined type enum_rdf_ref for enum declarations
 *
 * 	- Not used in GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION
 */

enum rdf_ref_enum {
	/*Option for calculating the cylindrical RDF where the reference
	 * line is the axis of a cylinder.*/

	AXIS_OF_CYLINDER = 0,
	/*Option for calculating the spherical RDF where the reference point is the center
	 * of a sphere (atom or spherical macroion).*/
	CENTER_OF_SPHERE = 1
};

typedef rdf_ref_enum enum_rdf_ref;

/**
 * @brief Enumerates the format of solute molecule coordinate files.
 *
 * Use enum_molecule_coordinate_file_type for enum declarations
 *
 *	- Used for entering the type of solute's charge/atom-coordinate file
 * 	- Used as value for the keyword SOLUTE_ALLATOM_COORDINATE_FILE_TYPE in the input parameter file.
 */

enum molecule_coordinate_file_type_enum {
	XYZR = 0, PQR = 1

};
typedef enum molecule_coordinate_file_type_enum enum_molecule_coordinate_file_type;
#endif /* ENUM_TYPES_HPP_ */
