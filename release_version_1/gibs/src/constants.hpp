/*
 * constants.hpp
 *
 *  Created on: Apr 13, 2015
 *      Author:  Dennis G. Thomas
 *
 *
 *  @brief Defines the constants used in GIBS
 *
 *  @file constants.hpp
 *  @author Dennis G. Thomas
 */

#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_


/* define constants */

#define MAX_NIONTYPE 4		// Maximum number of ion types allowed in a simulation
#define MAX_NIONPAIR 10		// n*(n+1)/2, where n denotes the max_niontype
#define PI  3.1415926536	// value of pi
#define VACDIEL 1.0		// Vacuum dielectric constant
#define KB	1.3806620*pow(10,-23)		// Boltzmann's constant in J/K
#define FUNDCHG 1.602176565*pow(10,-19)		// Fundamental charge in Coulomb unit
#define VACPERMIT 8.854187*pow(10,-12)  // Vacuum permittivity in F/m or C^2/(N m^2)
#define AVOGNUM 6.0220450*pow(10,23) // Avogadro's number
#define CONVJTOCAL 0.239005736		// Conversion factor to convert joules to calories
#define NION_MIN  1	// Minimum number of particles of a given ion type required in the box
#define MAXBINNUM 300 // Maximum number of bins allowed in the RDF calculation.

#define MAX_NUM_SOLUTE_ATOMS 10000 // Maximum allowed number of atoms in a solute
#define RDF_SCALING_FACTOR 0.0001

#endif /* CONSTANTS_HPP_ */
