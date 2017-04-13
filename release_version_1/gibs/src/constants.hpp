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
