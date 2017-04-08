/**
 * solvent_type.cpp
 *
 *  Created on: Aug 5, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file solvent_type.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Definitions for member functions of class CSolvent
 *
 */

#include "solvent_type.hpp"

/*
 * @brief CSolvent constructor
 *
 * Sets default values for the solvent parameters
 *
 * @param[in] box_volume(double) Box volume in Angstrom cube units
 *
 */
CSolvent::CSolvent(double box_vol){
	charge = 0.0; //
	radius = 1.4; // radius in Angstrom
	label = "Water";
	packing_fraction = 0.001;
	dielectric = 78.5;
	/*target_num = (int) (box_vol*packing_fraction*3/
			(PI*4*radius*radius*radius));*/
	target_num = (int) (box_vol*packing_fraction/
				(pow(2.0*radius,3)));
	conc = numToConcMolarity(target_num,box_vol);
	mu_ex = 0.0;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Sets the solvent parameters with values specified by the user in
 * 			in the input parameter file
 *
 * Total charge is always set to zero
 *
 * @param[in]	parameters(CInputParameters_t)	Input Parameters specified by user
 * @param[in]	box_vol(double)		Box volume in Angstrom cube
 *
 * @return void
 */
void CSolvent::setParameters(const CInputParameters_t &parameters,double box_vol){

	radius = parameters.solvent_radius;
	label = parameters.solvent_label;
	packing_fraction = parameters.solvent_packing_fraction;
	dielectric = parameters.solvent_dielectric;
	charge = 0;
	/*target_num = (int) (box_vol*packing_fraction*6/
				(PI*8*pow(radius,3)));*/
	target_num = (int) (box_vol*packing_fraction/
				(pow(2.0*radius,3)));
	conc = numToConcMolarity(target_num,box_vol);
	mu_ex = parameters.solvent_mu_ex;

}
