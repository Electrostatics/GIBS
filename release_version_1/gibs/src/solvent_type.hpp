/**
 * solvent_type.hpp
 *
 *  Created on: Apr 13, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file solvent_type.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Class declarations and function prototypes/definitions
 *  		for solvent type
 */

#ifndef SOLVENT_TYPE_HPP_
#define SOLVENT_TYPE_HPP_

#include "generic.hpp"
#include "constants.hpp"
#include "input_parameters.hpp"
#include "convert_functions.hpp"

/**
 * @brief Class declaration for storing solvent parameters
 *
 */
class CSolvent{

	/* Floating point variable to store the total charge (in e units) on
	 * solvent molecule, if represented as hard spheres
	 */
	double charge;
	/* Floating point variable to store the radius (in Angstrom) of
	 * solvent hard sphere
	 */
	double radius;
	/*
	 * String variable to store the label used for
	 * identifying a particle as solvent
	 */
	std::string label;
	/*
	 * Floating point variable to store the solvent packing fraction
	 */
	double packing_fraction;
	/*
	 * Floating point variable to store the solvent dielectric constant
	 */
	double dielectric;
	/*
	 * Integer variable to store the target number of solvent hard
	 * spheres in SPM simulations
	 */
	int target_num;
	/*
	 * Floating point variable to store the solvent concentration (M)
	 */
	double conc;
	/*
	 * Floating point variable to store the excess chemical potential
	 * of solvent hard spheres.
	 */
	double mu_ex;

public:
	// constructor
	CSolvent(double box_vol);

	// destructor
	~CSolvent(){}

	void setParameters(const CInputParameters_t &parameters,double box_vol);

	double getCharge() const {return charge;}
	double getRadius()const {return radius;}
	std::string getLabel() const {return label;}
	double getPackingFraction()const {return packing_fraction;}
	double getDielectricConstant() const {return dielectric;}
	double getTargetNumber() const {return target_num;}
	double getExcessChemicalPotential()const {return mu_ex;}
	double getConcentration() const {return conc;}// in Molar units
};





typedef class CSolvent CSolvent_t;

#endif /* SOLVENT_TYPE_HPP_ */
