/**
 * ion_type.hpp
 *
 *  Created on: Apr 13, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file ion_type.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains class declaration and function prototypes/definitions for ion types
 */

#ifndef ION_TYPE_HPP_
#define ION_TYPE_HPP_


#include "generic.hpp"
#include "input_parameters.hpp"
#include "iontype_structs.hpp"
#include "convert_functions.hpp"

/**
 * @brief Class for storing parameters and chemical
 * potential values of each ion type
 *
 * Use defined type CIonType_t for class declarations
 */
class CIonType{

	// Vector variable of type struct_iontype
	std::vector<struct_iontype> iontypes;

public:
	// constructor
	CIonType(int n){iontypes.resize(n);}

	// destructor
	~CIonType(){}

/**
 * @brief Function to add parameters for each ion type
 *
 */
void addIonType(int i,const CInputParameters_t &parameters,double box_vol);

/**
 * @brief Function to update the number of an ion species during simualation
 *
 */
void updateIonTypeNum(int i,int num){iontypes[i].num = num;}

// get methods
struct_iontype getIonType(int i) const {return iontypes[i];}

/**
 * @brief Get the ion concentration
 *
 * @param[in] i	The index (integer) corresponding to the position of
 * 					an iontype in the iontypes vector
 *
 * @return The concentration (double) of the ion species in Molar units
 */
double getIonConc(int i) const{return iontypes[i].conc;}

};

typedef class CIonType CIonType_t;

#endif /* ION_TYPE_HPP_ */
