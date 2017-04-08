/**
 * ion_type.cpp
 *
 *  Created on: Aug 5, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file ion_type.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains definitions for member functions of class CIonType
 *
 */


#include "ion_type.hpp"

/**
 * @brief Adds parameter values for each ion type
 *
 * @param[in]		i(int) 		Index of an element of iontypes struct vector
 * @param[in,out]	parameters(CInputParameters_t)	Input parameters
 * @param[in]		box_volume(double)	Box volume in Angstrom cube unit
 *
 * @return void
 */
void CIonType::addIonType(int i,const CInputParameters_t &parameters,double box_vol){
	struct_iontype iont;

	iont.radius = parameters.iontypes[i].radius;
			iont.conc = parameters.iontypes[i].conc;
			iont.num = 0;//parameters.iontypes[i].num;
			iont.charge = parameters.iontypes[i].charge;
			iont.itype = parameters.iontypes[i].itype;
			iont.mu_ex = parameters.iontypes[i].mu_ex;
			iont.label = parameters.iontypes[i].label;
			//iont.lj_epsilon = parameters.iontypes[i].lj_epsilon;
			//iont.lj_sigma = parameters.iontypes[i].lj_sigma;


			iont.mu_old = kTLnConc(iont.conc, parameters.temperature) + iont.mu_ex;
						iont.mu_new = iont.mu_old;
						iont.nion_avg = 0;//concMolarityToNum(iont.conc, box_vol);
						iont.cion_avg = 0;//iont.conc;

	this->iontypes[i]=iont;
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
