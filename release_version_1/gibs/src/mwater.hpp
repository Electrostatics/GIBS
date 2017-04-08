/*
 * mwater.hpp
 *
 *  Created on: Mar 1, 2016
 *      Author:  Dennis G. Thomas
 *
 *  @file mwater.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Class definitions and parameters of Monoatomic water (MWater) model
 *
 *
 * Reference for MWater model:
 * 	Valeria Molinero* and Emily B. Moore, "Water Modeled As an Intermediate Element
 * 	between Carbon and Silicon", J. Phys. Chem. B, 2009, 113 (13), pp 4008–4016
 */

#ifndef MWATER_HPP_
#define MWATER_HPP_


#include "generic.hpp"


class CMWater{

public:
	double A;
	double B;
	int p;
	int q;
	double gamma;
	double a;
	double costheta0;
	double lamda;
	double epsilon;	// kcal/mol
	double sigma;	// Angstrom
	double asigma; // asigma = a x sigma
	double AepsilonB; // AepsilonB = A x epsilon x B
	double Aepsilon; // Aepsilon = A x epsilon
	double lamdaepsilon; // lamdaepsilon = lamda x epsilon
	double gammasigma; // gammasigma  = gamma x sigma

	std::vector<int> solvent_indices;

	CMWater(){

		this->A = 7.049556277;
		this->B = 0.6022245584;
		this->p = 4;
		this->q = 0;
		this->gamma = 1.2;
		this->a = 1.8;
		this->costheta0 =-0.333333333333; // theta0 = 109.47; costheta0 = -1/3.
		this->lamda = 23.15;
		this->epsilon = 6.189;   // units: kcal/mol
		this->sigma = 2.3925;
		this->asigma = a * sigma;
		this->AepsilonB = A * epsilon * B;
		this->Aepsilon = A * epsilon;
		this->lamdaepsilon = lamda * epsilon;
		this->gammasigma = gamma * sigma;

		this->solvent_indices.resize(1000);


	}

	~CMWater(){}

};


#endif /* MWATER_HPP_ */
