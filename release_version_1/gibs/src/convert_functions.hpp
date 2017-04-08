/*
 * convert_functions.hpp
 *
 *  Created on: Aug 12, 2015
 *      Author:  Dennis G. Thomas
 *
 *
 *  @file convert_functions.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains functions for carrying out conversions from one quantity to another
 *
 */

#ifndef CONVERT_FUNCTIONS_HPP_
#define CONVERT_FUNCTIONS_HPP_

/**
 * @brief Computes the ideal gas chemical potential given the concentration and temperature
 *
 *
 * @param[in]	conc(double)	Concentration in Molar units
 * @param[in]	temp(double)	Temperature in Kelvin units
 *
 * @return Ideal gas chemical potential (double) in kcal/mol
 */
inline double kTLnConc(double conc, double temp) {

	//mu = k_B*avognum*convJtoCal*0.001*temp*log(conc);

	//mu = 0.001985862*temp*log(conc); // kcal/mol (until version 3_17_2013)
	return (0.00198719137 * temp * log(conc)); // kcal/mol

}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/**
 * @brief Converts number of particles (e.g, ions) to molar concentration
 *
 * @param[in]	num(double)		Number of particles
 * @param[in]	vol(double)		Simulation box volume in Angstrom cube
 *
 * @return Molar concentration (double)
 */

inline double numToConcMolarity(double num, double vol) {
	return (num * 10000 / (6.022045 * vol));
}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/**
 * @brief Converts molarity to number
 *
 * @param[in]	conc(double)	Particle concentration in Molar units
 * @param[in]	vol(double)		Simulation box volume in Angstrom cube
 *
 * @return Number of particles (double)
 */

inline double concMolarityToNum(double conc, double vol) {
	return (conc * 6.022045 * vol * 0.0001);
}

/**
 * @brief Converts number to mass density
 *
 * @param[in]	num(num_t)		Number of particles
 * @param[in]	vol(double)	Simulation box volume in Angstrom cube
 *
 * @return Mass density (double) in g/cm^3 units
 */

/**
 * @brief Converts number to mass concentration (density)
 *
 * @param[in]	num(num_t) number of particles
 * @param[in]	vol(double) volume of simulation box
 * @param[in]	molmass (double) molar mass of particle type in g/mol
 *
 * @return Mass density or concentration (double) in g/cm^3 units
 */

template <typename num_t>
inline double numToDensity(num_t num, double vol,double molmass){
	// density = n*18 x 10^-3/(NA x 10e-27 x V)
	// mass density = num x molwt/(vol x 10^(-24) x NA)
	//				 = num x molwt x 10/(vol x 6.022045)

//return (29.8901785*num/vol);
return (num*molmass*10/(vol*6.022045));

}

/**
 * @brief Converts molar concentration to mass concentration (density)
 *
 * @param[in]	conc(double)	Molar concentration of particles
 * @param[in]	molmass(double)	Molar mass of particle type (g/mol)
 *
 * @return Mass density or concentration (double) in g/cm^3 units
 */

inline double concToDensity(double conc, double molmass){
	// density = n*18 x 10^-3/(NA x 10e-27 x V)
	// mass density = conc x molwt / 1000
//return (29.8901785*conc * 6.022045 * 0.0001);
return (conc * molmass / 1000.0);
}

#endif /* CONVERT_FUNCTIONS_HPP_ */
