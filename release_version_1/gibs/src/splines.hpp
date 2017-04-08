/**
 * splines.hpp
 *
 *  Created on: Aug 11, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file splines.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Data structures, inline functions and function prototypes
 *  		for spline interpolation
 */

#ifndef SPLINES_HPP_
#define SPLINES_HPP_

#include "generic.hpp"

struct spline_coeff_struct {
	double cubic[4];
	double **bicubic;
};
typedef struct spline_coeff_struct struct_spline_coeff;



inline double cspline_interp(double a[], double x) {
	/**
	 *cspline_interp:  cubic spline interpolation
	 * Function interpolates the value at any point (x) between 0 and 1.
	 *
	 */
	return(a[0] + a[1] * x + a[2] * x * x + a[3] * x * x * x);

}
//int binarySearch(std::vector<double> &array, int n, double value);
void cspline1D_coeff(double a[], double p[]);
void create_splineCoeff(std::vector<double> &data, std::vector<spline_coeff_struct> &coeff, int n);

#endif /* SPLINES_HPP_ */
