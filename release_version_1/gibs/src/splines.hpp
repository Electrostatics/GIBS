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
