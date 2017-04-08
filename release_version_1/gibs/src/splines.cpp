/*
 * splines.cpp
 *
 *  Created on: Aug 11, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file splines.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief	Function definitions for spline interpolation
 */

#include "splines.hpp"
/**
 * @brief Performs a binary search for the index of an element with a value in a vector
 * 			of floating point values
 */
/*
int binarySearch(std::vector<double> &array, int n, double value) {

	int i, index, imid, imin, imax;
	int diff;

	index = -1;
	if (value == array[n - 1]) {
		index = n - 1;
	} else {
		imid = (n - 1) / 2;
		diff = n - 1 - imid;
		imin = 0;
		imax = n - 1;
		while (diff > 1) {
			if (value >= array[imid] && value < array[imax]) {
				imin = imid;
				imid = int((imid + imax) / 2);
				diff = imax - imin;
			} else if (value < array[imid] && value >= array[imin]) {
				imax = imid;
				imid = int((imid + imin) / 2);
				diff = imax - imin;
			}
		}
		index = imin;
	}

	return index;
}
*/




/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/*
 * @brief Calculates the coefficients of cubic polynomial that approximates a function
 * 			in the interval, x_{i} = 0 to x_{i+1} = 1.
 *
 *
 * f(x) = a0 + a1*x + a2*x^2 + a3*x^3
 * Coefficient formulas taken from reference: http://www.paulinternet.nl/?page=bicubic
 * 	x = x_{i-1}, x_{i}, x_{i+1}, and x_{i+2}
 * 	e.g., f(x_{i}) = p{i}
 *
 * @param[in,out] 	a(double[])  One dimensional array of length 4, which stores the coefficients of the cubic
 * 				polynomial
 * @param[in]	p(double[])  One dimensional array of length 4, which stores the values at the grid points,
 *
 * @return void
 */
void cspline1D_coeff(double a[], double p[]) {

	a[0] = p[1];
	a[1] = 0.5 * (p[2] - p[0]);
	a[2] = p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3];
	a[3] = -0.5 * p[0] + 1.5 * p[1] - 1.5 * p[2] + 0.5 * p[3];

}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/**
 * @brief Interpolates the value at any point (x,y) between (0,0) and (1,1)
 * 		  using bicubic spline interpolation
 *
 * @param[in] a(double**) 4 x 4 array, which stores the coefficients of the bicubic polynomial
 * @param[in] x(double) x-coordinate in Angstrom units
 * @param[in] y(double) y-coordinate in Angstrom units
 *
 * @return Interpolated value (double) at (x,y)
 */
double bcspline_interp(double **a, double x, double y) {

	int i, j;

	double value = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			value += a[i][j] * pow(x, i) * pow(y, j);
		}
	}

	return value;
}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */


/**
 * @brief Calculates the coefficients of bicubic polynomial that approximates a function
 * 		in the interval, x = 0 to x = 1, and y = 0 to y = 1.
 *
 * f(x,y) = sum [a_ij (x_i)^i (y_j)^j] : i=0->3;j=0->3
 * Coefficient formulas taken from reference: http://www.paulinternet.nl/?page=bicubic
 *
 * @param[in,out] 	a(double**)  4x4 array, which stores the coefficients of the bicubic polynomial
 * @param[in]	p(double**)  4x4 array, which stores the values at the grid points,
 *
 * @return A zero integer value
 */

int cspline2D_coeff(double **a, double **p) {

	a[0][0] = p[1][1];
	a[0][1] = -0.5 * p[1][0] + 0.5 * p[1][2];
	a[0][2] = p[1][0] - 2.5 * p[1][1] + 2.0 * p[1][2] - 0.5 * p[1][3];
	a[0][3] = -0.5 * p[1][0] + 2.5 * p[1][1] - 2.5 * p[1][2] + 0.5 * p[1][3];

	a[1][0] = -0.5 * p[0][1] + 0.5 * p[2][1];
	a[1][1] = 0.25 * p[0][0] - 0.25 * p[0][2] - 0.25 * p[2][0] + 0.25 * p[2][2];
	a[1][2] = -0.5 * p[0][0] + 1.25 * p[0][1] - p[0][2] + 0.25 * p[0][3]
			+ 0.5 * p[2][0] - 1.25 * p[2][1] + p[2][2] - 0.25 * p[2][3];
	a[1][3] = 0.25 * p[0][0] - 0.75 * p[0][1] + 0.75 * p[0][2] - 0.25 * p[0][3]
			- 0.25 * p[2][0] + 0.75 * p[2][1] - 0.75 * p[2][2] + 0.25 * p[2][3];

	a[2][0] = p[0][1] - 2.5 * p[1][1] + 2.0 * p[2][1] - 0.5 * p[3][1];
	a[2][1] = -0.5 * p[0][0] + 0.5 * p[0][2] + 1.25 * p[1][0] - 1.25 * p[1][2]
			- p[2][0] + p[2][2] + 0.25 * p[3][0] - 0.25 * p[3][2];
	a[2][2] = p[0][0] - 2.5 * p[0][1] + 2.0 * p[0][2] - 0.5 * p[0][3]
			- 2.5 * p[1][0] + 6.25 * p[1][1] - 5.0 * p[1][2] + 1.25 * p[1][3]
			+ 2.0 * p[2][0] - 5.0 * p[2][1] + 4.0 * p[2][2] - p[2][3]
			- 0.5 * p[3][0] + 1.25 * p[3][1] - p[3][2] + 0.25 * p[3][3];

	a[2][3] = -0.5 * p[0][0] + 1.5 * p[0][1] - 1.5 * p[0][2] + 0.5 * p[0][3]
			+ 1.25 * p[1][0] - 3.75 * p[1][1] + 3.75 * p[1][2] - 1.25 * p[1][3]
			- p[2][0] + 3.0 * p[2][1] - 3.0 * p[2][2] + p[2][3] + 0.25 * p[3][0]
			- 0.75 * p[3][1] + 0.75 * p[3][2] - 0.25 * p[3][3];

	a[3][0] = -0.5 * p[0][1] + 1.5 * p[1][1] - 1.5 * p[2][1] + 0.5 * p[3][1];
	a[3][1] = 0.25 * p[0][0] - 0.25 * p[0][2] - 0.75 * p[1][0] + 0.75 * p[1][2]
			+ 0.75 * p[2][0] - 0.75 * p[2][2] - 0.25 * p[3][0] + 0.25 * p[3][2];

	a[3][2] = -0.5 * p[0][0] + 1.25 * p[0][1] - p[0][2] + 0.25 * p[0][3]
			+ 1.5 * p[1][0] - 3.75 * p[1][1] + 3.0 * p[1][2] - 0.75 * p[1][3]
			- 1.5 * p[2][0] + 3.75 * p[2][1] - 3.0 * p[2][2] + 0.75 * p[2][3]
			+ 0.5 * p[3][0] - 1.25 * p[3][1] + p[3][2] - 0.25 * p[3][3];
	a[3][3] = 0.25 * p[0][0] - 0.75 * p[0][1] + 0.75 * p[0][2] - 0.25 * p[0][3]
			- 0.75 * p[1][0] + 2.25 * p[1][1] - 2.25 * p[1][2] + 0.75 * p[1][3]
			+ 0.75 * p[2][0] - 2.25 * p[2][1] + 2.25 * p[2][2] - 0.75 * p[2][3]
			- 0.25 * p[3][0] + 0.75 * p[3][1] - 0.75 * p[3][2] + 0.25 * p[3][3];

	return 0;
}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxv */
/**
 * @brief Sets up 4x1 array of grid point values for computing
 * 			the coefficients of 1D cubic polynomial
 *
 *
 * @param[in] data(std::vector<double>)	Vector of data points
 * @param[in,out] coeff(std::vector<spline_coeff_struct>)	coefficients
 * @param[in] n(int)	Number of data points(data)
 *
 * @return void
 */

void create_splineCoeff(std::vector<double> &data, std::vector<spline_coeff_struct> &coeff, int n) {


	double p[4];

	for (int ind1 = 0; ind1 < n - 1; ind1++) {
		p[1] = data[ind1];
		p[2] = data[ind1 + 1];

		if (ind1 > 0 && ind1 < n - 2) {
			p[0] = data[ind1 - 1];
			p[3] = data[ind1 + 2];
		} else if (ind1 == 0) {
			p[0] = 2 * p[1] - p[2];
			p[3] = data[ind1 + 2];

		} else if (ind1 == n - 2) {
			p[0] = data[ind1 - 1];
			p[3] = 2.0 * p[2] - p[1];
		}

		cspline1D_coeff(coeff[ind1].cubic, p);
	}

}
