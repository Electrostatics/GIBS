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
 * box.hpp
 *
 *  Created on: Apr 7, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @brief	Class definitions and function prototypes for simulation box
 *
 *  @file box.hpp
 *  @author Dennis G. Thomas
 */

#ifndef BOX_HPP_
#define BOX_HPP_

#include "generic.hpp"
#include "input_parameters.hpp"

/**
 * @brief Class declaration for storing and accessing the parameters of a simulation box (rectangular)
 *
 * The x and y axes are considered to be in the plane of the paper. The
 * z-axis is directed outward, perpendicular from the plane of the paper.
 *
 * Use defined type	CBox_t for class declarations
 */
class CBox{
	// Integer variables to store the number of grid points along x, y, and z dimensions
	int nx, ny, nz;
	/*Floating point variables to store the uniform grid spacing
	along x, y, and z dimensions, in Angstrom units */
	double hx,hy,hz;
	// Floating point variables to store the box lengths along x, y, and z in Angstrom
	double xlen, ylen, zlen;
	// Floating point variable to store the box volume in Angstrom cube units
	double vol;
	// Floating point variables to store the minimum value of x,y, and z in Angstrom units
	double xmin, ymin, zmin;
	// Floating point variables to store the (x,y,z) coordinates of the box center, in Angstrom units
	double xcenter,ycenter,zcenter;
	// Floating point variables to store the maximum value of x, y, and z in Angstrom units
	double xmax,ymax,zmax;

public:

	/********** set methods **********/
	/**
	 * @brief Function to set box dimensions
	 *
	 * @param[in]	parameters(CInputParameters)	Input parameters
	 *
	 * @return void
	 */
	void setBoxDimensions(const CInputParameters &parameters);

	/********** get methods **********/
	/**
	 * @brief Function to get the number of grid points along x-axis
	 *
	 * @param[out] nx(int) Number of grid points along x-axis
	 *
	 * @return	Number of grid points (integer) along x-axis
	 */
	int getNx() const {return nx;}
	/**
	 * @brief Function to get the number of grid points along y-axis
	 *
	 * @param[out] ny(int) Number of grid points along y-axis
	 *
	 * @return	Number of grid points (integer) along y-axis
	 */
	int getNy()const {return ny;}
	/**
	 * @brief Function to get the number of grid points along z-axis
	 *
	 * @param[out] ny(int) Number of grid points along z-axis
	 *
	 * @return	Number of grid points (integer) along z-axis
	 */

	int getNz()const {return nz;}
	/**
	 * @brief Function to get the x-axis grid spacing in Angstrom units
	 *
	 * @param[out] hx(double) x-axis grid spacing in Angstrom units
	 *
	 * @return 	x-axis grid spacing (double) in Angstrom units
	 */
	double getBoxHx() const {return hx;}
	/**
	 * @brief Function to get the y-axis grid spacing in Angstrom units
	 *
	 * @param[out] hy(double) y-axis grid spacing in Angstrom units
	 *
	 * @return 	y-axis grid spacing (double) in Angstrom units
	 */
	double getBoxHy()const {return hy;}
	/**
	 * @brief Function to get the z-axis grid spacing in Angstrom units
	 *
	 * @param[out] hz(double) z-axis grid spacing in Angstrom units
	 *
	 * @return 	z-axis grid spacing (double) in Angstrom units
	 */
	double getBoxHz()const {return hz;}
	/**
	 * @brief Function to get the box length along x-axis
	 *
	 * @param[out] xlen(double) x-axis box length in Angstrom units
	 *
	 * @return 	x-axis box length (double) in Angstrom units
	 *
	 */
	double getBoxXLen()const {return xlen;}

	/**
	 * @brief Function to get the box length along y-axis
	 *
	 *
	 * @param[out] ylen(double) y-axis box length in Angstrom units
	 *
	 * @return y-axis box length (double) in Angstrom units
	 */

	double getBoxYLen()const {return ylen;}
	/**
	 * @brief Function to get the box length along z-axis
	 *
	 *
	 * @param[out] zlen(double) z-axis box length in Angstrom units
	 *
	 * @return z-axis box length (double) in Angstrom units
	 */
	double getBoxZLen()const {return zlen;}

	/**
	 * @brief Function to get the box volume
	 *
	 * @param[out] vol(double) box volume in Angstrom cube units
	 *
	 * @return box volume (double) in Angstrom cube units
	 */

	double getBoxVol()const {return vol;}

	/**
	 * Function to get the lowest coordinate value along x-axis box dimension
	 *
	 *
	 * @param[out] xmin(double) lowest x-coordinate value in Angstrom units
	 *
	 * @return lowest x-coordinate value (double) in Angstrom units
	 */


	double getBoxXMin()const {return xmin;}

	/**
	 * @brief Function to get the lowest coordinate value along y-axis box dimension
	 *
	 * @param[out] ymin(double) lowest y-coordinate value in Angstrom units
	 *
	 * @return lowest y-coordinate value (double) in Angstrom units
	 */


	double getBoxYMin()const {return ymin;}
	/**
	 * @brief Function to get the lowest coordinate value along z-axis box dimension
	 *
	 *
	 * @param[out] zmin(double) lowest x-coordinate value in Angstrom units
	 *
	 * @return lowest x-coordinate value (double) in Angstrom units
	 */
	double getBoxZMin()const {return zmin;}
	/**
	 * @brief Function to get the highest coordinate value along x-axis box dimension
	 *
	 *
	 * @param[out] xmax(double) highest x-coordinate value in Angstrom units
	 *
	 * @return highest x-coordinate value (double) in Angstrom units
	 */
	double getBoxXMax()const {return xmax;}
	/**
	 * @brief Function to get the highest coordinate value along y-axis box dimension
	 *
	 *
	 * @param[out] ymax(double) highest y-coordinate value in Angstrom units
	 *
	 * @return highest y-coordinate value (double) in Angstrom units
	 */
	double getBoxYMax()const {return ymax;}
	/**
	 * @brief Function to get the highest coordinate value along z-axis box dimension
	 *
	 * @param[out] zmax(double) highest z-coordinate value in Angstrom units
	 *
	 * @return highest z-coordinate value (double) in Angstrom units
	 */
	double getBoxZMax()const {return zmax;}
	/**
	 * @brief Function to get the x-axis coordinate value at the box center
	 *
	 * @param[out] xcenter(double) x-axis coordinate value at the box center in Angstrom units
	 *
	 * @return x-axis coordinate value at the box center (double) in Angstrom units
	 */
	double getBoxXCenter()const {return xcenter;}
	/**
	 * @brief Function to get the y-axis coordinate value at the box center
	 *
	 * @param[out] ycenter(double) y-axis coordinate value at the box center in Angstrom units
	 *
	 * @return y-axis coordinate value at the box center (double) in Angstrom units
	 */

	double getBoxYCenter()const {return ycenter;}
	/**
	 * @brief Function to get the z-axis coordinate value at the box center
	 *
	 *
	 * @param[out] zcenter(double) z-axis coordinate value at the box center in Angstrom units
	 *
	 * @return z-axis coordinate value at the box center (double) in Angstrom units
	 */

	double getBoxZCenter()const {return zcenter;}


	/********** write methods **********/

	/**
	 * @brief Function to write out the box parameters on the screen
	 *
	 * @return void
	 */
	void writeBoxDimensions();

};
typedef class CBox CBox_t;


#endif /* BOX_HPP_ */
