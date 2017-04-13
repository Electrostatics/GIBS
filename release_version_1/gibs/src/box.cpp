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
 * box.cpp
 *
 *  Created on: Apr 7, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @brief	CBox class member function definitions
 *
 *  @file box.cpp
 *  @author Dennis G. Thomas
 */

#include "box.hpp"

/**
 * @brief	Sets the box dimensions
 *
 * @param[in] parameters Input parameters (instance of class CInputParameters)
 *
 * @return void
 *
 */
void CBox::setBoxDimensions(const CInputParameters &parameters){
	hx = parameters.hx;
	hy = parameters.hy;
	hz = parameters.hz;


	xmin=parameters.xmin;
	ymin=parameters.ymin;
	zmin=parameters.zmin;

	xlen = parameters.xlen;
	ylen = parameters.ylen;
	zlen = parameters.zlen;


	if (fmod(xlen,hx) ==0)
		nx = xlen/hx+1;
	else {
		nx = int(xlen/hx)+2;
		xlen = hx*(nx-1);
	}

	if(fmod(ylen,hy) == 0)
		ny = ylen/hy + 1;

	else {
		ny = int(ylen/hy)+2;
		ylen = hy*(ny-1);
	}

	if (fmod(zlen,hz) ==0)
		nz =zlen/hz + 1;

	else {
		nz = int(zlen/hz)+2;
		zlen = hz*(nz-1);
	}

	vol = xlen*ylen*zlen;

	xcenter = xmin + 0.5*xlen;
	ycenter = ymin + 0.5*ylen;
	zcenter = zmin + 0.5*zlen;

	xmax = xmin + xlen;
	ymax = ymin + ylen;
	zmax = zmin + zlen;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief	Writes out the box parameters on the screen
 *
 * @return void
 *
 */
void CBox::writeBoxDimensions(){
	std::cout << "[box.cpp: writeBoxDimensions] box.hx, box.xlen, box.nx = " << hx <<" "<< xlen << " " << nx << std::endl;
	std::cout << "[box.cpp: writeBoxDimensions] box.hy, box.ylen, box.ny = " << hy <<" "<< ylen << " " << ny << std::endl;
	std::cout << "[box.cpp: writeBoxDimensions] box.hz, box.zlen, box.nz = " << hz <<" "<< zlen << " " << nz << std::endl;

	std::cout << "[box.cpp: writeBoxDimensions] box.xmin, box.ymin, box.zmin = " << xmin <<" "<< ymin << " " << zmin << std::endl;
	std::cout << "[box.cpp: writeBoxDimensions] box.xcenter, box.ycenter, box.zcenter = " << xcenter <<" "<< ycenter << " " << zcenter << std::endl;

	std::cout <<"[box.cpp: writeBoxDimensions] box volume=" << vol << std::endl;
}
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

