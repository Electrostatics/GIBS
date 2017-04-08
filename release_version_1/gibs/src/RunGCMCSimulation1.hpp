/*
 * RunGCMCSimulation1.hpp
 *
 *  Created on: Apr 7, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file RunGCMCSimulation1.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains the prototypes of functions defined in RunGCMCSimulation1.cpp
 *
 */

#ifndef RUNGCMCSIMULATION1_HPP_
#define RUNGCMCSIMULATION1_HPP_

#include "generic.hpp"
#include "box.hpp"
#include "ion_type.hpp"
#include "particle_type.hpp"
#include "solvent_type.hpp"
#include "input_parameters.hpp"
#include "cavity_grid.hpp"
#include "particle_state.hpp"
#include "particlepair_type.hpp"
#include "rdf_methods.hpp"

#include "solute_type.hpp"

#include "gcmc_calc1.hpp"

#include "mwater.hpp"
/**
 * The driver program for the GIBS simulations.
 */
void runSimulation_1();

#endif /* RUNGCMCSIMULATION1_HPP_ */
