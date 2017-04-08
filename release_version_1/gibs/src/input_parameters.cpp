/**
 * input_parameters.cpp
 *
 *  Created on: Apr 13, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file input_parameters.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains definitions for member functions of class CInputParameters
 */

#include "input_parameters.hpp"

/**
 * @brief Reads the simulation run parameters from /inputfiles/inputparameters.in
 *
 *
 * Functions:
 * 1. Reads the parameter values from inputparameters.in
 * 2. Computes values for the following parameters:
 * 	- number of all possible ion pairs (nionpair)
 *  - rdf bin size (delr_rdf) if it is not read or if its value is 0.0.
 *  - rdf bin number (nrdfbin) if a non-zero bin size is read from the input parameter file.
 *
 *
 * @return void
 */
void CInputParameters::readFromFile() {

	this->nion_type = 0; // initialize number

	double delr_rdf;  // rdf spacing

	std::ifstream paramfile;
	const char* fname = "./inputfiles/inputparameters.in";

	std::stringstream ss_suffix;
	std::string suffix;

	this->delr_rdf = 0.0;

	std::string line;
	std::stringstream ss;
	std::string *sub_strings;
	sub_strings = new std::string[10];

	/* open the input parameter file 'simparamfile.txt'	 */
	paramfile.open(fname);

	if (paramfile.is_open()) {
		std::cout << "[read_parameter_file(readFromFile)] File " << fname << " is opened." << std::endl;
		while (getline(paramfile, line)) {

			ss.clear();
			ss << "";
			ss << line;
			int index = 0;

			while (ss >> sub_strings[index]) {
				++index;
			}

			if (sub_strings[0]=="SIMULATION_TYPE") {
				if (sub_strings[1]==
						"GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION") {
					this->simulation_type =
							GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION;
					this->run_type = 1;
				}
				if (sub_strings[1]=="GCMC_WITH_SOLUTE_ALLATOM_MODEL") {
					this->simulation_type = GCMC_WITH_SOLUTE_ALLATOM_MODEL;
					this->run_type = 2;
				}
				if (sub_strings[1]=="GCMC_WITH_SPHERICAL_MACROION") {
					this->simulation_type = GCMC_WITH_SPHERICAL_MACROION;
					this->run_type = 3;
				}
				if (sub_strings[1]=="GCMC_WITH_CYLINDRICAL_POLYION") {

					this->simulation_type = GCMC_WITH_CYLINDRICAL_POLYION;
					this->run_type = 4;
				}

			}else if(sub_strings[0]=="SOLUTE_MODEL_TYPE"){
				if(sub_strings[1]=="none"){
					this->solute_model = none;
				}else if(sub_strings[1]=="all_atom"){
					this->solute_model = all_atom;
				}else if(sub_strings[1]=="cylindrical_polyion"){
					this->solute_model = cylindrical_polyion;
				}else if(sub_strings[1]=="spherical_macroion"){
					this->solute_model = spherical_macroion;
				}

			}else if(sub_strings[0]=="USE_SPM"){
				if (sub_strings[1]=="YES"){
					this->use_spm = true;
				}else if (sub_strings[1]=="NO"){
					this->use_spm = false;
				}else{
					this->use_spm = false; // default value
				}
			}else if(sub_strings[0]=="USE_MWATER"){
				if (sub_strings[1]=="YES"){
					this->use_mwater = true;
				}else if (sub_strings[1]=="NO"){
					this->use_mwater = false;
				}else{
					this->use_mwater = false; // default value
				}
			}else if(sub_strings[0]=="SOLVENT_ION_ATTRACTION"){
				if (sub_strings[1]=="YES"){
					this->solvent_ion_attraction = true;
				}else if (sub_strings[1]=="NO"){
					this->solvent_ion_attraction = false;
				}else{
					this->solvent_ion_attraction = false; // default value
				}
			}else if(sub_strings[0]=="USE_PID"){
				if (sub_strings[1]=="YES"){
					this->use_pid = true;
				}else if (sub_strings[1]=="NO"){
					this->use_pid = false;
				}else{
					this->use_pid = false; // default value
				}
			}else if (sub_strings[0]=="CAVITY_GRID_SPACING") {
				this->cavity_grid_spacing = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="CAVITY_GRID_SEGMENTS") {
				this->num_cavitygrid_segments = atof(sub_strings[1].c_str());
			} else if(sub_strings[0]=="SOLVENT_LABEL"){
				this->solvent_label = sub_strings[1].c_str();
			}else if(sub_strings[0]=="LENNARD_JONES_CUTOFF"){
				this->lj_cutoff = atof(sub_strings[1].c_str());
			}else if(sub_strings[0]=="SOLVENT_RADIUS"){
				this->solvent_radius = atof(sub_strings[1].c_str());
			}else if(sub_strings[0]=="SOLVENT_PACKING_FRACTION"){
				this->solvent_packing_fraction = atof(sub_strings[1].c_str());
			}else if(sub_strings[0]=="SOLVENT_CONCENTRATION"){
				this->solvent_conc = atof(sub_strings[1].c_str());
			}else if(sub_strings[0]=="SOLVENT_EXCESS_CHEMICAL_POTENTIAL"){
				this->solvent_mu_ex = atof(sub_strings[1].c_str());
			}	else if(sub_strings[0]=="MAXIMUM_NUMBER_OF_PARTICLES"){
				this->maximum_particle_number = atof(sub_strings[1].c_str());
			}	else if (sub_strings[0]=="NUM_ITER") {
				this->niter = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="NUM_STEPS") {
				this->nsteps = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="NUM_GCMCCYCL") {
				this->gcmc_cycl = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="NUM_MOVCYCL") {
				this->disp_cycl = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="NUM_EQSTEPS") {
				this->neqlb_steps = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="MIN_PTCL_SPACING") {
				this->min_ptcl_spacing = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="INITIAL_PTCL_NUMBER") {
				this->npart_init = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="START_STATE") {
				this->start_state = atoi(sub_strings[1].c_str());
			}else if (sub_strings[0]=="RECORD_PARTICLETYPE_COUNT") {
				if (sub_strings[1]=="YES") {
					this->record_particletype_count = true;
				} else if (sub_strings[1]=="NO") {
					this->record_particletype_count = false;
				}
			}else if (sub_strings[0]=="RECORD_MOVE_ACCEPTANCE_RATE") {
				if (sub_strings[1]=="YES") {
					this->record_acceptance_rate = true;
				} else if (sub_strings[1]=="NO") {
					this->record_acceptance_rate = false;
				}
			}else if (sub_strings[0]=="RECORD_SYSTEM_ENERGY") {
				if (sub_strings[1]=="YES") {
					this->record_system_energy = true;
				} else if (sub_strings[1]=="NO") {
					this->record_system_energy = false;
				}
			}else if (sub_strings[0]=="RECORD_EVERY_N_STEPS") {
				this->record_every_n_steps = atoi(sub_strings[1].c_str());
			}	else if (sub_strings[0]=="NUM_RECORD_STATES") {
				this->num_record_states = atoi(sub_strings[1].c_str());
			}  else if (sub_strings[0]=="LOOKUP_PMF_DELR") {
				this->lookup_pmf_delr = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="LOOKUP_PMF_NUMR") {
				this->lookup_pmf_numr = atoi(sub_strings[1].c_str());
			}else if (sub_strings[0]=="SLOTH_SORENSEN_CORRECTION") {
				this->icorr = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="CALCULATE_SOLUTE_PARTICLE_RDF") {
				if (sub_strings[1]=="YES") {
					this->calculate_solute_particle_rdf = true;
				} else if (sub_strings[1]=="NO") {
					this->calculate_solute_particle_rdf = false;
				}
			} else if (sub_strings[0]=="CALCULATE_PARTICLE_PAIR_RDF") {
				if (sub_strings[1]=="YES") {
					this->calculate_particle_pair_rdf = true;
				} else if (sub_strings[1]=="NO") {
					this->calculate_particle_pair_rdf = false;
				}
			} else if (sub_strings[0]=="RDF_BIN_SIZE") {
				this->delr_rdf = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="NUM_RDF_BINS") {
				this->nrdfbin = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="RDF_MAX_DISTANCE") {
				this->rdf_max_distance = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="RDF_REF") {
				if (sub_strings[1]=="AXIS_OF_CYLINDER") {
					this->rdf_ref = AXIS_OF_CYLINDER;
				} else if (sub_strings[1]=="CENTER_OF_SPHERE") {
					this->rdf_ref = CENTER_OF_SPHERE;
				}
			} else if (sub_strings[0]=="RDF_REF_AXIS_LENGTH_FRACTION") {
				this->rdf_ref_axis_length_fraction = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="ION") {
				struct_iontype ion;
				ion.label = sub_strings[1];
				ion.charge = atof(sub_strings[2].c_str() );
				ion.radius = atof(sub_strings[3].c_str() );
				ion.conc = atof(sub_strings[4].c_str() );
				ion.mu_ex = atof(sub_strings[5].c_str() );


				++this->nion_type;
				ion.itype = this->nion_type;

				this->iontypes.push_back(ion);
			}  else if (sub_strings[0]=="X_LEN") {
				this->xlen = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="Y_LEN") {
				this->ylen = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="Z_LEN") {
				this->zlen = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="HX") {
				this->hx = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="HY") {
				this->hy = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="HZ") {
				this->hz = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="X_MIN") {
				this->xmin = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="Y_MIN") {
				this->ymin = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="Z_MIN") {
				this->zmin = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="TEMPERATURE") {
				this->temperature = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="SOLVENT_DIELECTRIC") {
				this->solvent_dielectric = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="SOLUTE_DIELECTRIC") {
				this->solute_dielectric = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="MACROION_SPHERE_CHARGE") {
				this->sphere_charge = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="MACROION_SPHERE_RADIUS") {
				this->sphere_radius = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="POLYION_CYLINDER_MONOMER_LENGTH") {
				this->cylinder_monomer_length = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="POLYION_CYLINDER_NUM_MONOMERS") {
				this->cylinder_num_monomers = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="POLYION_CYLINDER_CHARGE_PER_MONOMER") {
				this->cylinder_charge_per_monomer = atof(sub_strings[1].c_str()); // (in e units)
			} else if (sub_strings[0]=="POLYION_CYLINDER_RADIUS") {
				this->cylinder_radius = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="SELECT_POLYION_CYLINDER_POTENTIAL") {
				this->select_cylinder_polyion_potential = atoi(sub_strings[1].c_str());
			} else if (sub_strings[0]=="POLYION_CYLINDER_LENGTH") {
				this->cylinder_length = atof(sub_strings[1].c_str());
			} else if (sub_strings[0]=="ELECTROSTATIC_POTENTIAL_DX_MAP") {
				this->pot_dxmap_file = sub_strings[1].c_str();
			} else if (sub_strings[0]=="SOLUTE_ALLATOM_COORDINATE_FILETYPE") {
				if (sub_strings[1]=="XYZR") {
					this->solute_allatom_coordinate_filetype = XYZR;
				} else if (sub_strings[1]=="PQR") {
					this->solute_allatom_coordinate_filetype = PQR;
				}
			} else if (sub_strings[0]=="SOLUTE_ALLATOM_COORDINATE_FILENAME") {
				this->solute_allatom_coordinate_filename = sub_strings[1].c_str();
			} else if (sub_strings[0]=="CYLINDRICAL_RDF_AXIS_DIRECTION") {
				if (sub_strings[1]=="X_AXIS") {
					this->cylindrical_rdf_axis_direction = X_AXIS;
				} else if (sub_strings[1]=="Y_AXIS") {
					this->cylindrical_rdf_axis_direction = Y_AXIS;
				} else if (sub_strings[1]=="Z_AXIS") {
					this->cylindrical_rdf_axis_direction = Z_AXIS;
				}

			} else if (sub_strings[0]=="#") {
				// do nothing
			} else {
				std::cout << "[input_parameters.cpp(readFromFile)]: Keyword,"
						<< sub_strings[0].c_str() << " is not recognized." << std::endl;
				exit(0);
			}

		}


	} // end of IF loop

	else{
		std::cout
		<< "[input_parameters.cpp(readFromFile)]Error opening file: "
		<< fname << std::endl;
	}
	paramfile.close();

	delete[] sub_strings;
	sub_strings = NULL;


	/* Calculating the value of nionpair and delr_rdf */

	int nrdfbin = this->nrdfbin;

	// bin size takes precedance over number of bins, unless its value is set to zero.
	if (this->delr_rdf == 0.0) {
		// Calculate delr_rdf based on maximum radial distance (rdf_len) and number of bins (nrdfbin).

		if (this->run_type == 1) {
			delr_rdf = this->rdf_max_distance / (2.0 * nrdfbin);
		} else if (this->run_type == 3) {

			delr_rdf = (0.5 * this->rdf_max_distance - this->sphere_radius) / nrdfbin;

		} else if (this->run_type == 4) {
			delr_rdf = (0.5 * this->rdf_max_distance - this->cylinder_radius) / nrdfbin;
		}


		this->delr_rdf = delr_rdf;
	} else {
		// Calculate number of bins based on maximum radial distance (rdf_len) and bin size (delr_rdf).
		nrdfbin = (int) (this->rdf_max_distance / this->delr_rdf);
		if (this->nrdfbin != nrdfbin) {  // number of bins will be over-written if not equal to input value.
			this->nrdfbin = nrdfbin;
		}
	}
	/* Compute the number of ion pairs based on the number of ion types */
	this->nionpair = nion_type * (nion_type + 1) / 2;

	if (nrdfbin > MAXBINNUM) {
		std::cout
		<< "[input_parameters.cpp (readFromFile)]Error: The number of rdf bins is greater than the maximum allowed value"
		<< MAXBINNUM << std::endl;
		std::cout << "Run stopped!" << std::endl;
		exit(0);
	}

	std::cout
	<< "[input_parameters.cpp (readFromFile)]:  Completed readFromFile routine."
	<< std::endl;

}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Prints out the values of input parameter struct members to screen
 *
 * @return void
 */

void CInputParameters::printOut() {

	int ind, niontype;

	std::string simulation_types[] = { "GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION",
			"GCMC_WITH_SOLUTE_ATOM_MODEL", "GCMC_WITH_SPHERICAL_MACROION",
			"GCMC_WITH_CYLINDRICAL_POLYION" };

	std::string solute_molecule_file_types[] = { "XYZR", "PQR" };

	std::string axis_direction[] = { "X_AXIS", "Y_AXIS", "Z_AXIS" };

	std::cout
	<< "[input_parameters.cpp(printOut)] input parameter values"
	<< std::endl;

	std::cout << "SIMULATION_TYPE = " << simulation_types[this->simulation_type]
														  << std::endl;
	std::cout << "RUN_TYPE = " << this->run_type << std::endl;

	std::cout << "CAVITY_GRID_SPACING = " << this->cavity_grid_spacing << std::endl;
	std::cout << "SOLVENT_LABEL: " << this->solvent_label << std::endl;
	std::cout << "SOLVENT_RADIUS (Angstrom) = " << this->solvent_radius << std::endl;
	std::cout << "SOLVENT_PACKING_FRACTION = " << this->solvent_packing_fraction << std::endl;
	std::cout << "SOLVENT_CONCENTRATION (M) = " << this->solvent_conc << std::endl;
	std::cout << "SOLVENT_EXCESS_CHEMICAL_POTENTIAL (kcal/mol) = "  << this->solvent_mu_ex << std::endl;
	std::cout << "MAXIMUM NUMBER OF PARTICLES (for allocating array and vector size): "  << this->maximum_particle_number << std::endl;

	std::cout << "NUM_ITER (niter) = " << this->niter << std::endl;
	std::cout << "NUM_STEPS (nsteps) = " << this->nsteps << std::endl;
	std::cout << "NUM_GCMCCYCL (gcmccycl) = " << this->gcmc_cycl << std::endl;
	std::cout << "NUM_MOVCYCL (movcycl) = " << this->disp_cycl << std::endl;
	std::cout << "NUM_EQSTEPS (eqsteps) = " << this->neqlb_steps << std::endl;
	niontype = this->nion_type;
	std::cout << "NUM_IONTYPES (niontype) = " << niontype << std::endl;
	for (ind = 0; ind < niontype; ind++) {
		std::cout << "ion label: " << this->iontypes[ind].label << std::endl;
		std::cout << "ion type = " << this->iontypes[ind].itype << std::endl;
		std::cout << "ion valency (charge in e units) = "
				<< this->iontypes[ind].charge << std::endl;
		std::cout << "ion radius (in Angstrom) = " << this->iontypes[ind].radius
				<< std::endl;
		std::cout << "ion conc (in M) = " << this->iontypes[ind].conc << std::endl;
		std::cout << "ion muex (in kcal/mol) = " << this->iontypes[ind].mu_ex
				<< std::endl;
		std::cout << "ion num = " << this->iontypes[ind].num << std::endl;
		std::cout << "LJ potential well depth (in kcal/mol) = "
				<< this->iontypes[ind].lj_epsilon << std::endl;
		std::cout << "LJ collision diameter (in Angstrom) = "
				<< this->iontypes[ind].lj_sigma << std::endl;

		std::cout << "---" << std::endl;

	}

	std::cout << "X_LEN (xlen, in Angstrom) = " << this->xlen << std::endl;
	std::cout << "Y_LEN (ylen, in Angstrom) = " << this->ylen << std::endl;
	std::cout << "Z_LEN (zlen, in Angstrom) = " << this->zlen << std::endl;
	std::cout << "HX (hx, in Angstrom) = " << this->hx << std::endl;
	std::cout << "HY (hy, in Angstrom) = " << this->hy << std::endl;
	std::cout << "HZ (hz, in Angstrom) = " << this->hz << std::endl;
	std::cout << "X_MIN (xmin, in Angstrom) = " << this->xmin << std::endl;
	std::cout << "Y_MIN (ymin, in Angstrom) = " << this->ymin << std::endl;
	std::cout << "Z_MIN (zmin, in Angstrom) = " << this->zmin << std::endl;
	std::cout << "TEMPERATURE (temp,in Kelvin) = " << this->temperature << std::endl;
	std::cout << "SOLVENT_DIELECTRIC (solvdiel) = " << this->solvent_dielectric
			<< std::endl;
	std::cout << "SOLUTE_DIELECTRIC (solutediel) = " << this->solute_dielectric
			<< std::endl;

	std::cout << "SLOTH_SORENSON_CORRECTION (icorr) = " << this->icorr << std::endl;


	if (this->calculate_solute_particle_rdf == true) {
		std::cout << "CALCULATE_SOLUTE_PARTICLE_RDF option is " << "YES" << std::endl;
	} else if (this->calculate_solute_particle_rdf == false) {
		std::cout << "CALCULATE_SOLUTE_PARTICLE_RDF option is " << "NO" << std::endl;
	}
	if (this->calculate_particle_pair_rdf == true) {
		std::cout << "CALCULATE_PARTICLE_PAIR_RDF option is " << "YES" << std::endl;
	} else if (this->calculate_particle_pair_rdf == false) {
		std::cout << "CALCULATE_PARTICLE_PAIR_RDF option is " << "NO" << std::endl;
	}
	std::cout << "NUM_RDF_BINS (nrdfbin) = " << this->nrdfbin << std::endl;
	std::cout << "RDF_BIN_SIZE (delr_rdf) = " << this->delr_rdf << std::endl;
	std::cout << "RDF_MAX_DISTANCE in Angstrom = " << this->rdf_max_distance
			<< std::endl;
	if (this->rdf_ref == AXIS_OF_CYLINDER) {
		std::cout << "RDF_REF option is " << "AXIS_OF_CYLINDER" << std::endl;
	} else if (this->rdf_ref == CENTER_OF_SPHERE) {
		std::cout << "RDF_REF option is " << "CENTER_OF_SPHERE" << std::endl;
	}
	std::cout << "ELECTROSTATIC_POTENTIAL_DX_MAP (pot_dxmap_file) = "
			<< this->pot_dxmap_file << std::endl;

	if (this->run_type == 2) {
		std::cout << "SOLUTE_ALLATOM_COORDINATE_FILE_TYPE = "
				<< solute_molecule_file_types[this->solute_allatom_coordinate_filetype]
											  << std::endl;
		std::cout << "SOLUTE_ALLATOM_COORDINATE_FILE_NAME = "
				<< this->solute_allatom_coordinate_filename << std::endl;
		std::cout << "CYLINDRICAL_RDF_AXIS_DIRECTION = "
				<< axis_direction[this->cylindrical_rdf_axis_direction]
								  << std::endl;
		std::cout << "RDF_REF_AXIS_LENGTH_FRACTION= "
				<< this->rdf_ref_axis_length_fraction << std::endl;
		std::cout
		<< "NUM_RECORD_STATES (number of states to record after equilibration) = "
		<< this->num_record_states << std::endl;

	} else if (this->run_type == 3) {
		std::cout << "MACROION_SPHERE_CHARGE (sphere_charge, e units) = "
				<< this->sphere_charge << std::endl;
		std::cout << "MACROION_SPHERE_RADIUS (sphere_radius, in Angstrom units) = "
				<< this->sphere_radius << std::endl;

	} else if (this->run_type == 4) {
		std::cout
		<< "POLYION_CYLINDER_MONOMER_LENGTH (cylinder_monomer_length, in Angstrom units) = "
		<< this->cylinder_monomer_length << std::endl;
		std::cout << "POLYION_CYLINDER_NUM_MONOMERS (cylinder_num_monomers)  = "
				<< this->cylinder_num_monomers << std::endl;
		std::cout
		<< "POLYION_CYLINDER_CHARGE_PER_MONOMER (cylinder_charge_per_monomer, in e units) = "
		<< this->cylinder_charge_per_monomer << std::endl;
		std::cout
		<< "POLYION_CYLINDER_RADIUS (cylinder_radius, in Angstrom units) = "
		<< this->cylinder_radius << std::endl;
		std::cout
		<< "SELECT_POLYION_CYLINDER_POTENTIAL (select_cylinder_polyion_potential) = "
		<< this->select_cylinder_polyion_potential << std::endl;
		std::cout << "POLYION_CYLINDER_LENGTH (cylinder_length (in A)) = "
				<< this->cylinder_length << std::endl;
		std::cout << "Number of RDF segments for the cylinder polyion = "
				<< this->rdf_cylinder_num_segments << std::endl;
	}

	std::cout
	<< "[input_parameters.cpp(printOut)]: MIN_PTCL_SPACING (minimum spacing between particles placed initially in the box, in Angstrom units) = "
	<< this->min_ptcl_spacing << std::endl;

	std::cout << "[input_parameters.cpp(printOut)]: INITIAL_PTCL_NUMBER (ipart) = " << this->npart_init << std::endl;
	std::cout << "[input_parameters.cpp(printOut)]: START_STATE (start state) = " << this->start_state << std::endl;

	std::cout << "[input_parameters.cpp(printOut)]: Minimum image convention is used. "<< std::endl;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
