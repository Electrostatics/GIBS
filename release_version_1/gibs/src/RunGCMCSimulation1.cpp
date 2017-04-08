/*
 * RunGCMCSimulation1.cpp
 *
 *  Created on: Apr 7, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file RunGCMCSimulation1.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Contains the driver program for the GIBS simulations
 *
 */

#include "RunGCMCSimulation1.hpp"

/**
 * @brief Driver program for GCMC simulations
 *
 * @return void
 */
void runSimulation_1() {

	/*-------Read input parameters for simulation-------*/

	/*A variable of type CInputParameters_t to store the simulation run parameter values*/

	CInputParameters_t parameters;
	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: READING INPUT "
			"PARAMETER FILE ---------->" << std::endl;

	/* Read parameters from input file */
	parameters.readFromFile();
	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: FINISHED!"
			<< std::endl;

	/*------- Set box dimensions ------- */
	/* A variable of type CBox_t to store the box dimensions */
	CBox_t box;
	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: "
			"SETTING BOX DIMENSIONS ---------->" << std::endl;

	/*Set the box dimensions from 'parameters' */
	box.setBoxDimensions(parameters);
	/* Write out the box dimensions to screen */
	box.writeBoxDimensions();


	/*-------Set up solute model (should be set before setting up particle state vector,
	cavity lists, and particle neighbor count lists -------*/

	// A variable of type CSoluteModel_t to store the solute structure data
	CSoluteModel_t solute(parameters,box);


	/*-------Set up cavity grid dimensions -------*/
	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: Setting up particle grid dimensions and "
			" cell index maps (3D to 1D, and vice-versa) ---------->" << std::endl;

	// A variable of type CCavityGrid_t to store the cavity grid information
	CCavityGrid_t cavity_grid(box,parameters.cavity_grid_spacing,parameters.num_cavitygrid_segments);

	/*------- Set ion type parameters -------*/
	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: SETTING UP PARAMETERS FOR EACH ION TYPE ---------->" << std::endl;
	// Integer variable for the number of ion types
	int niontype = parameters.nion_type;
	/* A class variable of type CIonType_t to store the parameters,
	 * chemical potential and the simulation number of each ion species
	 */
	CIonType_t iontypes(niontype);

	for(int i=0;i<niontype;i++){
		iontypes.addIonType(i,parameters,box.getBoxVol());
	}

	/*------- Set solvent parameters -------*/

	/*Declare a variable of type CSolvent_t to store the solvent parameters */
	CSolvent_t solvent_type(box.getBoxVol());
	solvent_type.setParameters(parameters,box.getBoxVol());
	int nsolvent = 0;
	if (parameters.use_spm == true) {
		/*
		 * Represent solvent molecules as neutral hard spheres.
		 * Particles can be ions or solvent molecule.
		 */
		nsolvent = 1;
	} else {
		/*
		 * Represent solvent as a continuum dielectric medium.
		 * All particles are ions.
		 */
		nsolvent = 0;
	}

	/*------- Set parameters of all particle types (ions, solvent) ------- */

	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: SETTING PARAMETERS OF ALL PARTICLE TYPES ---------->" << std::endl;
	/*Create list of particle indices by count and grid cell index.
	 Define grid variable for total number of neighboring particles whose centers are within a hard-sphere distance
	 from the lattice site.*/

	int nparticle_types = niontype + nsolvent; // number of particle types
	std::cout << "[RunGCMCSimulation1.cpp (runSimulation_1]: Number of particle types = " << nparticle_types << std::endl;

	/*-------Set up the parameters and cavity grid data for each particle type-------*/
	CParticleType_t particle_types(nparticle_types,cavity_grid);
	particle_types.setParameters(iontypes,solvent_type,niontype,nsolvent,parameters.temperature,box.getBoxVol());

	/*--------- Read molar mass and set additional parameters for each particle type------- */
	particle_types.readMolarMassAndSetParameters("./inputfiles/molar_mass.in",nparticle_types);

	/*-------Read and set the parameters for using the PID method-------*/
	if(parameters.use_pid == true){
		particle_types.readAndSetPIDInitValues("./inputfiles/pid_init.in");
	}

	/*-------Set particle pair labels and create mappings between vector indices and labels-------*/
	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]:Setting particle pair labels and mappings between indices and labels ---------->" << std::endl;
	/*Calculate the number of particle pair types based on the number of particle types*/
	CParticlePairType_t particlepair_types(nparticle_types);
	int nparticle_pairs = particlepair_types.getNumPairs();  // get the number of particle pair types
	/*Set up the labels for each particle pair type */
	particlepair_types.setTypeAndLabels(nparticle_types,particle_types);

	for (int i=0; i<nparticle_pairs; i++){
		std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: particle index = " <<
				i << ", " <<	particlepair_types.getLabelFromIndex(i) << std::endl;
	}

	/*-------Read the types of interaction potentials to use for each pair of particle types -------*/
	particlepair_types.readAndSetInteractionTypes("./inputfiles/particle_pair_interactions.in");

	/*------- Read the hard sphere cutoff distances for each particle pair types -------*/
	particlepair_types.readAndSetHardSphereCutoffs("./inputfiles/hard_sphere_cutoffs.in");

	/*------- Set boolean variables to determine whether to read the Lennard-Jones parameters,
	 * square well potential parameters and/or PMF look up values from respective input files-------*/
	bool read_ljfile = false;
	bool read_sqrwfile = false;
	bool read_lookup_table = false;
	for(int i=0;i<nparticle_pairs;i++){

		if(particlepair_types.getType(i).lennard_jones_potential==true){
			read_ljfile = true;
			particle_types.setLJSwitch(particlepair_types.getTypeIndex1(i)-1,true);
			particle_types.setLJSwitch(particlepair_types.getTypeIndex2(i)-1,true);

		}else if(particlepair_types.getType(i).square_well_potential==true){
			read_sqrwfile = true;
			if(particlepair_types.getType(i).label1 == parameters.solvent_label ||
					particlepair_types.getType(i).label2 == parameters.solvent_label){
				parameters.solvent_ion_attraction = true;
			}
		}else if(particlepair_types.getPMFLookTable(i)==true){
			read_lookup_table = true;
		}
	}



	std::cout << "read_ljfile = " << read_ljfile << std::endl;
	std::cout << "read_sqrwfile = " << read_sqrwfile << std::endl;

	/*------- Read the Lennard-Jones parameters if true for at least one particle pair type*/
	if(read_ljfile==true){
		// set up the cell index maps for Lennard Jones cells
		cavity_grid.setLJCellIndexMaps(parameters.lj_cutoff);

		// read LJ well depth and collision diameter from file
		particle_types.readAndSetLennardJonesParameters("./inputfiles/lennard_jones_parameters.in",nparticle_types);

		// reserve space for the LJ cell cavity indices list
		particle_types.initializeLJCellCavityIndices(nparticle_types,cavity_grid);

		particlepair_types.setLennardJonesParameters(particle_types);

		particlepair_types.setLJCutoffs(parameters.lj_cutoff);
	}

	/*------- Read the Square well potential parameters if true for at least one particle pair type*/
	if(read_sqrwfile==true){
		particlepair_types.readAndSetSquareWellParameters("./inputfiles/square_well_parameters.in");
	}

	/*------- Read the PMF look up table if true for at least one particle pair type*/

	if(read_lookup_table==true){
		particlepair_types.readPMFLookupFile(parameters);
		particlepair_types.writePMFLookupTableToFile();
	}

	/*------- Write out the parameter values to screen ------- */
	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: WRITING OUT THE SIMULATION PARAMETERS ---------->" << std::endl;
	parameters.printOut();


	// Set monoatomic water model parameters-------------------------------------------------->
	/*------- Set Monatomic water model parameters if the MWater model with SPM hard spheres are used*/
	CMWater mwater;
	if(parameters.use_mwater==true && parameters.use_spm == true){
		particlepair_types.setMWaterPotentialCutoff(mwater.asigma);
		particlepair_types.setMWaterCellIndexOffsetMatrix(cavity_grid);

	}

	/*------- Create a matrix of cell index offset values for determining the cells
	 * that are within the hard sphere distance of each particle type in a particle pair-------*/

	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: CONSTRUCTING NEIGHBOR GRID CELL OFFSET MATRIX FOR EACH PARTICLE PAIR TYPE"
			" BASED ON HARD-SPHERE DISTANCE---------->" << std::endl;

	for(int i=0;i<nparticle_pairs;i++){

		if(particlepair_types.getType(i).hard_sphere_repulsion==false){
			particlepair_types.changeHardSphereCutoff(particlepair_types.getType(i).label12,0.5);
			particlepair_types.changeHardSphereCutoff(particlepair_types.getType(i).label21,0.5);

		}
	}
	particlepair_types.setHardSphereCellIndexOffsetMatrix(cavity_grid);

	/*------- Create a matrix of cell index offset values for determining the cells
	 * that are within the Lennard-Jones potential cutoff distance of each particle
	 * type in a particle pair-------*/

	particlepair_types.setSqrWellWidthCellIndexOffsetMatrix(cavity_grid);


	/*------- Set up the particle state vector, cavity grid lists, and particle neighbor
	 * count lists -------*/

	/*Cavity list: For each particle type, we create a cavity (lattice site) list, which
	 * is a vector of 1D  particle grid cell indices.Initial cavity list consists of all the
	 *  particle grid cell indices for each particle type.*/

	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: Setting up particle state vector"
			" cavity lists, and particle neighbor count lists." << std::endl;
	CParticleState_t particle_state(nparticle_types,cavity_grid,particlepair_types,parameters.maximum_particle_number);


	/*------- Update cavity list of each particle type based on the solute -------*/
	if(parameters.solute_model==all_atom ){
		std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: Updating cavity lists based on solute atoms." << std::endl;
		particle_state.updateCavityGridVariables_SoluteAllAtom(particle_types,cavity_grid,solute);
	}

	/*------- Set the number and positions of particles for the initial state of the system -------*/

	int target_num[nparticle_types];

	double energy_old = 0;
	switch (parameters.start_state) {

	case 0:  /* Create initial state */
	{
		/*Set the number of ion types as the target number*/

		for(int i=0;i<niontype;i++){
			/* Set the target number of ions and solvent */
			target_num[i] = concMolarityToNum(iontypes.getIonConc(i),box.getBoxVol());
			std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: target number of particle type " << particle_types.getParticleType(i).label
					<< " = " << target_num[i] << ", chemical potential (kcal/mol) = " << particle_types.getParticleType(i).mu_old << std::endl;
		}
		if(nparticle_types-niontype==1){  // nsolvent = 1
			target_num[nparticle_types-1] = solvent_type.getTargetNumber(); //solvent_type.target_num;
			std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: target number of particle type " << particle_types.getParticleType(nparticle_types-1).label
					<< " = " << target_num[nparticle_types-1] <<
					", chemical potential (kcal/mol) = " << particle_types.getParticleType(nparticle_types-1).mu_old << std::endl;

		}
		std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1()]: SETTING UP INITIAL STATE OF SYSTEM USING GRID-BASED INSERTION---------->" << std::endl;

		//if(parameters.solute_model==all_atom ){

		//	particle_state.readInitializeParticleState_WithGridBasedInsertion(particle_types,cavity_grid,
		//	particlepair_types,0, 4, "./inputfiles/",target_num);


		//}else{
		particle_state.initializeParticleState_WithGridBasedInsertion(particle_types,cavity_grid,particlepair_types,target_num);
		//}

		std::cout
		<< "[RunGCMCSimulation1s.cpp (runSimulation_1)]: initial number of particles in the "
		<< "box (after placement and assignment of particle types = "
		<< particle_state.getParticleNum() << std::endl;

		if(parameters.use_pid == true && parameters.simulation_type == GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION){
			for (int j=0;j<nparticle_types;j++){
				particle_types.updateChemPotentialPID(j,box.getBoxVol(),parameters.temperature);
			}
		}/*else if(parameters.use_pid == true && parameters.simulation_type == GCMC_WITH_SOLUTE_ALLATOM_MODEL){
			for (int j=0; j< nparticle_types; j++){
			particle_types.targetChemPotentialPID(j,box.getBoxVol(),parameters.temperature);
			}
			}*/

		break;
	}
	case 1: /* Read initial state from file */
	{
		particle_state.readFromFile(0, 4,energy_old, "./inputfiles/");
		std::cout
		<< "[RunGCMCSimulation1.cpp (runSimulation_1)] Energy of the initial system = "
		<< energy_old << std::endl;


		std::cout
		<< "[RunGCMCSimulation1.cpp (runSimulation_1)] initial particle number read from file:"
		<< particle_state.getParticleNum() << std::endl;

		std::cout
		<< "[RunGCMCSimulation1.cpp (runSimulation_1)] initial energy of the system (kcal/mol:"
		<< energy_old << std::endl;

		/*Count the number of particles of each type in the simulation box */

		particle_state.countParticleTypeFromState(particle_types);


		/* Update cavity list and set mappings between particle indices and grid cell indices. */

		/*Update grid variables and number of particles of each type */
		particle_state.updateCavityGridVariables_AllParticles(particle_types,cavity_grid,particlepair_types);

		if(parameters.use_pid == true && parameters.simulation_type == GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION){
			for (int j=0; j< nparticle_types; j++){
				particle_types.updateChemPotentialPID(nparticle_types-1,box.getBoxVol(),parameters.temperature);
			}
		}else if(parameters.use_pid == true && parameters.simulation_type == GCMC_WITH_SOLUTE_ALLATOM_MODEL){
			for (int j=0; j< nparticle_types; j++){
				particle_types.targetChemPotentialPID(j,box.getBoxVol(),parameters.temperature);
			}
		}

		break;
	}
	default:
		std::cout
		<< "[RunGCMCSimulation1.cpp (runSimulation_1)]Error:  initial state of system not set to 0 or 1"
		<< std::endl;
		exit(0);

	}


	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1]: WRITING OUT INTIAL STATE OF PARTICLES---------->" << std::endl;
	particle_state.writeToFile(-1,-1,box,energy_old,"./outputfiles/state_");

	std::cout <<"[RunGCMCSimulation1.cpp(runSimulation_1]: FINISHED!" << std::endl;


	/*------- write out initial values of excess chemical potential,
	 * concentration, and number of particles of each type -------*/

	particle_types.writeChemicalPotentialToFile(-1);

	/*------- Set up the RDF data structures -------*/
	CRdf_t rdf(parameters,particle_types,particlepair_types,solute,box.getBoxVol());

	/*------- Start GCMC simulations for calculating excess chemical potential -------*/
	for (int ind1 = 0; ind1 < parameters.niter; ind1++) {
		gcmcCalc1(parameters,box,rdf,particle_types,particlepair_types,mwater,particle_state,cavity_grid,solute,energy_old);

		/*Update chemical potential and average number of each particle type, and
		 * write out the results
		 */

		if (parameters.simulation_type
				== GCMC_FOR_CHEMICAL_POTENTIAL_CALCULATION) {
			particle_types.updateChemicalPotential(ind1,parameters,box);
			particle_types.writeChemicalPotentialToFile(ind1);

			/* Replace the old chemical potential with the new one*/

			particle_types.replaceOldChemicalPotentialWithNew();
		}
		/*Write out the final state (x,y,z, radius, type label) of all the particles in
		 * the system after each iteration */
		particle_state.writeToFile(ind1,ind1,box,energy_old,"./outputfiles/finalstate_iter");

		particle_types.resetNumConcAvg();
		rdf.resetData();
	}
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
