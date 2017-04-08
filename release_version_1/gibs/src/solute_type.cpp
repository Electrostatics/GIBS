/*
 * solute_type.cpp
 *
 *  Created on: Aug 12, 2015
 *      Author:  Dennis G. Thomas
 *
 *   @brief	CSoluteModel class constructor and member function definitions
 *
 *   @file solute_type.cpp
 *   @author Dennis G. Thomas
 *
 */

#include "solute_type.hpp"

/**
 * @brief	Constructs the solute model type
 *
 * Sets up solute molecule data structure with atom coordinates and potential
 * data information
 *
 * @param[in] parameters(CInputParameters_t) Simulation input parameters
 * @param[in] box(CBox_t)	Simulation box parameters
 *
 */
CSoluteModel::CSoluteModel(const CInputParameters_t &parameters,const CBox_t &box){
	this->solute_model = parameters.solute_model;

	switch (this->solute_model){
	case none:

		break;
	case all_atom:
		solute_allatom = new struct_solute_allatom;


		// set up coordinate data
		this->readSoluteAllAtomCoordsFromFile(parameters.solute_allatom_coordinate_filename,
				parameters.solute_allatom_coordinate_filetype);


		// set up potential grid data
		potential_grid = new struct_allatom_potentialgrid;
		potential_grid->nx = box.getNx();
		potential_grid->ny = box.getNy();
		potential_grid->nz = box.getNz();
		potential_grid->hx = box.getBoxHx();
		potential_grid->hy = box.getBoxHy();
		potential_grid->hz = box.getBoxHz();

		potential_grid->xgrid = new double[potential_grid->nx];
		potential_grid->ygrid = new double[potential_grid->ny];
		potential_grid->zgrid = new double[potential_grid->nz];

		for(int i=0;i<potential_grid->nx;i++){
			potential_grid->xgrid[i] = box.getBoxXMin()+double(i)*potential_grid->hx;
		}
		for(int i=0;i<potential_grid->ny;i++){
			potential_grid->ygrid[i] = box.getBoxYMin()+double(i)*potential_grid->hy;
		}
		for(int i=0;i<potential_grid->nz;i++){
			potential_grid->zgrid[i] = box.getBoxZMin() + double(i)*potential_grid->hz;
		}


		potential_grid->data = new double**[potential_grid->nx];

		for (int ind1 = 0; ind1 < potential_grid->nx; ind1++) {
			potential_grid->data[ind1] = new double*[potential_grid->ny];
			for (int ind2 = 0; ind2 < potential_grid->ny; ind2++) {
				potential_grid->data[ind1][ind2] = new double[potential_grid->nz];
			}
		}

		potential_grid->dxmap_filename = parameters.pot_dxmap_file;


		this->readSolutePotentialFromFile(box,parameters.temperature);

		break;

	case cylindrical_polyion:

		break;

	case spherical_macroion:

		break;

	default:
		std::cerr << "[solute_type.cpp(CSoluteModel::CSoluteModel)] : solute model type is not recognized." <<std::endl;

		break;
	}

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief	Reads the atom coordinates, radius, and charge values from a file
 * 			of an all-atom solute
 *
 *
 * @param[in] file_name(string) File name of the solute atom coordinate file
 * @param[in] file_type(enum_molecule_coordinate_file_type)	The type of solute coordinate file
 *
 * @return void
 */
void CSoluteModel::readSoluteAllAtomCoordsFromFile(std::string file_name,
		enum_molecule_coordinate_file_type file_type) {

	int string_count = 0;
	int index = 0;
	int atom_line_number = -1;

	std::string line;
	std::vector<std::string> sub_strings;
	std::string sub_string;
	std::stringstream ss;

	std::string filepath = "";
	std::ifstream molecule_file;

	enum_molecule_coordinate_file_type molecule_file_type;

	const char* fname = "";

	const char* field_name="";
	const char* atom_number = "";
	const char* atom_name = "";
	const char* residue_name = "";
	const char* chain_id = "";
	const char* residue_number = "";

	const char* x = "";
	const char* y = "";
	const char* z = "";
	const char* radius = "";
	const char* charge = "";
	const char* lj_epsilon = "";

	double total_charge = 0.0;

	molecule_file_type = file_type;

	filepath.append(file_name);
	fname = filepath.c_str();

	// molecule file name is taken as the molecule name
	solute_allatom->molecule_name = file_name;

	molecule_file.open(fname);

	if (molecule_file.is_open()) {
		std::cout << "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] File " << fname
				<< " is opened." << std::endl;

		while (getline(molecule_file, line)) {
			ss.clear();
			ss << "";
			ss << line;
			string_count = 0;

			while (ss >> sub_string) {
				++string_count;
			}

			// store the sub strings of a line into an array of strings
			sub_strings.resize(string_count);

			index = 0;
			ss.clear();
			ss << "";
			ss << line;
			while (ss >> sub_strings[index]) {
				++index;
			}

			switch (molecule_file_type) {

			case 0:
				/* XYZR file type
				 * no. of columns is 6.
				 *
				 * values read are x, y, z, radius, charge, and LJ epsilon values
				 */

				if (sub_strings[0] != "") {
					atom_line_number = atom_line_number + 1;

					x = sub_strings[0].c_str();
					y = sub_strings[1].c_str();
					z = sub_strings[2].c_str();
					radius = sub_strings[3].c_str();
					charge = sub_strings[4].c_str();
					lj_epsilon = sub_strings[5].c_str();

					solute_allatom->atoms[atom_line_number].x = atof(x);
					solute_allatom->atoms[atom_line_number].y = atof(y);
					solute_allatom->atoms[atom_line_number].z = atof(z);
					solute_allatom->atoms[atom_line_number].radius = atof(radius);
					solute_allatom->atoms[atom_line_number].charge = atof(charge);
					solute_allatom->atoms[atom_line_number].lj_epsilon = atof(
							lj_epsilon);

					total_charge += solute_allatom->atoms[atom_line_number].charge;

				}

				break;
			case 1:
				/* PQR file type
				 * no. of columns is either 10 or 11.
				 */

				if ((sub_strings[0] == "ATOM") && (string_count >= 10)) {
					atom_line_number = atom_line_number + 1;

					//	cout << "[read_molecule_file.cpp (readMoleculeFile)]: atom line number = " << atom_line_number << endl;
					field_name = sub_strings[0].c_str();
					atom_number = sub_strings[1].c_str();
					atom_name = sub_strings[2].c_str();
					residue_name = sub_strings[3].c_str();
					if (string_count == 11) {
						chain_id = sub_strings[4].c_str();
						residue_number = sub_strings[5].c_str();
						index = 5;
					} else if (string_count == 10) {
						residue_number = sub_strings[4].c_str();
						index = 4;
					} else {
						std::cerr
						<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)]: no. of columns in pqr file is neither 10 nor 11."
						<< std::endl;
						exit(0);
					}

					x = sub_strings[index + 1].c_str();
					y = sub_strings[index + 2].c_str();
					z = sub_strings[index + 3].c_str();
					charge = sub_strings[index + 4].c_str();
					radius = sub_strings[index + 5].c_str();

					solute_allatom->atoms[atom_line_number].atom_number = atoi(
							atom_number);
					solute_allatom->atoms[atom_line_number].atom_name = atom_name;
					solute_allatom->atoms[atom_line_number].residue_name =
							residue_name;
					solute_allatom->atoms[atom_line_number].chain_id = chain_id;
					solute_allatom->atoms[atom_line_number].residue_number = atoi(
							residue_number);

					solute_allatom->atoms[atom_line_number].x = atof(x);
					solute_allatom->atoms[atom_line_number].y = atof(y);
					solute_allatom->atoms[atom_line_number].z = atof(z);
					solute_allatom->atoms[atom_line_number].radius = atof(radius);
					solute_allatom->atoms[atom_line_number].charge = atof(charge);

					total_charge += solute_allatom->atoms[atom_line_number].charge;

					if (atom_line_number <= -1) {
						std::cout
						<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] atom line number = "
						<< atom_line_number << std::endl;
						std::cout
						<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] x position (A) = "
						<< solute_allatom->atoms[atom_line_number].x << std::endl;
						std::cout
						<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] y position (A) = "
						<< solute_allatom->atoms[atom_line_number].y << std::endl;
						std::cout
						<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] z position (A) = "
						<< solute_allatom->atoms[atom_line_number].z << std::endl;
						std::cout
						<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] radius (A) = "
						<< solute_allatom->atoms[atom_line_number].radius
						<< std::endl;
						std::cout
						<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] charge (e) = "
						<< solute_allatom->atoms[atom_line_number].charge
						<< std::endl;
						std::cout << " ----" << std::endl;

					}

				}

				break;

			default:
				std::cerr
				<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] Unrecognized solute molecule file type. "
				<< std::endl;
				exit(0);

			}

			sub_strings.clear();

		}
		solute_allatom->num_atoms = atom_line_number + 1;
		solute_allatom->total_charge = total_charge;

	} else {
		std::cerr << "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)]Error opening file: "
				<< fname << std::endl;
		exit(0);
	}

	molecule_file.close();


	std::cout
	<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] total number of atoms in molecule = "
	<< solute_allatom->num_atoms << "." << std::endl;
	std::cout
	<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] total charge on the molecule (e units) = "
	<< solute_allatom->total_charge << "." << std::endl;
	std::cout << "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)] File " << fname
			<< " is closed." << std::endl;


	/* find minimum and maximum values of x, y, and z in the molecule coordinate file */

	solute_allatom->min_x = solute_allatom->atoms[0].x;
	solute_allatom->max_x = solute_allatom->atoms[0].x;
	solute_allatom->min_y = solute_allatom->atoms[0].y;
	solute_allatom->max_y = solute_allatom->atoms[0].y;
	solute_allatom->min_z = solute_allatom->atoms[0].z;
	solute_allatom->max_z = solute_allatom->atoms[0].z;

	for (int ind1 = 1; ind1 < solute_allatom->num_atoms; ind1++) {
		if (solute_allatom->min_x >= solute_allatom->atoms[ind1].x) {
			solute_allatom->min_x = solute_allatom->atoms[ind1].x;
		}
		if (solute_allatom->max_x <= solute_allatom->atoms[ind1].x) {
			solute_allatom->max_x = solute_allatom->atoms[ind1].x;
		}
		if (solute_allatom->min_y >= solute_allatom->atoms[ind1].y) {
			solute_allatom->min_y = solute_allatom->atoms[ind1].y;
		}
		if (solute_allatom->max_y <= solute_allatom->atoms[ind1].y) {
			solute_allatom->max_y = solute_allatom->atoms[ind1].y;
		}

		if (solute_allatom->min_z >= solute_allatom->atoms[ind1].z) {
			solute_allatom->min_z = solute_allatom->atoms[ind1].z;
		}
		if (solute_allatom->max_z <= solute_allatom->atoms[ind1].z) {
			solute_allatom->max_z = solute_allatom->atoms[ind1].z;
		}

	}

	std::cout
	<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)]: lowest and highest values of an atom's position along x-axis (A): "
	<< solute_allatom->min_x << " and " << solute_allatom->max_x << std::endl;
	std::cout
	<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)]: lowest and highest values of an atom's position along y-axis (A): "
	<< solute_allatom->min_y << " and " << solute_allatom->max_y << std::endl;
	std::cout
	<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)]: lowest and highest values of an atom's position along z-axis (A): "
	<< solute_allatom->min_z << " and " << solute_allatom->max_z << std::endl;
	std::cout << "------" << std::endl;
	std::cout
	<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)]: maximum length of solute in x-axis direction (A): "
	<< solute_allatom->max_x - solute_allatom->min_x << std::endl;
	std::cout
	<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)]: maximum length of solute in y-axis direction (A): "
	<< solute_allatom->max_y - solute_allatom->min_y << std::endl;
	std::cout
	<< "[solute_type.cpp(CSoluteModel::readSoluteAllAtomCoordsFromFile)]: maximum length of solute in z-axis direction (A): "
	<< solute_allatom->max_z - solute_allatom->min_z << std::endl;
	std::cout << "------" << std::endl;


}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief	Reads the electrostatic potential from a .dx file
 *
 * The electrostatic potential file is in .DX format. The potential valueas are read
 * in kT/e units and converted to kcal/(mol e)
 *
 * @param[in] box(CBox_t) Simulation box parameters
 * @param[in] temp(double)	Temperature in Kelvin
 *
 * @return void
 */
void CSoluteModel::readSolutePotentialFromFile(const CBox_t &box,double temp){

	std::ifstream readdxfile;
	const char* fname = "";

	int nrows_skip = 11;
	int ncols = 3;


	std::string sline;
	std::string filename = "";
	filename.append(potential_grid->dxmap_filename);

	fname = potential_grid->dxmap_filename.c_str();
	readdxfile.open(fname);

	double **preadvalue;

	int npoints = potential_grid->nx*potential_grid->ny*potential_grid->nz;
	int nrows;

	if(npoints % ncols == 0)
		nrows = npoints/ncols;
	else
		nrows = npoints/ncols+1;

	std::cout << "[solute_type.cpp (CSoluteModel::readSolutePotentialFromFile)]: "
			"no. of lines read from file " << potential_grid->dxmap_filename << " = " << nrows << std::endl;

	preadvalue = new double*[nrows];
	for(int ind=0;ind<nrows;ind++){
		preadvalue[ind] = new double[ncols];
	}


	if(readdxfile.is_open()){

		for(int ind=0;ind<nrows_skip;ind++) {
			getline(readdxfile,sline);

		}

		if(npoints % ncols ==0){
			for(int ind=0;ind<nrows;ind++){

				readdxfile >> preadvalue[ind][0] >> preadvalue[ind][1] >> preadvalue[ind][2];
			}

		}
		else {
			for(int ind=0;ind<nrows-1;ind++){
				readdxfile >> preadvalue[ind][0] >> preadvalue[ind][1] >> preadvalue[ind][2];
			}

			readdxfile >> preadvalue[nrows-1][0] >> preadvalue[nrows-1][1];
		}

	}

	else {
		std::cout << "[solute_type.cpp (CSoluteModel::readSolutePotentialFromFile)]:"
				" Error opening file: " << fname << std::endl;
	}
	readdxfile.close();

	for(int ix=0;ix<potential_grid->nx;ix++){
		for(int iy=0;iy<potential_grid->ny;iy++){
			for(int iz=0;iz<potential_grid->nz;iz++){

				int count = ix*potential_grid->ny*potential_grid->nz+iy*potential_grid->nz+iz+1;
				int ir,ic;
				if(count % ncols == 0) {
					ir = count/ncols;
				}
				else {
					ir=count/ncols+1;
				}
				ic = count-(ir-1)*ncols;
				potential_grid->data[ix][iy][iz] = preadvalue[ir-1][ic-1];

			}
		}
	}

	for (int ind=0;ind<nrows;ind++){
		delete[] preadvalue[ind];
		preadvalue[ind] = NULL;
	}

	delete[] preadvalue;
	preadvalue = NULL;



	/* convert units of potential from kT/e to kcal/(mol e) */
	for (int ind1 = 0; ind1 < potential_grid->nx; ind1++) {
		for (int ind2 = 0; ind2 < potential_grid->ny; ind2++) {
			for (int ind3 = 0; ind3 < potential_grid->nz; ind3++) {

				/* conversion factor = k_B*avognum*convJtoCal*0.001*temperature = 0.001987188 * temperature */

				potential_grid->data[ind1][ind2][ind3] *= 0.001987188 * temp;
			}
		}
	}
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief	Interpolates the electrostatic potential from potential grid data
 * 			at any location in the simulation box
 *
 * The potential is obtained using trilinear interpolation method
 *
 * @param[in] x(double) x coordinate
 * @param[in] y(double) y coordinate
 * @param[in] z(double) z coordinate
 * @param[in] box(CBox_t) Simulation box parameters
 *
 * @return Interpolated value of the electrostatic potential(double) at (x,y,z) in kcal/(mol e) units
 *
 */
double CSoluteModel::getSoluteAllAtom_Potential(double x, double y, double z,const CBox_t box) const{

	double c0, c1, c2, c3, c4, c5, c6, c7;
	double p000, p100,p010,p110,p001, p011, p101,p111;



	double	dx0 = floor((x - box.getBoxXMin())/box.getBoxHx());
	double dy0 = floor((y - box.getBoxYMin())/box.getBoxHy());
	double	dz0 = floor((z - box.getBoxZMin())/box.getBoxHz());

	int x0 = (int)dx0;
	int y0 = (int)dy0;
	int z0 = (int)dz0;

	double delx = (x-potential_grid->xgrid[x0])/potential_grid->hx;
	double dely = (y-potential_grid->ygrid[y0])/potential_grid->hy;
	double delz = (z-potential_grid->zgrid[z0])/potential_grid->hz;

	p000 = potential_grid->data[x0][y0][z0];


	//
	if(delx!=0.0){
		p100 = potential_grid->data[x0+1][y0][z0];

	}
	else{
		p100 = p000;
	}
	//
	if(dely!=0.0){
		p010 = potential_grid->data[x0][y0+1][z0];
	}
	else{
		p010 = p000;

	}
	//
	if(delz!=0.0){
		p001 = potential_grid->data[x0][y0][z0+1];
	}
	else{
		p001 = p000;

	}
	//
	if((delx!=0.0) && (dely!=0.0)){
		p110 = potential_grid->data[x0+1][y0+1][z0];
	}
	else if((delx==0.0) && (dely!=0.0)){
		p110 = p010;
	}
	else if((delx!=0.0) && (dely==0.0)){
		p110 = p100;
	}

	//
	if((dely!=0.0)&&(delz!=0.0)){
		p011 = potential_grid->data[x0][y0+1][z0+1];
	}
	else if((dely==0.0) &&(delz!=0.0)){
		p011 = p001;
	}
	else if((dely!=0.0)&&(delz==0.0)){
		p011 = p010;
	}


	//
	if((delx!=0.0)&&(delz!=0.0)){
		p101 = potential_grid->data[x0+1][y0][z0+1];
	}
	else if((delx==0.0) &&(delz!=0.0)){
		p101 = p001;
	}
	else if((delx!=0.0)&&(delz==0.0)){
		p101 = p100;
	}


	//
	if((delx!=0.0) &&(dely!=0.0) && (delz!=0.0)){
		p111 = potential_grid->data [x0+1][y0+1][z0+1];
	}
	else if((delx==0.0) &&(dely!=0.0) && (delz!=0.0)){
		p111 = p011;
	}
	else if((delx!=0.0) &&(dely!=0.0) && (delz==0.0)){
		p111 = p110;
	}
	else if((delx!=0.0) &&(dely==0.0) && (delz!=0.0)){
		p111 = p101;
	}
	else if((delx==0.0) &&(dely==0.0) && (delz!=0.0)){
		p111 = p001;
	}
	else if((delx==0.0) &&(dely!=0.0) && (delz==0.0)){
		p111 = p010;
	}
	else if((delx!=0.0) &&(dely==0.0) && (delz==0.0)){
		p111 = p100;
	}


	c0 = p000;
	c1 = p100 - p000;
	c2 = p010 - p000;
	c3 = p001 - p000;
	c4 = p110 - p010 - p100 + p000;
	c5 = p011 - p001 - p010 + p000;
	c6 = p101 - p001 - p100 + p000;
	c7 = p111 - p011 - p101 - p110 + p100 + p001 + p010 - p000;

	double value = c0 + c1*delx + c2*dely + c3*delz + c4*delx*dely
			+ c5*dely*delz + c6*delz*delx + c7*delx*dely*delz;

	return value;
}
