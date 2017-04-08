/*
 * readStateFromFile.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file readStateFromFile.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief CParticleState class member function definitions for reading the state of the system
 *  		(coordinates, charge, radius, particle type label of all the particles) from a file.
 *
 */

#include "particle_state.hpp"

/**
 * @brief Reads the atom coordinates, radius, charge and atom type of particles from a file
 *
 * @param[in]	file_num(int)	Suffix number before the extension of the file name
 * @param[in]	num_lines_skip(int)	Number of lines to skip before starting to store the data
 * @param[out]	energy(double)		Energy of the system read from file in kcal/mol
 * @param[in]	dir_name(string)	Directory name containing the file to be read
 *
 * @return void
 */
void CParticleState::readFromFile(int file_num, int num_lines_skip, double &energy, std::string dir_name) {

	std::ifstream statefile;

	int i;

	int col_index;

	int string_count;
	int line_num;

	std::stringstream ss_suffix;
	std::string statefile_name;
	std::string suffix;

	std::string line;
	std::string* sub_strings;
	std::string sub_string;
	std::stringstream ss;

	struct_particle ptcl;
	const char* x = "";
	const char* y = "";
	const char* z = "";
	const char* radius = "";
	const char* charge = "";
	const char* chgcoeff = "";
	const char* ptype = "";

	sub_strings = NULL;

	ss_suffix << file_num;
	suffix = ss_suffix.str();

	statefile_name = dir_name;
	statefile_name.append("state_");
	statefile_name.append(suffix);

	line_num = 0;
	i = -1;

	statefile.open(statefile_name.c_str());

	if (statefile.is_open()) {
		while (getline(statefile, line)) {
			ss.clear();
			ss << "";
			ss << line;
			string_count = 0;

			while (ss >> sub_string) {
				++string_count;
			}
			// store the sub strings of a line into an array of strings
			sub_strings = new std::string[string_count];
			col_index = 0;
			ss.clear();
			ss << "";
			ss << line;
			while (ss >> sub_strings[col_index]) {
				++col_index;
			}

			++line_num ;

			if (line_num == 3){
				energy = atof(sub_strings[1].c_str());
			}
			if (line_num > num_lines_skip) {

				++i;

				if (string_count == 8) {
					x = sub_strings[1].c_str();
					y = sub_strings[2].c_str();
					z = sub_strings[3].c_str();
					radius = sub_strings[4].c_str();
					chgcoeff = sub_strings[5].c_str();
					ptype = sub_strings[6].c_str();

					ptcl.x = atof(x);
					ptcl.y = atof(y);
					ptcl.z = atof(z);
					ptcl.radius = atof(radius);
					ptcl.charge = atof(chgcoeff);
					ptcl.ptype = atoi(ptype);
					ptcl.label = sub_strings[7];

					particles[i]=ptcl;

				} else {
					std::cout
					<< "[readStateFromFile.cpp (readFromFile)] Error reading file: "
					<< statefile_name << ", at line number " << line_num
					<< std::endl;
					std::cout << "No. of columns is not equal to 8." << std::endl;
					exit(0);
				}
			}

			delete[] sub_strings;
			sub_strings = NULL;
		}

	} else {
		std::cout << "[readStateFromFile.cpp (readFromFile)] Error opening file: "
				<< statefile_name << std::endl;
		exit(0);
	}

	statefile.close();

	npart = i + 1; /* total number of particles (ions) in the state file */

	std::cout << "[readStateFromFile.cpp (readFromFile)] finished reading file: "
			<< statefile_name << std::endl;

}
