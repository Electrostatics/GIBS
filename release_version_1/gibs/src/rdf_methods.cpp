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
 * rdf_methods.cpp
 *
 *  Created on: Jul 9, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file rdf_methods.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Member function definitions of class CRdf for computing radial distribution functions
 *
 *
 */

#include "rdf_methods.hpp"

/**
 * @brief Constructor of class CRdf
 *
 * @param[in]	parameters(CInputParameters_t)
 * @param[in]	particle_types(CParticleType_t)
 * @param[in]	particlepair_types(CParticlePairType_t)
 * @param[in]	solute_model(CSoluteModel_t)
 * @param[in]	box_vol(double)
 *
 */
CRdf::CRdf(const CInputParameters_t &parameters, const CParticleType_t &particle_types,
		const CParticlePairType_t &particlepair_types,const CSoluteModel_t &solute_model,double box_vol){
	this->nparticle_types = particle_types.getNumParticleTypes();
	this->nparticle_pairs = particlepair_types.getNumPairs();

	nrdfbin = parameters.nrdfbin;
	delr_rdf = parameters.delr_rdf;
	xlen = parameters.xlen;
	ylen = parameters.ylen;
	zlen = parameters.zlen;
	xmin = parameters.xmin;
	ymin = parameters.ymin;
	zmin = parameters.zmin;
	this->box_vol = box_vol;

	// set up particle pair RDF
	particlepair_rdf= new struct_particlepair_rdf[nparticle_pairs];

	for (int i=0;i<nparticle_pairs;i++){

		particlepair_rdf[i].particlepair_label = particlepair_types.getLabelFromIndex(i);
		particlepair_rdf[i].ref_ptcl_type = particlepair_types.getTypeIndex1(i);
		particlepair_rdf[i].rdf_ptcl_type = particlepair_types.getTypeIndex2(i);

		std::cout <<"[rdf_methods.cpp (setUpParticlePairRDF]: particle pair labels:" <<
				particlepair_rdf[i].particlepair_label << std::endl;

		// get mid point radial distance value of the first bin

		particlepair_rdf[i].rmid_firstbin = particle_types.getRadius(particlepair_rdf[i].ref_ptcl_type-1) +
				particle_types.getRadius(particlepair_rdf[i].rdf_ptcl_type-1);

		std::cout << "[rdf.h (class CRdf)]: allocating memory for rdf data of particle pair type, " <<
				particlepair_rdf[i].particlepair_label << std::endl;

		particlepair_rdf[i].data = new struct_rdf_data;
		particlepair_rdf[i].rdist = new double[nrdfbin];


		particlepair_rdf[i].data->grdf = new double[nrdfbin];
		particlepair_rdf[i].data->nbin = new double[nrdfbin];
		particlepair_rdf[i].data->nbin_avg = new double[nrdfbin];
		particlepair_rdf[i].data->nbin_avg_sqr = new double[nrdfbin];
		particlepair_rdf[i].data->ndist = new double[nrdfbin];
		particlepair_rdf[i].data->ndist_avg = new double[nrdfbin];
		particlepair_rdf[i].data->ndist_avg_sqr = new double[nrdfbin];


		particlepair_rdf[i].data->nbin_prev = new double[nrdfbin];
		particlepair_rdf[i].data->ndist_prev = new double[nrdfbin];

		/* initialize values*/
		particlepair_rdf[i].ndist_icount = 0;
		for (int j = 0; j < nrdfbin; j++) {
			particlepair_rdf[i].data->ndist[j] = 0;
			particlepair_rdf[i].data->ndist_avg[j] = 0.0;
			particlepair_rdf[i].data->ndist_avg_sqr[j] = 0.0;
			particlepair_rdf[i].data->nbin[j] = 0;
			particlepair_rdf[i].data->nbin_avg[j] = 0.0;
			particlepair_rdf[i].data->nbin_avg_sqr[j] = 0.0;

			particlepair_rdf[i].data->ndist_prev[j] = 0;
			particlepair_rdf[i].data->nbin_prev[j] = 0;

			particlepair_rdf[i].data->grdf[j] = 0.0;

			particlepair_rdf[i].rdist[j] = particlepair_rdf[i].rmid_firstbin + j * delr_rdf;
		}

	}

	this->calculate_particlepair_rdf = parameters.calculate_particle_pair_rdf;
	// declare all atom solute particle rdf

	this->solute_model_type =parameters.solute_model;
	this->calculate_solute_particle_rdf = parameters.calculate_solute_particle_rdf;
	this->solute_total_charge = 0.0;
	if(solute_model_type==all_atom && this->calculate_solute_particle_rdf ==true){
		this->soluteparticle_rdfref = parameters.rdf_ref;
		this->solute_total_charge = solute_model.getSoluteAllAtom_TotalCharge();
		soluteallatom_particle_rdf = new struct_soluteallatom_particle_rdf[nparticle_types];

		this->setUpSoluteAllAtomParticleRdf(parameters,particle_types,solute_model);

	}
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Calculates particle pair RDFs
 *
 * @param[in]	particle_types(CParticleType_t)	Particle type parameters and cavity grid data
 *
 * @return void
 */
void CRdf::calcParticlePairRDF(const CParticleType_t &particle_types) {



	double bin_volume;

	for (int ind1 = 0; ind1 < nparticle_pairs; ind1++) {

		int ptype = particlepair_rdf[ind1].rdf_ptcl_type;
		int ref_ptype = particlepair_rdf[ind1].ref_ptcl_type;

		double ptcl_num_density = concMolarityToNum(particle_types.getParticleType(ptype - 1).conc, box_vol)
																																																																/ box_vol; //  ion number density used to normalize the rdf
		for (int ind2 = 0; ind2 < nrdfbin; ind2++) {

			if (ind2 == 0) {
				bin_volume = 4.0 * PI
						* (pow(particlepair_rdf[ind1].rdist[ind2] + 0.5 * delr_rdf, 3)
								- pow(particlepair_rdf[ind1].rdist[ind2], 3)) / 3.0;
			} else {
				bin_volume = 4.0 * PI
						* (pow(particlepair_rdf[ind1].rdist[ind2] + 0.5 * delr_rdf, 3)
								- pow(particlepair_rdf[ind1].rdist[ind2] - 0.5 * delr_rdf,
										3)) / 3.0;
			}


			particlepair_rdf[ind1].data->grdf[ind2] = particlepair_rdf[ind1].data->nbin_avg[ind2]
			                                                                                /(bin_volume * particle_types.getParticleType(ref_ptype - 1).num * ptcl_num_density);
		}
	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Counts the number of each particle type within a radial distance from a selected particle
 *
 * @param[in]	state_change(enum_state_change)	Type of state change
 * @param[in]	particle_types(CParticleType_t)	Particle type parameters and cavity grid data
 * @param[in]	particlepair_types(CParticlePairType_t)	Particle type parameters
 * @param[in]	particle_state(CParticleState_t)	Particle state data
 * @param[in]	ptcl(struct_particle)	Struct data of selected particle
 * @param[in]	ptclindex(int)		Array/vector index of the selected particle in the simulation box
 * @param[in]	ptcltype_index(int)	Array/vevtor index of the selected particle type
 *
 * @return void
 */

void CRdf::calcParticlePairBinCountsPerStateChange(enum_state_change state_change,const CParticleType_t &particle_types,
		const CParticlePairType_t &particlepair_types,const CParticleState_t &particle_state,
		const struct_particle &ptcl,int ptclindex, int ptcltype_index){

	int nparticle_types = particle_types.getNumParticleTypes();

	for (int ind1=0;ind1<nparticle_pairs;ind1++){
		for(int ind2=0; ind2<nrdfbin; ind2++){
			particlepair_rdf[ind1].data->ndist[ind2] = particlepair_rdf[ind1].data->ndist_prev[ind2];
			particlepair_rdf[ind1].data->nbin[ind2] = particlepair_rdf[ind1].data->nbin_prev[ind2] ;
		}
	}

	for (int ind1=0; ind1< nparticle_types;ind1++){
		std::string label12 = ptcl.label;

		label12.append("_");
		label12.append(particle_types.getLabel(ind1));
		int particlepair_index = particlepair_types.getIndexFromLabel(label12);



		/*for (int ind2 = 0; ind2 < nrdfbin; ind2++) {
			particlepair_rdf[particlepair_index].data->ndist[ind2] = particlepair_rdf[particlepair_index].data->ndist_prev[ind2];
			particlepair_rdf[particlepair_index].data->nbin[ind2] = particlepair_rdf[particlepair_index].data->nbin_prev[ind2] ;
		}*/

		int ptype = particlepair_rdf[particlepair_index].rdf_ptcl_type;
		int ref_ptype = particlepair_rdf[particlepair_index].ref_ptcl_type;

		int reftype_num = particle_types.getNum(ref_ptype-1);
		int ptype_num = particle_types.getNum(ptype-1);

		double val = 1.0*RDF_SCALING_FACTOR;
		if (ref_ptype == ptype){
			val = 2.0*RDF_SCALING_FACTOR;
		}
		particlepair_rdf[particlepair_index].ref_ptcl_count = reftype_num;

		int loop_ptcltype_index = ind1;
		int loop_ptcltype_num = particle_types.getNum(ind1);

		/*		// The particle type of ptcl - is it the reftype?
		if(ref_ptype == ptcl.ptype){
			// then loop over particles of type ptype
			loop_ptcltype_index = ptype - 1;
			loop_ptcltype_num = ptype_num;
		}else{
			loop_ptcltype_index = ref_ptype - 1;
			loop_ptcltype_num = reftype_num;
		}*/

		double	xi = ptcl.x;
		double	yi = ptcl.y;
		double	zi = ptcl.z;

		for (int ind4 = 0; ind4 < loop_ptcltype_num; ind4++) {


			int ptclindex_j = particle_state.getParticleIndicesByCount(loop_ptcltype_index,ind4);
			if(ptclindex != ptclindex_j){
				struct_particle ptclj = particle_state.getParticle(ptclindex_j);

				double	xj = ptclj.x;
				double	yj = ptclj.y;
				double	zj = ptclj.z;

				// Find the nearest image
				if (xj < xi - xlen * 0.5) {
					xj = xj + xlen;
				} else if (xj > xi + xlen * 0.5) {
					xj = xj - xlen;
				}
				if (yj < yi - ylen * 0.5) {
					yj = yj + ylen;
				} else if (yj > yi + ylen * 0.5) {
					yj = yj - ylen;
				}
				if (zj < zi - zlen * 0.5) {
					zj = zj + zlen;
				} else if (zj > zi + zlen * 0.5) {
					zj = zj - zlen;
				}

				double xij = xi - xj;
				double yij = yi - yj;
				double zij = zi - zj;

				double	rij = sqrt(xij * xij + yij * yij + zij * zij);

				int pos_bin = int(
						(rij - (ptcl.radius + ptclj.radius - 0.5*delr_rdf)) / delr_rdf);
				if(pos_bin < nrdfbin){
					if(state_change==insertion){

						particlepair_rdf[particlepair_index].data->nbin[pos_bin] += val;



					}else if(state_change==deletion){

						particlepair_rdf[particlepair_index].data->nbin[pos_bin] -= val;
					}
				}

			}
		}


		/*particlepair_rdf[particlepair_index].data->nbin[0] = particlepair_rdf[particlepair_index].data->ndist[0];

				for (int pos_bin = 1; pos_bin < nrdfbin; pos_bin++) {
					particlepair_rdf[particlepair_index].data->nbin[pos_bin] =
							particlepair_rdf[particlepair_index].data->ndist[pos_bin]
							                                   - particlepair_rdf[particlepair_index].data->ndist[pos_bin-1];
				}*/
		//particlepair_rdf[particlepair_index].data->ndist_prev= particlepair_rdf[particlepair_index].data->ndist;
		//	particlepair_rdf[particlepair_index].data->nbin_prev = particlepair_rdf[particlepair_index].data->nbin;



		for (int pos_bin = 0; pos_bin < nrdfbin; pos_bin++) {
			if(pos_bin==0){
				particlepair_rdf[ind1].data->ndist[pos_bin] = particlepair_rdf[ind1].data->nbin[pos_bin];
			}else{
				particlepair_rdf[ind1].data->ndist[pos_bin] = particlepair_rdf[ind1].data->ndist[pos_bin-1]
				                                                                                 + particlepair_rdf[ind1].data->nbin[pos_bin];
			}

			particlepair_rdf[particlepair_index].data->nbin_prev[pos_bin]= particlepair_rdf[particlepair_index].data->nbin[pos_bin];
			particlepair_rdf[particlepair_index].data->ndist_prev[pos_bin]= particlepair_rdf[particlepair_index].data->ndist[pos_bin];
		}

	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Calculate the radial distribution of the number of all particle pair types in the simulation box
 *
 * @param[in]	particle_types(CParticleType_t)	Particle type parameters and cavity grid data
 * @param[in]	particle_state(CParticleState_t)	Particle state data
 *
 * @return void
 */
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void CRdf::calcParticlePairBinCounts(const CParticleType_t &particle_types,const CParticleState_t &particle_state) {

	int ind4_start = 0;
	int ind3_end_offset = 0;
	for (int ind1 = 0; ind1 < nparticle_pairs; ind1++) {

		for (int ind2 = 0; ind2 < nrdfbin; ind2++) {
			particlepair_rdf[ind1].data->ndist[ind2] = 0;
			particlepair_rdf[ind1].data->nbin[ind2] = 0;
		}

		int ptype = particlepair_rdf[ind1].rdf_ptcl_type;
		int ref_ptype = particlepair_rdf[ind1].ref_ptcl_type;

		int reftype_num = particle_types.getNum(ref_ptype-1);
		int ptype_num = particle_types.getNum(ptype-1);
		particlepair_rdf[ind1].ref_ptcl_count = reftype_num;

		double val = 1.0;
		if(ptype == ref_ptype){
			ind3_end_offset = -1;
			val = 2.0*RDF_SCALING_FACTOR;
		} else{
			ind3_end_offset = 0;
			val = 1.0*RDF_SCALING_FACTOR;
		}


		for (int ind3 = 0; ind3 < reftype_num+ind3_end_offset; ind3++) {
			int ptclindex_i = particle_state.getParticleIndicesByCount(ref_ptype-1,ind3);
			struct_particle ptcli = particle_state.getParticle(ptclindex_i);
			if(ptype == ref_ptype){
				ind4_start = ind3 + 1;
			} else{
				ind4_start = 0;
			}
			for (int ind4 = ind4_start; ind4 < ptype_num; ind4++) {
				int ptclindex_j = particle_state.getParticleIndicesByCount(ptype-1,ind4);

				struct_particle ptclj = particle_state.getParticle(ptclindex_j);

				double	xi = ptcli.x;
				double	yi = ptcli.y;
				double	zi = ptcli.z;

				double	xj = ptclj.x;
				double	yj = ptclj.y;
				double	zj = ptclj.z;

				// Find the nearest image
				if (xj < xi - xlen * 0.5) {
					xj = xj + xlen;
				} else if (xj > xi + xlen * 0.5) {
					xj = xj - xlen;
				}
				if (yj < yi - ylen * 0.5) {
					yj = yj + ylen;
				} else if (yj > yi + ylen * 0.5) {
					yj = yj - ylen;
				}
				if (zj < zi - zlen * 0.5) {
					zj = zj + zlen;
				} else if (zj > zi + zlen * 0.5) {
					zj = zj - zlen;
				}

				double xij = xi - xj;
				double yij = yi - yj;
				double zij = zi - zj;

				double	rij = sqrt(xij * xij + yij * yij + zij * zij);
				//int pos_bin=int(rij/delr_rdf);
				int pos_bin = int(
						(rij - (ptcli.radius + ptclj.radius- 0.5*delr_rdf)) / delr_rdf);

				if (pos_bin < nrdfbin){
					particlepair_rdf[ind1].data->nbin[pos_bin] += val; // val = 2 if the two particles are of the same type
				}

			}

		}
		for (int pos_bin = 0; pos_bin < nrdfbin; pos_bin++) {

			if(pos_bin==0){
				particlepair_rdf[ind1].data->ndist[pos_bin] = particlepair_rdf[ind1].data->nbin[pos_bin];
			}else{
				particlepair_rdf[ind1].data->ndist[pos_bin] = particlepair_rdf[ind1].data->ndist[pos_bin-1]
				                                                                                 + particlepair_rdf[ind1].data->nbin[pos_bin];
			}
			particlepair_rdf[ind1].data->nbin_prev[pos_bin]= particlepair_rdf[ind1].data->nbin[pos_bin];
			particlepair_rdf[ind1].data->ndist_prev[pos_bin]= particlepair_rdf[ind1].data->ndist[pos_bin];
		}
		//particlepair_rdf[ind1].data->nbin_prev = particlepair_rdf[ind1].data->nbin;
	}

	/* calculate number of ions in each bin for each segment */
	/*	for (int ind2 = 0; ind2 < nparticle_pairs; ind2++) {

		particlepair_rdf[ind2].data->nbin[0] = particlepair_rdf[ind2].data->ndist[0];

	for (int pos_bin = 1; pos_bin < nrdfbin; pos_bin++) {
			particlepair_rdf[ind2].data->nbin[pos_bin] =
					particlepair_rdf[ind2].data->ndist[pos_bin]
					                                   - particlepair_rdf[ind2].data->ndist[pos_bin - 1];
		}

		particlepair_rdf[ind2].data->ndist_prev= particlepair_rdf[ind2].data->ndist;
		particlepair_rdf[ind2].data->nbin_prev = particlepair_rdf[ind2].data->nbin;


	}*/
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Increment the total bin counts after sampling each state of the system for computing
 * 			an average value at the end
 *
 * @return void
 *
 */
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

void CRdf::incrementParticlePairBinAvgCounts(){

	for (int ind1 = 0; ind1 < nparticle_pairs; ind1++) {


		for (int ind2 = 0; ind2 < nrdfbin; ind2++) {
			particlepair_rdf[ind1].data->ndist_avg[ind2] += particlepair_rdf[ind1].data->ndist[ind2];
			particlepair_rdf[ind1].data->ndist_avg_sqr[ind2] += particlepair_rdf[ind1].data->ndist[ind2]*particlepair_rdf[ind1].data->ndist[ind2];


			particlepair_rdf[ind1].data->nbin_avg[ind2] +=
					particlepair_rdf[ind1].data->nbin[ind2];

			particlepair_rdf[ind1].data->nbin_avg_sqr[ind2] +=
					particlepair_rdf[ind1].data->nbin[ind2]*particlepair_rdf[ind1].data->nbin[ind2];
		}


	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Compute the average bin counts for each particle pair type
 *
 * @param[in]	state_count(int)	Number of states sampled for compute the RDFs
 *
 * @return void
 *
 */
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

void CRdf::finalizeParticlePairBinAvgCounts(int state_count){

	double inv_rdf_scaling_factor = 1.0/RDF_SCALING_FACTOR;

	for (int ind1 = 0; ind1 < nparticle_pairs; ind1++) {

		particlepair_rdf[ind1].ndist_icount = state_count;
		for (int ind2 = 0; ind2 < nrdfbin; ind2++) {
			particlepair_rdf[ind1].data->ndist_avg[ind2] =
					particlepair_rdf[ind1].data->ndist_avg[ind2]*inv_rdf_scaling_factor/double(state_count);

			particlepair_rdf[ind1].data->ndist_avg_sqr[ind2] =
					particlepair_rdf[ind1].data->ndist_avg_sqr[ind2]*inv_rdf_scaling_factor*inv_rdf_scaling_factor/double(state_count);

			particlepair_rdf[ind1].data->nbin_avg[ind2] =
					particlepair_rdf[ind1].data->nbin_avg[ind2]*inv_rdf_scaling_factor/double(state_count);

			particlepair_rdf[ind1].data->nbin_avg_sqr[ind2] =	particlepair_rdf[ind1].data->nbin_avg_sqr[ind2]
			                                                 	                                          *inv_rdf_scaling_factor*inv_rdf_scaling_factor
			                                                 	                                          /double(state_count);
		}
	}
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief	Write each particle pair RDF to a file
 *
 * The output file name of each particle pair type has the format
 * rdfpair_{particle pair type label}.out
 *
 * @return void
 *
 */
void CRdf::writeParticlePairRDFtoFile() {
	std::ofstream writefile;

	const char* fname = "";

	for (int ind = 0; ind < nparticle_pairs; ind++) {

		std::string filename = "./outputfiles/rdfpair_";
		filename.append(particlepair_rdf[ind].particlepair_label);
		filename.append(".out");

		fname = filename.c_str();
		//writefile.open(fname, std::ios::app);
		writefile.open(fname);

		/*writefile << "r (A) of bin mid-point" << '\t' << "rho/rho_bulk"
				<< '\t' << "cumulative sum (avg)" << '\t'
				<< "cumulative sum(s.d.)" << '\t' << "# of ions in bin (avg)"
				<< '\t' << "# ions in bin (s.d.)" << std::endl;
		 */


		writefile << "r_of_bin_mid-point_in_A" << '\t' << "rdf"
				<< '\t' << "avg_cumulative_sum" << '\t'
				<< "sd_cumulative_sum" << '\t' << "avg_particle_count_in_bin"
				<< '\t' << "sd_particle_count_in_bin" << std::endl;


		for (int ind1 = 0; ind1 < nrdfbin; ind1++) {
			writefile << particlepair_rdf[ind].rdist[ind1] << '\t'
					<< particlepair_rdf[ind].data->grdf[ind1] << '\t' <<
					particlepair_rdf[ind].data->ndist_avg[ind1] <<	'\t' <<
					sqrt(particlepair_rdf[ind].data->ndist_avg_sqr[ind1]-
							particlepair_rdf[ind].data->ndist_avg[ind1]* particlepair_rdf[ind].data->ndist_avg[ind1])
							<< '\t' << particlepair_rdf[ind].data->nbin_avg[ind1] << '\t'
							<< sqrt(particlepair_rdf[ind].data->nbin_avg_sqr[ind1] -
									particlepair_rdf[ind].data->nbin_avg[ind1]* particlepair_rdf[ind].data->nbin_avg[ind1])
									<< std::endl;

			/*writefile << particlepair_rdf[ind].rdist[ind1] << '\t'
					<< particlepair_rdf[ind].data->grdf[ind1] << std::endl;*/

		}

		writefile.close();
	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief	Reset the bin counts and RDF values to 0
 *
 * @return void
 */

void CRdf::resetData(){
	for (int i=0;i<nparticle_pairs;i++){
		particlepair_rdf[i].ndist_icount = 0;
		for (int j = 0; j < nrdfbin; j++) {
			particlepair_rdf[i].data->ndist[j] = 0;
			particlepair_rdf[i].data->ndist_avg[j] = 0.0;
			particlepair_rdf[i].data->ndist_avg_sqr[j] = 0.0;
			particlepair_rdf[i].data->nbin[j] = 0;
			particlepair_rdf[i].data->nbin_avg[j] = 0.0;
			particlepair_rdf[i].data->nbin_avg_sqr[j] = 0.0;

			particlepair_rdf[i].data->grdf[j] = 0.0;

			particlepair_rdf[i].data->ndist_prev[j] = 0;
			particlepair_rdf[i].data->nbin_prev[j] = 0;

		}
	}
	if(solute_model_type==all_atom && this->calculate_solute_particle_rdf==true){
		for (int i = 0; i < nparticle_types;i++){
			for (int j =0; j < nrdfbin; j++){

				soluteallatom_particle_rdf[i].data->ndist[j] = 0;
				soluteallatom_particle_rdf[i].data->ndist_avg[j] = 0.0;
				soluteallatom_particle_rdf[i].data->ndist_avg_sqr[j] = 0.0;
				soluteallatom_particle_rdf[i].data->nbin[j] = 0;
				soluteallatom_particle_rdf[i].data->nbin_avg[j] = 0.0;
				soluteallatom_particle_rdf[i].data->nbin_avg_sqr[j] = 0.0;

				soluteallatom_particle_rdf[i].data->grdf[j] = 0.0;

				soluteallatom_particle_rdf[i].data->ndist_prev[j] = 0;
				soluteallatom_particle_rdf[i].data->nbin_prev[j] = 0;
			}
		}


	}


}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Initialize and set up the arrays and parameters for computing the radial distribution of
 * 		each particle type around the solute
 *
 * @param[in]	parameters(CInputParameters_t)	Input parameters specified by the user
 * @param[in]	particle_types(CParticleType_t)	Parameters and data of all particle types
 * @param[in]	solute_model(CSoluteModel_t)	Data and parameters of solute model
 *
 * @return void
 *
 */
void CRdf::setUpSoluteAllAtomParticleRdf(const CInputParameters_t &parameters,
		const CParticleType_t &particle_types,const CSoluteModel_t &solute_model){

	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {


		soluteallatom_particle_rdf[ind1].data = new struct_rdf_data;
		soluteallatom_particle_rdf[ind1].rdist = new double[nrdfbin];
		/* allocate memory for the members of rdf data struct */

		soluteallatom_particle_rdf[ind1].data->grdf = new double[nrdfbin];
		soluteallatom_particle_rdf[ind1].data->nbin = new double[nrdfbin];
		soluteallatom_particle_rdf[ind1].data->nbin_avg = new double[nrdfbin];
		soluteallatom_particle_rdf[ind1].data->nbin_avg_sqr = new double[nrdfbin];
		soluteallatom_particle_rdf[ind1].data->ndist = new double[nrdfbin];
		soluteallatom_particle_rdf[ind1].data->ndist_avg = new double[nrdfbin];
		soluteallatom_particle_rdf[ind1].data->ndist_avg_sqr = new double[nrdfbin];

		soluteallatom_particle_rdf[ind1].data->ndist_prev = new double[nrdfbin];
		soluteallatom_particle_rdf[ind1].data->nbin_prev = new double[nrdfbin];

		/* initialize values*/

		/* specify the axis from which the rdfs have to be computed */
		soluteallatom_particle_rdf[ind1].cylindrical_rdf_axis_dir =
				parameters.cylindrical_rdf_axis_direction;

		if (parameters.rdf_ref == AXIS_OF_CYLINDER) {

			soluteallatom_particle_rdf[ind1].ref_xcoord = solute_model.getSoluteAllAtom_MinX()
																																											+ 0.5 * (solute_model.getSoluteAllAtom_MaxX() - solute_model.getSoluteAllAtom_MinX());
			soluteallatom_particle_rdf[ind1].ref_ycoord = solute_model.getSoluteAllAtom_MinY()
																																											+ 0.5 * (solute_model.getSoluteAllAtom_MaxY() - solute_model.getSoluteAllAtom_MinY());
			soluteallatom_particle_rdf[ind1].ref_zcoord = solute_model.getSoluteAllAtom_MinZ()
																																											+ 0.5 * (solute_model.getSoluteAllAtom_MaxZ() - solute_model.getSoluteAllAtom_MinZ());

		} else if (parameters.rdf_ref == CENTER_OF_SPHERE) {
			soluteallatom_particle_rdf[ind1].ref_xcoord = solute_model.getSoluteAllAtom_AtomCoordX(0);
			soluteallatom_particle_rdf[ind1].ref_ycoord = solute_model.getSoluteAllAtom_AtomCoordY(0);
			soluteallatom_particle_rdf[ind1].ref_zcoord = solute_model.getSoluteAllAtom_AtomCoordZ(0);

			soluteallatom_particle_rdf[ind1].atom_sphere_radius = solute_model.getSoluteAllAtom_AtomRadius(0);
		}

		std::cout << "[rdf_methods.cpp(CRdf::setUpSoluteAllAtomParticleRdf)] for particle type, "
				<< particle_types.getLabel(ind1) << std::endl;
		std::cout
		<< "[rdf_methods.cpp(CRdf::setUpSoluteAllAtomParticleRdf)]: rdf ref x coord (A) = "
		<< soluteallatom_particle_rdf[ind1].ref_xcoord << std::endl;
		std::cout
		<< "[rdf_methods.cpp(CRdf::setUpSoluteAllAtomParticleRdf)]: rdf ref y coord (A) = "
		<< soluteallatom_particle_rdf[ind1].ref_ycoord << std::endl;
		std::cout
		<< "[rdf_methods.cpp(CRdf::setUpSoluteAllAtomParticleRdf)]: rdf ref z coord (A) = "
		<< soluteallatom_particle_rdf[ind1].ref_zcoord << std::endl;
		std::cout << "-----" << std::endl;

		switch (soluteallatom_particle_rdf[ind1].cylindrical_rdf_axis_dir) {

		case X_AXIS:

			soluteallatom_particle_rdf[ind1].ref_axis_length =
					parameters.rdf_ref_axis_length_fraction
					* (solute_model.getSoluteAllAtom_MaxX() - solute_model.getSoluteAllAtom_MinX());

			soluteallatom_particle_rdf[ind1].coeff_x = 0.0;
			soluteallatom_particle_rdf[ind1].coeff_y = 1.0;
			soluteallatom_particle_rdf[ind1].coeff_z = 1.0;

			break;

		case Y_AXIS:

			soluteallatom_particle_rdf[ind1].ref_axis_length =
					parameters.rdf_ref_axis_length_fraction
					* (solute_model.getSoluteAllAtom_MaxY() - solute_model.getSoluteAllAtom_MinY());

			soluteallatom_particle_rdf[ind1].coeff_x = 1.0;
			soluteallatom_particle_rdf[ind1].coeff_y = 0.0;
			soluteallatom_particle_rdf[ind1].coeff_z = 1.0;
			break;

		case Z_AXIS:
			soluteallatom_particle_rdf[ind1].ref_axis_length =
					parameters.rdf_ref_axis_length_fraction
					* (solute_model.getSoluteAllAtom_MaxZ() - solute_model.getSoluteAllAtom_MinZ());

			soluteallatom_particle_rdf[ind1].coeff_x = 1.0;
			soluteallatom_particle_rdf[ind1].coeff_y = 1.0;
			soluteallatom_particle_rdf[ind1].coeff_z = 0.0;
			break;

		default:
			std::cerr
			<< "[rdf_methods.cpp(CRdf::setUpSoluteAllAtomParticleRdf)]: cylindrical rdf axis is not specified."
			<< std::endl;
			exit(0);
		}

		// set the parameters for the calculating the rdf of ions around a solute described in atomic detail,
		// its center fixed at the point (0,0,0).

		soluteallatom_particle_rdf[ind1].nrdfbin = nrdfbin;
		soluteallatom_particle_rdf[ind1].particletype_index = ind1;
		soluteallatom_particle_rdf[ind1].particletype_label = particle_types.getLabel(ind1);
		soluteallatom_particle_rdf[ind1].delr_rdf = delr_rdf;
		soluteallatom_particle_rdf[ind1].ndist_icount = 0;

		soluteallatom_particle_rdf[ind1].solute_x_min =  solute_model.getSoluteAllAtom_MinX();
		soluteallatom_particle_rdf[ind1].solute_y_min =  solute_model.getSoluteAllAtom_MinY();
		soluteallatom_particle_rdf[ind1].solute_z_min =  solute_model.getSoluteAllAtom_MinZ();

		for (int ind2 = 0; ind2 < nrdfbin; ind2++) {
			soluteallatom_particle_rdf[ind1].data->ndist[ind2] = 0;
			soluteallatom_particle_rdf[ind1].data->ndist_avg[ind2] = 0.0;
			soluteallatom_particle_rdf[ind1].data->ndist_avg_sqr[ind2] = 0.0;
			soluteallatom_particle_rdf[ind1].data->nbin[ind2] = 0;
			soluteallatom_particle_rdf[ind1].data->nbin_avg[ind2] = 0.0;
			soluteallatom_particle_rdf[ind1].data->nbin_avg_sqr[ind2] = 0.0;
			soluteallatom_particle_rdf[ind1].data->grdf[ind2] = 0.0;

			soluteallatom_particle_rdf[ind1].data->ndist_prev[ind2] = 0;
			soluteallatom_particle_rdf[ind1].data->ndist_prev[ind2] = 0;

			/* radial distance is calculated from the cylindrical axis
			 * rdist array elements refer to the bin mid-points.
			 */
			switch (parameters.rdf_ref) {

			case AXIS_OF_CYLINDER:
				soluteallatom_particle_rdf[ind1].rdist[ind2] = particle_types.getRadius(ind1)
				+ ind2 * delr_rdf;
				break;

			case CENTER_OF_SPHERE:
				soluteallatom_particle_rdf[ind1].rdist[ind2] =
						soluteallatom_particle_rdf[ind1].atom_sphere_radius - 4.0 * delr_rdf
						+ particle_types.getRadius(ind1) + ind2 * delr_rdf;
				break;
			}

		}

	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Calculate the radial distance of the selected particle from the solute's reference axis or point,
 * 			and update the bin count after a state change
 *
 * @param[in]	state_change(enum_state_change)	Type of state change
 * @param[in]	particle_types(CParticleType_t)	Particle type parameters and cavity grid data
 * @param[in]	particle_state(CParticleState_t)	Particle state data
 * @param[in]	ptcl(struct_particle)	Struct data of selected particle
 * @param[in]	ptclindex(int)		Array/vector index of the selected particle in the simulation box
 * @param[in]	ptcltype_index(int)	Array/vevtor index of the selected particle type
 *
 * @return void
 */
void CRdf::calcSoluteParticleBinCountsPerStateChange(enum_state_change state_change,const CParticleType_t &particle_types,
		const CParticleState_t &particle_state,const struct_particle &ptcl,int ptclindex, int ptcltype_index){

	double xij, yij, zij, rij;
	int ind1 = ptcltype_index;  // particle type index


	double nrdfbin = soluteallatom_particle_rdf[ind1].nrdfbin;  // number of bins
	double delr_rdf = soluteallatom_particle_rdf[ind1].delr_rdf; // bin size

	switch (soluteparticle_rdfref){ // select cylindrical axis or point as reference
	case AXIS_OF_CYLINDER:  // cylinder axis as reference
	{

		double axial_dist;
		double half_height = 0.5 * this->soluteallatom_particle_rdf[ind1].ref_axis_length;

		// copy previous particle count in each bin to current count
		for (int ind2 = 0; ind2 < nrdfbin; ind2++) {
			soluteallatom_particle_rdf[ind1].data->nbin[ind2] = soluteallatom_particle_rdf[ind1].data->nbin_prev[ind2];
			soluteallatom_particle_rdf[ind1].data->ndist[ind2] = soluteallatom_particle_rdf[ind1].data->ndist_prev[ind2];
		}


		xij = soluteallatom_particle_rdf[ind1].coeff_x
				* (ptcl.x - soluteallatom_particle_rdf[ind1].ref_xcoord);
		yij = soluteallatom_particle_rdf[ind1].coeff_y
				* (ptcl.y - soluteallatom_particle_rdf[ind1].ref_ycoord);
		zij = soluteallatom_particle_rdf[ind1].coeff_z
				* (ptcl.z - soluteallatom_particle_rdf[ind1].ref_zcoord);


		switch (soluteallatom_particle_rdf[ind1].cylindrical_rdf_axis_dir) {
		case X_AXIS:

			/* calculate axial distance of particle from the mid point of the solute length along x-axis direction */
			axial_dist = ptcl.x - soluteallatom_particle_rdf[ind1].ref_xcoord;
			rij = sqrt(yij * yij + zij * zij);
			break;
		case Y_AXIS:
			/* calculate axial distance of particle from the mid point of the solute length along y-axis direction */
			axial_dist = ptcl.y - soluteallatom_particle_rdf[ind1].ref_ycoord;
			rij = sqrt(xij * xij + zij * zij);
			break;
		case Z_AXIS:
			/* calculate axial distance of particle from the mid point of the solute length along z-axis direction */
			axial_dist = ptcl.z - soluteallatom_particle_rdf[ind1].ref_zcoord;
			rij = sqrt(xij * xij + yij * yij);
			break;

		default:

			std::cerr
			<< "[rdf_methods.cpp(CRdf::calcSoluteParticleBinCountsPerStateChange)] cylindrical rdf axis is not specified."
			<< std::endl;

			exit(0);
		}


		if (axial_dist < 0.0) {
			axial_dist = -1.0 * axial_dist;
		}

		double val = 1.0*RDF_SCALING_FACTOR;
		if (axial_dist <= half_height) {
			double rdist = soluteallatom_particle_rdf[ind1].rdist[0] - 0.5 * delr_rdf;
			int pos_bin = int(	(rij - rdist) / delr_rdf);
			if(pos_bin < nrdfbin){
				if(state_change==insertion){

					soluteallatom_particle_rdf[ind1].data->nbin[pos_bin] += val;

					for (int ind2 = pos_bin; ind2 < nrdfbin; ind2++) {
						soluteallatom_particle_rdf[ind1].data->ndist[ind2]+= val;
						soluteallatom_particle_rdf[ind1].data->ndist_prev[ind2] =soluteallatom_particle_rdf[ind1].data->ndist[ind2] ;
					}

				}else if(state_change==deletion){

					soluteallatom_particle_rdf[ind1].data->nbin[pos_bin] -= val;

					for (int ind2 = pos_bin; ind2 < nrdfbin; ind2++) {
						soluteallatom_particle_rdf[ind1].data->ndist[ind2] -= val;
						soluteallatom_particle_rdf[ind1].data->ndist_prev[ind2] =soluteallatom_particle_rdf[ind1].data->ndist[ind2] ;
					}

				}

				soluteallatom_particle_rdf[ind1].data->nbin_prev[pos_bin]  = soluteallatom_particle_rdf[ind1].data->nbin[pos_bin];

			}

			/*for (int pos_bin = 0; pos_bin < nrdfbin; pos_bin++) {
						double rdist = soluteallatom_particle_rdf[ind1].rdist[pos_bin] + 0.5 * delr_rdf;


						if (rij < rdist) {
							++soluteallatom_particle_rdf[ind1].data->ndist[pos_bin];
						}
					}*/

		}

	}
		break;

		case CENTER_OF_SPHERE:
		{

			// copy previous particle count in each bin to current count
			for (int ind2 = 0; ind2 < nrdfbin; ind2++) {
				soluteallatom_particle_rdf[ind1].data->nbin[ind2] = soluteallatom_particle_rdf[ind1].data->nbin_prev[ind2];
				soluteallatom_particle_rdf[ind1].data->ndist[ind2] = soluteallatom_particle_rdf[ind1].data->ndist_prev[ind2];

			}

			xij = ptcl.x - soluteallatom_particle_rdf[ind1].ref_xcoord;
			yij = ptcl.y - soluteallatom_particle_rdf[ind1].ref_ycoord;
			zij = ptcl.z - soluteallatom_particle_rdf[ind1].ref_zcoord;

			rij = sqrt(xij * xij + yij * yij + zij * zij);

			double val= 1.0*RDF_SCALING_FACTOR;
			double rdist = soluteallatom_particle_rdf[ind1].rdist[0] - 0.5 * delr_rdf;
			int pos_bin = int(	(rij - rdist) / delr_rdf);
			if(pos_bin < nrdfbin){
				if(state_change==insertion){

					soluteallatom_particle_rdf[ind1].data->nbin[pos_bin] += val;
					for (int ind2 = pos_bin; ind2 < nrdfbin; ind2++) {
						soluteallatom_particle_rdf[ind1].data->ndist[ind2] += val;
						soluteallatom_particle_rdf[ind1].data->ndist_prev[ind2] =soluteallatom_particle_rdf[ind1].data->ndist[ind2] ;
					}

				}else if(state_change==deletion){

					soluteallatom_particle_rdf[ind1].data->nbin[pos_bin]-= val;
					for (int ind2 = pos_bin; ind2 < nrdfbin; ind2++) {
						soluteallatom_particle_rdf[ind1].data->ndist[ind2] -= val;
						soluteallatom_particle_rdf[ind1].data->ndist_prev[ind2] =soluteallatom_particle_rdf[ind1].data->ndist[ind2] ;
					}
				}

				soluteallatom_particle_rdf[ind1].data->nbin_prev[pos_bin]  = soluteallatom_particle_rdf[ind1].data->nbin[pos_bin];

			}

			/* cumulative sum of ion numbers within a certain radial distance*/

			/*for (int pos_bin = 0; pos_bin < nrdfbin; pos_bin++) {
						double	rdist = soluteallatom_particle_rdf[ind1].rdist[pos_bin] + 0.5 * delr_rdf;
						//rdist = (pos_bin+1)*delr_rdf;

						if (rij < rdist) {
							++soluteallatom_particle_rdf[ind1].data->ndist[pos_bin];

						}
					}*/


		}
			break;
	}


	/* calculate number of ions in each bin for each segment */
	/*	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {

		soluteallatom_particle_rdf[ind1].data->nbin[0] =  soluteallatom_particle_rdf[ind1].data->ndist[0];

		for (int pos_bin = 1; pos_bin < nrdfbin; pos_bin++) {
			soluteallatom_particle_rdf[ind1].data->nbin[pos_bin] =
					soluteallatom_particle_rdf[ind1].data->ndist[pos_bin]
					                                             -  soluteallatom_particle_rdf[ind1].data->ndist[pos_bin - 1];
		}

	}*/


}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Compute the number of each particle type in the simulation box, within each
 * 		bin's radial distance range from the solute's reference axis or point
 *
 * @param[in]	particle_types(CParticleType_t)	Particle type parameters and cavity grid data
 * @param[in]	particle_state(CParticleState_t)	Particle state data
 *
 * @return void
 */
void CRdf::calcSoluteParticleBinCounts(const CParticleType_t &particle_types,
		const CParticleState_t &particle_state){

	switch (soluteparticle_rdfref){
	case AXIS_OF_CYLINDER:
	{

		double axial_dist;
		double rij;
		double half_height;
		for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
			half_height = 0.5 * this->soluteallatom_particle_rdf[ind1].ref_axis_length;

			double nrdfbin = soluteallatom_particle_rdf[ind1].nrdfbin;
			double delr_rdf = soluteallatom_particle_rdf[ind1].delr_rdf;

			for (int ind2 = 0; ind2 < nrdfbin; ind2++) {
				soluteallatom_particle_rdf[ind1].data->ndist[ind2] = 0.0;
				soluteallatom_particle_rdf[ind1].data->nbin[ind2] = 0.0;
			}

			for (int ind2=0;ind2 < particle_types.getNum(ind1); ind2++){

				int ptclindex = particle_state.getParticleIndicesByCount(ind1,ind2);

				struct_particle ptcl = particle_state.getParticle(ptclindex);

				double xij = soluteallatom_particle_rdf[ind1].coeff_x
						* (ptcl.x - soluteallatom_particle_rdf[ind1].ref_xcoord);
				double	yij = soluteallatom_particle_rdf[ind1].coeff_y
						* (ptcl.y - soluteallatom_particle_rdf[ind1].ref_ycoord);
				double		zij = soluteallatom_particle_rdf[ind1].coeff_z
						* (ptcl.z - soluteallatom_particle_rdf[ind1].ref_zcoord);


				switch (soluteallatom_particle_rdf[ind1].cylindrical_rdf_axis_dir) {
				case X_AXIS:

					/* calculate axial distance of particle from the mid point of the solute length along x-axis direction */
					axial_dist = ptcl.x - soluteallatom_particle_rdf[ind1].ref_xcoord;
					rij = sqrt(yij * yij + zij * zij);
					break;
				case Y_AXIS:
					/* calculate axial distance of particle from the mid point of the solute length along y-axis direction */
					axial_dist = ptcl.y - soluteallatom_particle_rdf[ind1].ref_ycoord;
					rij = sqrt(xij * xij + zij * zij);
					break;
				case Z_AXIS:
					/* calculate axial distance of particle from the mid point of the solute length along z-axis direction */
					axial_dist = ptcl.z - soluteallatom_particle_rdf[ind1].ref_zcoord;
					rij = sqrt(xij * xij + yij * yij);
					break;

				default:

					std::cerr
					<< "[rdf_methods.cpp(CRdf::calcSoluteParticleBinCounts)] cylindrical rdf axis is not specified."
					<< std::endl;

					exit(0);
				}


				if (axial_dist < 0.0) {
					axial_dist = -1.0 * axial_dist;
				}

				double val= 1.0*RDF_SCALING_FACTOR;
				if (axial_dist <= half_height) {
					double rdist =soluteallatom_particle_rdf[ind1].rdist[0] - 0.5*delr_rdf;
					int pos_bin = int((rij - rdist) / delr_rdf);
					if (pos_bin < nrdfbin){
						soluteallatom_particle_rdf[ind1].data->nbin[pos_bin] += val;



						/*		for (int pos_bin = 0; pos_bin < nrdfbin; pos_bin++) {
						double rdist = soluteallatom_particle_rdf[ind1].rdist[pos_bin] + 0.5 * delr_rdf;


						if (rij < rdist) {
							++soluteallatom_particle_rdf[ind1].data->ndist[pos_bin];
						}
					}*/

					}

				}

			}

		}
	}

			break;

	case CENTER_OF_SPHERE: //[This needs to be checked and corrected]
	{
		for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
			double nrdfbin = soluteallatom_particle_rdf[ind1].nrdfbin;
			double delr_rdf = soluteallatom_particle_rdf[ind1].delr_rdf;

			for (int ind2 = 0; ind2 < nrdfbin; ind2++) {
				soluteallatom_particle_rdf[ind1].data->ndist[ind2] = 0.0;
				soluteallatom_particle_rdf[ind1].data->nbin[ind2] = 0.0;
			}

			double val = 1.0*RDF_SCALING_FACTOR;
			for (int ind2=0;ind2 < particle_types.getNum(ind1); ind2++){

				int ptclindex = particle_state.getParticleIndicesByCount(ind1,ind2);

				struct_particle ptcl = particle_state.getParticle(ptclindex);
				double 	xij = ptcl.x - soluteallatom_particle_rdf[ind1].ref_xcoord;
				double	yij = ptcl.y - soluteallatom_particle_rdf[ind1].ref_ycoord;
				double zij = ptcl.z - soluteallatom_particle_rdf[ind1].ref_zcoord;

				double	rij = sqrt(xij * xij + yij * yij + zij * zij);

				/* cumulative sum of ion numbers within a certain radial distance*/


				for (int pos_bin = 0; pos_bin < nrdfbin; pos_bin++) {
					double	rdist = soluteallatom_particle_rdf[ind1].rdist[pos_bin] + 0.5 * delr_rdf;
					//rdist = (pos_bin+1)*delr_rdf;

					if (rij < rdist) {
						soluteallatom_particle_rdf[ind1].data->ndist[pos_bin] += val;

					}
				}
			}
		}

	}
		break;
		}

		/* calculate number of ions in each bin for each segment */
		for (int ind1 = 0; ind1 < nparticle_types; ind1++) {

			soluteallatom_particle_rdf[ind1].data->ndist[0] =  soluteallatom_particle_rdf[ind1].data->nbin[0];


			// save a copy
			soluteallatom_particle_rdf[ind1].data->nbin_prev[0] = soluteallatom_particle_rdf[ind1].data->nbin[0];
			soluteallatom_particle_rdf[ind1].data->ndist_prev[0] = soluteallatom_particle_rdf[ind1].data->ndist[0];

			for (int pos_bin = 1; pos_bin < nrdfbin; pos_bin++) {


				soluteallatom_particle_rdf[ind1].data->ndist[pos_bin] =
						soluteallatom_particle_rdf[ind1].data->ndist[pos_bin-1]
						                                             + soluteallatom_particle_rdf[ind1].data->nbin[pos_bin];

				// save a copy
				soluteallatom_particle_rdf[ind1].data->nbin_prev[pos_bin] = soluteallatom_particle_rdf[ind1].data->nbin[pos_bin];
				soluteallatom_particle_rdf[ind1].data->ndist_prev[pos_bin] = soluteallatom_particle_rdf[ind1].data->ndist[pos_bin];


			}


		}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Increment the total bin counts after sampling each state of the system for computing
 * 			the average number of each particle type in each bin
 *
 * @return void
 *
 */
void CRdf::incrementSoluteParticleBinAvgCounts(){
	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {

		int nbin = soluteallatom_particle_rdf[ind1].nrdfbin;
		for (int ind2 = 0; ind2 < nbin; ind2++) {

			soluteallatom_particle_rdf[ind1].data->ndist_avg[ind2] +=
					soluteallatom_particle_rdf[ind1].data->ndist[ind2];
			soluteallatom_particle_rdf[ind1].data->ndist_avg_sqr[ind2] +=
					soluteallatom_particle_rdf[ind1].data->ndist[ind2]
					                                             * soluteallatom_particle_rdf[ind1].data->ndist[ind2];
			soluteallatom_particle_rdf[ind1].data->nbin_avg[ind2] +=
					soluteallatom_particle_rdf[ind1].data->nbin[ind2];
			soluteallatom_particle_rdf[ind1].data->nbin_avg_sqr[ind2]+=
					soluteallatom_particle_rdf[ind1].data->nbin[ind2]
					                                            * soluteallatom_particle_rdf[ind1].data->nbin[ind2];

		}

	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Compute the average bin counts for each solute-particle type
 *
 * @param[in]	state_count(int)	Number of states sampled for compute the RDFs
 *
 * @return void
 *
 */
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void CRdf::finalizeSoluteParticleBinAvgCounts(int state_count){

	double inv_rdf_scaling_factor = 1.0/RDF_SCALING_FACTOR;
	double val = inv_rdf_scaling_factor/double(state_count);

	double val2 = val*inv_rdf_scaling_factor;

	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {


		int nbin = soluteallatom_particle_rdf[ind1].nrdfbin;

		soluteallatom_particle_rdf[ind1].ndist_icount = state_count;
		/* --- compute the radial distribution profile of ions around each segment --- */



		for (int ind2 = 0; ind2 < nbin; ind2++) {
			soluteallatom_particle_rdf[ind1].data->ndist_avg[ind2] *= val;
			soluteallatom_particle_rdf[ind1].data->ndist_avg_sqr[ind2] *= val2;
			soluteallatom_particle_rdf[ind1].data->nbin_avg[ind2] *= val;
			soluteallatom_particle_rdf[ind1].data->nbin_avg_sqr[ind2]*= val2;


		}

	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief	Calculate the solute-particle type RDF
 *
 * @param[in]	 @param[in]	particle_types(CParticleType_t)	Particle type parameters and cavity grid data
 *
 * @return void
 *
 */
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void CRdf::calcSoluteParticleRDF(const CParticleType_t &particle_types) {
	double bin_volume;


	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {
		double delr = soluteallatom_particle_rdf[ind1].delr_rdf;// bin size is same for all ion types
		int nbin = soluteallatom_particle_rdf[ind1].nrdfbin;

		double ptcl_num_density = concMolarityToNum(particle_types.getConc(ind1), box_vol)/ box_vol; // reference ion number density used to normalize the rdf
		double height = 0;
		switch(soluteparticle_rdfref){


		case AXIS_OF_CYLINDER:
			height = soluteallatom_particle_rdf[ind1].ref_axis_length;

			for (int ind2 = 0; ind2 < nbin; ind2++) {

				if (ind2 == 0) {
					bin_volume = PI * height
							* (pow(soluteallatom_particle_rdf[ind1].rdist[ind2] + 0.5 * delr, 2)
									- pow(soluteallatom_particle_rdf[ind1].rdist[ind2], 2));
				} else {
					bin_volume = PI* height
							* (pow(soluteallatom_particle_rdf[ind1].rdist[ind2] + 0.5 * delr, 2)
									- pow(soluteallatom_particle_rdf[ind1].rdist[ind2]
									                                             - 0.5 * delr, 2));
				}

				soluteallatom_particle_rdf[ind1].data->grdf[ind2] =
						soluteallatom_particle_rdf[ind1].data->nbin_avg[ind2]
						                                                / (bin_volume * ptcl_num_density);


			}
			break;

		case CENTER_OF_SPHERE:

			for (int ind2 = 0; ind2 < nbin; ind2++) {

				if (ind2 == 0) {
					bin_volume = 4.0 * PI
							* (pow(soluteallatom_particle_rdf[ind1].rdist[ind2] + 0.5 * delr, 3)
									- pow(soluteallatom_particle_rdf[ind1].rdist[ind2], 3)) / 3.0;
				} else {
					bin_volume = 4.0 * PI
							* (pow(soluteallatom_particle_rdf[ind1].rdist[ind2] + 0.5 * delr, 3)
									- pow(
											soluteallatom_particle_rdf[ind1].rdist[ind2]
											                                       - 0.5 * delr, 3)) / 3.0;
				}

				soluteallatom_particle_rdf[ind1].data->grdf[ind2] =
						soluteallatom_particle_rdf[ind1].data->nbin_avg[ind2]
						                                                / (bin_volume * ptcl_num_density);


			}

			break;


		}

	}

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Write solute particle type RDF to file
 *
 * The file name has the format rdf_solute_{particle type label}.out
 *
 * @return void
 */
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void CRdf::writeSoluteParticleRDFtoFile() {
	std::ofstream writefile;

	const char* fname = "";

	for (int ind1 = 0; ind1 < nparticle_types; ind1++) {

		std::string filename = "./outputfiles/rdf_solute_";


		filename.append(soluteallatom_particle_rdf[ind1].particletype_label);
		filename.append(".out");
		fname = filename.c_str();
		//writefile.open(fname, std::ios::app);
		writefile.open(fname);
		if(soluteparticle_rdfref==AXIS_OF_CYLINDER){
			writefile << "solute axis length (A) = "
					<< soluteallatom_particle_rdf[ind1].ref_axis_length << std::endl;
		}else if(soluteparticle_rdfref == CENTER_OF_SPHERE){
			writefile << "solute (atom) radius (A) = "
					<< soluteallatom_particle_rdf[ind1].atom_sphere_radius << std::endl;
		}
		writefile << "solute charge (e units) = " << solute_total_charge
				<< std::endl;


		writefile << "r_of_bin_mid-point_in_A" << '\t' << "rdf"
				<< '\t' << "avg_cumulative_sum" << '\t'
				<< "sd_cumulative_sum" << '\t' << "avg_particle_count_in_bin"
				<< '\t' << "sd_particle_count_in_bin" << std::endl;

		int nbin = soluteallatom_particle_rdf[ind1].nrdfbin;
		for (int ind2 = 0; ind2 < nbin; ind2++) {
			writefile << soluteallatom_particle_rdf[ind1].rdist[ind2] << '\t'
					<< soluteallatom_particle_rdf[ind1].data->grdf[ind2] << '\t'
					<< soluteallatom_particle_rdf[ind1].data->ndist_avg[ind2] << '\t'
					<< sqrt(soluteallatom_particle_rdf[ind1].data->ndist_avg_sqr[ind2]
					                                                             - soluteallatom_particle_rdf[ind1].data->ndist_avg[ind2]
					                                                                                                                * soluteallatom_particle_rdf[ind1].data->ndist_avg[ind2])
			 << '\t' << soluteallatom_particle_rdf[ind1].data->nbin_avg[ind2] << '\t'<<
			 sqrt(soluteallatom_particle_rdf[ind1].data->nbin_avg_sqr[ind2]
			 	   - soluteallatom_particle_rdf[ind1].data->nbin_avg[ind2] * soluteallatom_particle_rdf[ind1].data->nbin_avg[ind2]) << std::endl;
		}
		writefile << std::endl;
		writefile.close();
	}



}
