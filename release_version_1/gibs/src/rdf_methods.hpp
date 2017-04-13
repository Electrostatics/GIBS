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
 * rdf_methods.hpp
 *
 *  Created on: Jul 9, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file rdf_methods.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Class definition for radial distribution functions
 *
 */

#ifndef RDF_METHODS_HPP_
#define RDF_METHODS_HPP_

#include "generic.hpp"
#include "input_parameters.hpp"
#include "particle_type.hpp"
#include "particlepair_type.hpp"
#include "rdf_structs.hpp"
#include "constants.hpp"
#include "particle_state.hpp"
#include "convert_functions.hpp"

#include "solute_type.hpp"
/**
 * @brief Class definition storing RDF data and parameters
 *
 * Use defined type CRdf_t for class declarations
 */
class CRdf {

	enum_solute_models solute_model_type;
	enum_rdf_ref soluteparticle_rdfref;

	int nparticle_types;
	int nparticle_pairs;
	int nrdfbin;
	double delr_rdf;
	double xlen;
	double ylen;
	double zlen;
	double xmin;
	double ymin;
	double zmin;
	double box_vol;
	double solute_total_charge;
	bool calculate_particlepair_rdf;
	bool calculate_solute_particle_rdf;

	struct_particlepair_rdf *particlepair_rdf;
	struct_soluteallatom_particle_rdf *soluteallatom_particle_rdf;
public:
	// constructor
	CRdf(const CInputParameters_t &parameters, const CParticleType_t &particle_types,
			const CParticlePairType_t &particlepair_types,const CSoluteModel_t &solute_model,double box_vol);

	// set methods
	void setUpSoluteAllAtomParticleRdf(const CInputParameters_t &parameters,
			const CParticleType_t &particle_types,const CSoluteModel_t &solute_model);

	// destructor
	~CRdf(){

		for(int i=0;i<nparticle_pairs;i++){
			delete[] particlepair_rdf[i].data->grdf;
			delete[] particlepair_rdf[i].data->nbin;
			delete[] particlepair_rdf[i].data->nbin_avg;
			delete[] particlepair_rdf[i].data->nbin_avg_sqr;
			delete[] particlepair_rdf[i].data->ndist;
			delete[] particlepair_rdf[i].data->ndist_avg;
			delete[] particlepair_rdf[i].data->ndist_avg_sqr;
			delete[] particlepair_rdf[i].rdist;

			delete[] particlepair_rdf[i].data->nbin_prev;
			delete[] particlepair_rdf[i].data->ndist_prev;

			std::cout << "[rdf.hpp (class CRdf)]: deallocated memory for rdf data of particle pair type, " <<
					particlepair_rdf[i].particlepair_label << std::endl;

		}

		delete[] particlepair_rdf;

		if(solute_model_type==all_atom && this->calculate_solute_particle_rdf==true){
			for (int i = 0; i < nparticle_types; i++){
				delete[] soluteallatom_particle_rdf[i].data->nbin;
				delete[] soluteallatom_particle_rdf[i].data->nbin_avg;
				delete[] soluteallatom_particle_rdf[i].data->nbin_avg_sqr;
				delete[] soluteallatom_particle_rdf[i].data->ndist;
				delete[] soluteallatom_particle_rdf[i].data->ndist_avg;
				delete[] soluteallatom_particle_rdf[i].data->ndist_avg_sqr;
				delete[] soluteallatom_particle_rdf[i].rdist;

				delete[] soluteallatom_particle_rdf[i].data->ndist_prev;
				delete[] soluteallatom_particle_rdf[i].data->nbin_prev;
			}
			delete[] soluteallatom_particle_rdf;
			soluteallatom_particle_rdf = NULL;
		}
	}

	// calculation methods for particle pair RDFs
	void calcParticlePairRDF(const CParticleType_t &particle_types);
	void calcParticlePairBinCounts(const CParticleType_t &particle_types,const CParticleState_t &particle_state);
	void incrementParticlePairBinAvgCounts();
	void finalizeParticlePairBinAvgCounts(int state_count);

	void calcParticlePairBinCountsPerStateChange(enum_state_change state_change,const CParticleType_t &particle_types,
			const CParticlePairType_t &particlepair_types,
			const CParticleState_t &particle_state, const struct_particle &ptcl,int ptclindex, int ptcltype_index);

	// for solute-particle RDFs
	void calcSoluteParticleBinCountsPerStateChange(enum_state_change state_change,const CParticleType_t &particle_types,
			const CParticleState_t &particle_state,const struct_particle &ptcl,int ptclindex, int ptcltype_index);
	void calcSoluteParticleBinCounts(const CParticleType_t &particle_types,
			const CParticleState_t &particle_state);
	void incrementSoluteParticleBinAvgCounts();
	void finalizeSoluteParticleBinAvgCounts(int state_count);
	void calcSoluteParticleRDF(const CParticleType_t &particle_types);

	// get methods
	int getNParticleTypes()const {return nparticle_types;}
	int getNParticlePairs() const {return nparticle_pairs;}
	int getNumBins() const {return nrdfbin;}
	double getBinWidth() const{return delr_rdf;}

	// write methods
	void writeParticlePairRDFtoFile();
	void writeSoluteParticleRDFtoFile();

	// reset methods
	void resetData();

};

typedef class CRdf CRdf_t;

#endif /* RDF_METHODS_HPP_ */
