/*
 * energy.cpp
 *
 *  Created on: Jul 15, 2015
 *      Author: Dennis G. Thomas
 *
 *  @file energy.cpp
 *  @author Dennis G. Thomas
 *
 *  @brief Functions for calculating energies
 */

#include "energy.hpp"

/**
 * @brief Calculates the total interaction energy between particles of the same type
 *
 * @param[in]	state_change(bool)					True or False option to indicate whether the state of the system has changed or not
 * @param[in]	use_spm(bool)						True or False option to indicate whether to use SPM or not
 * @param[in]	use_mwater(bool)					True or False option to indicate whether to use MWater model
 * @param[in]	particle_types(CParticleType_t)		Class instance containing parameters of each particle type
 * @param[in]	particle_state(CParticleState_t)	Class instance containing particle state data
 * @param[in]	particlepair_types(CParticlePairType_t)	Class instance containing parameters of each particle pair type
 * @param[in]	mwater(CMWater)						Class instance containing parameters of MWater model
 * @param[in]	cavity_grid(CCavityGrid_t)		Class instance containing parameters of the cavity grid data
 * @param[in]	box	(CBox_t)					Class instance containing parameters of the simulation box
 * @param[in]	solvdiel(double)		Solvent dielectric constant
 * @param[in]	ptclindex(int)		Vector element index  of particle state vector
 * @param[in]	ptcltype(int)		Vector element index of particle type vector
 *
 * @return 	Energy(double) in kcal/mol
 */
double calcParticlePairEnergy_SameType(bool state_change,bool use_spm,bool use_mwater,const CParticleType_t &particle_types,
		const CParticleState_t &particle_state, const CParticlePairType_t &particlepair_types,CMWater &mwater,
		const CCavityGrid_t &cavity_grid,const CBox_t &box,
		double solvdiel, int ptclindex, int ptcltype){


	double energy = 0.0;

	double xlen = box.getBoxXLen();
	double ylen = box.getBoxYLen();
	double zlen = box.getBoxZLen();

	int npart_types = particle_types.getNumParticleTypes();

	if(use_spm == true){
		// Reduce number of particle types by 1 to exclued solvent-solvent electrostatics
		--npart_types;

	}

//	switch (state_change){
	//case false:
	/* Calculate interaction energy between all particles of the same type */
		if(state_change==false){
		for (int ind1=0; ind1<npart_types; ind1++){
			bool clmb_potential = particlepair_types.getCoulombPotential(ind1);
			bool lj_potential = particlepair_types.getLennardJonesPotential(ind1);
			bool pmf_lookup_table = particlepair_types.getPMFLookTable(ind1);

			int ind1_numtype = particle_types.getNum(ind1);
			for (int ind2 =0; ind2 < ind1_numtype-1; ind2++){
				int i = particle_state.getParticleIndicesByCount(ind1,ind2);

				double xi = particle_state.getParticlePositionX(i);
				double yi = particle_state.getParticlePositionY(i);
				double zi = particle_state.getParticlePositionZ(i);
				double zvi =particle_state.getParticleCharge(i);



				for (int ind3 =ind2+1; ind3 < ind1_numtype; ind3++){

					int j = particle_state.getParticleIndicesByCount(ind1,ind3);

					double xj = particle_state.getParticlePositionX(j);
					double	yj = particle_state.getParticlePositionY(j);
					double	zj = particle_state.getParticlePositionZ(j);

					// Find the nearest image
					if (xj < xi - xlen * 0.5) {
						xj += xlen;
					} else if (xj > xi + xlen * 0.5) {
						xj -= xlen;
					}
					if (yj < yi - ylen * 0.5) {
						yj += ylen;
					} else if (yj > yi + ylen * 0.5) {
						yj -= ylen;
					}
					if (zj < zi - zlen * 0.5) {
						zj += zlen;
					} else if (zj > zi + zlen * 0.5) {
						zj -= zlen;
					}

					double xij = xi - xj;
					double yij = yi - yj;
					double zij = zi - zj;

					double rij = sqrt(xij * xij + yij * yij + zij * zij);

					double	zvj = particle_state.getParticleCharge(j);



					if(clmb_potential==1 && lj_potential==1){
						// Assumption: pairs of similar particles are listed in the same order as the particle types


						double ljpot = particlepair_types.getAij(ind1)/pow(rij,12) -
								particlepair_types.getBij(ind1)/pow(rij,6); //units in kcal/mol

						double clmbpot = 332.06364 * zvi * zvj / (solvdiel * rij); // units in kcal/mol

						energy += ljpot;
						energy += clmbpot;

					}else if(clmb_potential==1){
						double clmbpot = 332.06364 * zvi * zvj / (solvdiel * rij); // units in kcal/mol

						energy += clmbpot;
						//	std::cout << "[energy.cpp(calcPPEnergy_SameType)]: coulomb potential, energy]"
						//		<< energy << ", between pair types " << particlepair_types[ind1].label12 << std::endl;
						/*std::cout << "[energy.cpp(calcPPEnergy_SameType)]: all-to-all, coulomb potential, energy] rij = " <<
												rij << ", coulomb energy (kcal/mol) = " << clmbpot << ", pair label " << particlepair_types[ind1].label12 <<
												" between particle indices "  << i << " and " << j << std::endl;
						 */
					}
					else if(lj_potential ==1){

						double ljpot = particlepair_types.getAij(ind1)/pow(rij,12) -
								particlepair_types.getBij(ind1)/pow(rij,6);
						energy += ljpot;  // units in kcal/mol
					}
					if(pmf_lookup_table==1){
						int n_lookup=	particlepair_types.getNLookup(ind1);
						if (rij < particlepair_types.getRLookup(ind1,0)) {
							// give a very large value, if the radial distance between the two ions is
							// less than the lowest value found in the lookup table.
							energy = 10000.0;

							return (energy);
						} else if (rij > particlepair_types.getRLookup(ind1,n_lookup-1)) {

							// interaction energy for particle pairs at distances greater than the maximum distance recorded in the lookup table
							// are set to zero.
							energy += 0;
						} else {


							int indm = int(
									floor(
											(rij - particlepair_types.getRLookup(ind1,0))
											/ (particlepair_types.getRLookup(ind1,1)
													- particlepair_types.getRLookup(ind1,0))));
							if (indm != n_lookup - 1) {
								double x = (rij - particlepair_types.getRLookup(ind1,indm))
																																																																																																								/ (particlepair_types.getRLookup(ind1,1)
																																																																																																										- particlepair_types.getRLookup(ind1,0));
								energy += cspline_interp(
										particlepair_types.getSplineCoeff(ind1,indm).cubic, x);


							} else if (indm == n_lookup - 1) {

								energy += particlepair_types.getPMFLookup(ind1,n_lookup - 1);
							}
						}
					}
				}
			}
		}

		// If monoatomic water model is used with spm.
		if(use_mwater == true && use_spm == true){
			std::string label12 = "Water_Water";
			int pair_index = particlepair_types.getIndexFromLabel(label12);
			//double	hs_cutoff =  particlepair_types.getRadius1(pair_index) + particlepair_types.getRadius2(pair_index);
			double hs_cutoff = particlepair_types.getHardSphereCutoff(label12);
			int numtype = particle_types.getNum(npart_types);
			// find grid cells of solvent spheres that are within the M water potential cut off
			double mw_cutoff = mwater.asigma;

			// calculate half box length
			double xlenh = xlen*0.5;
			double ylenh = ylen*0.5;
			double zlenh = zlen*0.5;

			double minlen = std::min(xlen,ylen);
			minlen=std::min(minlen,zlen);
			double boxh14=0.25*minlen;

			for (int ind1=0;ind1< numtype-1; ind1++){
				int i = particle_state.getParticleIndicesByCount(npart_types,ind1);

				double xi = particle_state.getParticlePositionX(i);
				double yi = particle_state.getParticlePositionY(i);
				double zi = particle_state.getParticlePositionZ(i);


				for (int ind2=ind1+1; ind2 < numtype; ind2++){

					int j = particle_state.getParticleIndicesByCount(npart_types,ind2);

					double xj = particle_state.getParticlePositionX(j);
					double	yj = particle_state.getParticlePositionY(j);
					double	zj = particle_state.getParticlePositionZ(j);

					double	xij = xi - xj;
					double	yij = yi - yj;
					double	zij = zi - zj;

					// Find the nearest image
					if (xij > xlenh) {
						xij -= xlen;

					} else if (xij < -xlenh){
						xij += xlen;
					}
					if (yij > ylenh){
						yij -= ylen;

					} else if (yij < -ylenh){
						yij += ylen;

					}
					if (zij > zlenh){
						zij -= zlen;

					} else if (zij < -zlenh){
						zij += zlen;
					}

					double rij = sqrt(xij * xij + yij * yij + zij * zij);

					if (rij < mw_cutoff){
						// 2-body potential
						//					energy += (mwater.AepsilonB*pow(mwater.sigma/rij,mwater.p) -
						//						mwater.Aepsilon*pow(mwater.sigma/rij,mwater.q))*exp(mwater.sigma/(rij-mwater.asigma));

						// q = 0 for M water. So, simplifying the term containing q.
						double argij = mwater.sigma/(rij - mw_cutoff);
						if(argij > -30){

							energy += (mwater.AepsilonB*pow(mwater.sigma/rij,mwater.p) -mwater.Aepsilon)*exp(argij);

						}
					}
					// 3-body potential
					double expij = 0.0;
					if (rij < mw_cutoff){
						double argij = mwater.gammasigma/(rij - mw_cutoff);
						if (argij > -30){
							expij = exp(argij);
						}else{
							expij = 0.0;
						}
					}

					if(ind2 < numtype-1){
						for (int ind3=ind2+1;ind3 < numtype;ind3++){


							double cosjik = 0.0;
							double cosijk = 0.0;
							double cosjki = 0.0;



							int k = particle_state.getParticleIndicesByCount(npart_types,ind3);

							double xk = particle_state.getParticlePositionX(k);
							double yk = particle_state.getParticlePositionY(k);
							double zk = particle_state.getParticlePositionZ(k);


							// Find the nearest image between i and k
							double xik = xi - xk;
							double yik = yi - yk;
							double zik = zi - zk;

							// Find the nearest image
							if (xik > xlenh) {
								xik -= xlen;

							} else if (xik < -xlenh){
								xik += xlen;

							}
							if (yik > ylenh){
								yik -= ylen;

							} else if (yik < -ylenh){
								yik += ylen;

							}
							if (zik > zlenh){
								zik -= zlen;
								//zj += zlen;
							} else if (zik < -zlenh){
								zik += zlen;

							}

							double rik = sqrt(xik * xik + yik * yik + zik * zik);

							// Find the nearest image between j and k

							double xjk = xj - xk;
							double yjk = yj - yk;
							double zjk = zj - zk;

							// Find the nearest image
							if (xjk > xlenh) {
								xjk -= xlen;

							} else if (xjk < -xlenh){
								xjk += xlen;

							}
							if (yjk > ylenh){
								yjk -= ylen;

							} else if (yjk < -ylenh){
								yjk += ylen;

							}
							if (zjk > zlenh){
								zjk -= zlen;
								//zj += zlen;
							} else if (zjk < -zlenh){
								zjk += zlen;

							}

							double rjk = sqrt(xjk * xjk + yjk * yjk + zjk * zjk);

							double expik = 0.0;

							if (rik < mw_cutoff){
								double argik = mwater.gammasigma/(rik - mw_cutoff);
								if (argik > -30){
									expik = exp(argik);
								}else{
									expik = 0.0;
								}

							}

							double expjk = 0.0;
							if (rjk < mw_cutoff){
								double argjk = mwater.gammasigma/(rjk - mw_cutoff);
								if (argjk > -30){
									expjk = exp(argjk);
								}else{
									expjk = 0.0;
								}

							}

							double expa = expij*expik;
							if (expa !=0.0){
								double cosjik = (xij*xik + yij*yik + zij*zik)/(rik*rij) - mwater.costheta0;
								energy += mwater.lamdaepsilon*cosjik*cosjik*expa;
							}

							expa = expij*expjk;
							if (expa !=0.0){
								double cosijk = -(xij*xjk + yij*yjk + zij*zjk)/(rjk*rij) - mwater.costheta0;
								energy += mwater.lamdaepsilon*cosijk*cosijk*expa;
							}

							expa = expik*expjk;
							if (expa !=0.0){
								double cosikj = (xik*xjk + yik*yjk + zik*zjk)/(rjk*rik) - mwater.costheta0;
								energy += mwater.lamdaepsilon*cosikj*cosikj*expa;
							}


						}
					} // ind3 for-loop


				} //ind2 for-loop
			}  //  ind1 for-loop

		} // if(use_mwater == true && use_spm == true)

		//break;

	//case true:
		//std::cout << "[energy.cpp(calcPPEnergy_SameType)]: Calculating interaction energy one particle and other particles of the same type." << std::endl;
		//int ind1, i;
		}else if(state_change==true) {
			/* Calculate interaction energy between 1 particle and other particles of its type*/
		int ind1 = ptcltype-1;
		int i = ptclindex;
		std::string ptcl_label = particle_types.getLabel(ind1);



		if(ptcl_label == "Water" && use_mwater == true && use_spm == true){

			int nsolvent_indices = 0;

			std::string label12 = "Water_Water";
			int pair_index = particlepair_types.getIndexFromLabel(label12);
			//double	hs_cutoff =  particlepair_types.getRadius1(pair_index) + particlepair_types.getRadius2(pair_index);
			double hs_cutoff = particlepair_types.getHardSphereCutoff(label12);
			double mw_cutoff = mwater.asigma;
			int nparticle_types = particle_types.getNumParticleTypes();
			// calculate half box length
			double xlenh = xlen*0.5;
			double ylenh = ylen*0.5;
			double zlenh = zlen*0.5;

			double minlen = std::min(xlen,ylen);
			double boxh14 = 0.25*std::min(minlen,zlen);

			double xi = particle_state.getParticlePositionX(i);
			double yi = particle_state.getParticlePositionY(i);
			double zi = particle_state.getParticlePositionZ(i);


			int lix = particle_state.getParticleGridCellIndexX(i);
			int liy = particle_state.getParticleGridCellIndexY(i);
			int liz = particle_state.getParticleGridCellIndexZ(i);

			struct_cellindex_offset cellindex_offset = particlepair_types.getMWaterCellIndexOffsetMatrix(label12);
			int ncx = cellindex_offset.ncx;
			int ncy = cellindex_offset.ncy;
			int ncz = cellindex_offset.ncz;
			//	int Ncx = cellindex_offset.Ncx;
			//	int Ncy = cellindex_offset.Ncy;
			//	int Ncz = cellindex_offset.Ncz;

			int ljx,ljy,ljz;

			for (int ind3=0;ind3<2*ncx+1;ind3++){
				ljx = lix + cellindex_offset.delxj[ind3] + cellindex_offset.bcmap_x[lix][ind3];

				for (int ind4 = 0; ind4<2*ncy+1; ind4++){
					ljy = liy + cellindex_offset.delyj[ind4] + cellindex_offset.bcmap_y[liy][ind4];

					for (int ind5 = 0; ind5 < 2*ncz+1; ind5++){
						ljz = liz + cellindex_offset.delzj[ind5] + cellindex_offset.bcmap_z[liz][ind5];

						int cellindex_1D =  cavity_grid.getCellIndex1D(ljx,ljy,ljz);
						int solvent_index = particle_state.getParticleIndicesByCellindex(nparticle_types-1,cellindex_1D);
						if (solvent_index!=-1 && solvent_index!=ptclindex && nsolvent_indices == 0){
							mwater.solvent_indices[0] = solvent_index; // save solvent_indices
							++nsolvent_indices;
						}
						if (solvent_index!=-1 && solvent_index!=ptclindex && nsolvent_indices > 0){
							bool index_present=false;
							int isolv = 0;

							while (isolv <nsolvent_indices && index_present==false){
								if(solvent_index==mwater.solvent_indices[isolv]){
									index_present=true;
								}
								++isolv;
							}
							if(index_present==false){
								mwater.solvent_indices[nsolvent_indices] = solvent_index; // save solvent_indices
								++nsolvent_indices;
							}
						}

					}
				}
			}

			int nsolvent_indices_1 = nsolvent_indices;
			if(nsolvent_indices_1 > 0){
				for(int solvent_index=0; solvent_index < nsolvent_indices_1; solvent_index++){


					int ljx = particle_state.getParticleGridCellIndexX(mwater.solvent_indices[solvent_index]);
					int ljy = particle_state.getParticleGridCellIndexY(mwater.solvent_indices[solvent_index]);
					int ljz = particle_state.getParticleGridCellIndexZ(mwater.solvent_indices[solvent_index]);


					int lkx,lky,lkz;

					for (int ind3j=0;ind3j<2*ncx+1;ind3j++){
						lkx = ljx + cellindex_offset.delxj[ind3j] + cellindex_offset.bcmap_x[ljx][ind3j];

						for (int ind4j = 0; ind4j<2*ncy+1; ind4j++){
							lky = ljy + cellindex_offset.delyj[ind4j] + cellindex_offset.bcmap_y[ljy][ind4j];

							for (int ind5j = 0; ind5j < 2*ncz+1; ind5j++){
								lkz = ljz + cellindex_offset.delzj[ind5j] + cellindex_offset.bcmap_z[ljz][ind5j];

								int cellindex_1Dk =  cavity_grid.getCellIndex1D(lkx,lky,lkz);
								int solvent_index_k = particle_state.getParticleIndicesByCellindex(nparticle_types-1,cellindex_1Dk);

								if (solvent_index_k!=-1 && solvent_index_k!=mwater.solvent_indices[solvent_index] && solvent_index_k!=ptclindex){
									bool index_present=false;
									int isolv = 0;
									while (isolv <nsolvent_indices && index_present == false){
										if(solvent_index_k==mwater.solvent_indices[isolv]){
											index_present=true;
										}
										++isolv;
									}
									if(index_present==false){
										mwater.solvent_indices[nsolvent_indices] = solvent_index_k; // save solvent_indices
										++nsolvent_indices;
									}

								}
							}
						}
					}

				}

			}
			if (nsolvent_indices >=1){
				for (int j=0; j < nsolvent_indices;j++){
					int solvent_index = mwater.solvent_indices[j];

					double xj = particle_state.getParticlePositionX(solvent_index);
					double yj = particle_state.getParticlePositionY(solvent_index);
					double zj = particle_state.getParticlePositionZ(solvent_index);

					double	xij = xi - xj;
					double	yij = yi - yj;
					double	zij = zi - zj;


					// Find the nearest image
					if (xij > xlenh) {
						xij -= xlen;

					} else if (xij < -xlenh){
						xij += xlen;

					}
					if (yij > ylenh){
						yij -= ylen;

					} else if (yij < -ylenh){
						yij += ylen;

					}
					if (zij > zlenh){
						zij -= zlen;

					} else if (zij < -zlenh){
						zij += zlen;

					}

					double rij = sqrt(xij * xij + yij * yij + zij * zij);
					/*std::cout << "[energy.cpp]: neighbor solvent_index = " << solvent_index
							<< ", xj = " << xj << ", yj = " << yj << ", zj = " << zj <<
							", rij = " << rij << std::endl;*/

					//if (rij < mw_cutoff && rij >= hs_cutoff)
					if (rij < mw_cutoff){
						// 2-body potential

						// q = 0 for M water. So, simplifying the term containing q.
						double argij = mwater.sigma/(rij - mw_cutoff);
						if(argij > -30){

							//	double ener =(mwater.AepsilonB*pow(mwater.sigma/rij,mwater.p) -mwater.Aepsilon)*exp(argij);
							//		std::cout << "ener = " << ener << std::endl;
							energy += (mwater.AepsilonB*pow(mwater.sigma/rij,mwater.p) -mwater.Aepsilon)*exp(argij);
							//std::cout << "[energy.cpp()]: Energy (kcal/mol) = " << energy << std::endl;
						}
					}


					// 3-body potential

					if(nsolvent_indices >=2 && j < nsolvent_indices-1){

						double expij = 0.0;
						if (rij < mw_cutoff){
							double argij = mwater.gammasigma/(rij - mw_cutoff);
							if (argij > -30){
								expij = exp(argij);
							}else{
								expij = 0.0;
							}
						}


						for (int k=j+1;k<nsolvent_indices;k++){
							int solvent_index_k = mwater.solvent_indices[k];


							double xk = particle_state.getParticlePositionX(solvent_index_k);
							double yk = particle_state.getParticlePositionY(solvent_index_k);
							double zk = particle_state.getParticlePositionZ(solvent_index_k);



							// Find the nearest image between i and k
							double xik = xi - xk;
							double yik = yi - yk;
							double zik = zi - zk;

							// Find the nearest image
							if (xik > xlenh) {
								xik -= xlen;

							} else if (xik < -xlenh){
								xik += xlen;

							}
							if (yik > ylenh){
								yik -= ylen;

							} else if (yik < -ylenh){
								yik += ylen;

							}
							if (zik > zlenh){
								zik -= zlen;
								//zj += zlen;
							} else if (zik < -zlenh){
								zik += zlen;

							}

							double rik = sqrt(xik * xik + yik * yik + zik * zik);

							// Find the nearest image between j and k

							double xjk = xj - xk;
							double yjk = yj - yk;
							double zjk = zj - zk;


							// Find the nearest image
							if (xjk > xlenh) {
								xjk -= xlen;

							} else if (xjk < -xlenh){
								xjk += xlen;

							}
							if (yjk > ylenh){
								yjk -= ylen;

							} else if (yjk < -ylenh){
								yjk += ylen;

							}
							if (zjk > zlenh){
								zjk -= zlen;
								//zj += zlen;
							} else if (zjk < -zlenh){
								zjk += zlen;

							}

							double rjk = sqrt(xjk * xjk + yjk * yjk + zjk * zjk);

							double expik = 0.0;

							if(rij < boxh14 && rjk < boxh14 && rik < boxh14){
								if (rik < mw_cutoff){
									double argik = mwater.gammasigma/(rik - mw_cutoff);
									if (argik > -30){
										expik = exp(argik);
									}else{
										expik = 0.0;
									}

								}

								double expjk = 0.0;
								if (rjk < mw_cutoff){
									double argjk = mwater.gammasigma/(rjk - mw_cutoff);
									if (argjk > -30){
										expjk = exp(argjk);
									}else{
										expjk = 0.0;
									}

								}
								double expa = expij*expik;


								if (expa !=0.0){
									double cosjik = (xij*xik + yij*yik + zij*zik)/(rik*rij) - mwater.costheta0;
									energy += mwater.lamdaepsilon*cosjik*cosjik*expa;
								}

								expa = expij*expjk;
								if (expa !=0.0){
									double cosijk = -(xij*xjk + yij*yjk + zij*zjk)/(rjk*rij) - mwater.costheta0;
									energy += mwater.lamdaepsilon*cosijk*cosijk*expa;
								}

								expa = expik*expjk;
								if (expa !=0.0){
									double cosikj = (xik*xjk + yik*yjk + zik*zjk)/(rjk*rik) - mwater.costheta0;
									energy += mwater.lamdaepsilon*cosikj*cosikj*expa;
								}

							}
						}

					}
				}
			}

			/*std::cout << "[energy.cpp()]: Energy (kcal/mol) = " << energy << ", nsolvent_indices = "<<
					nsolvent_indices << ", ptclindex = " << ptclindex
					<< ", xi = " << xi << ", yi = " << yi << ", zi = " << zi << std::endl;*/

		} else if (ptcl_label != "Water"){ // only ions

			bool clmb_potential = particlepair_types.getCoulombPotential(ind1);
			bool lj_potential = particlepair_types.getLennardJonesPotential(ind1);
			bool pmf_lookup_table = particlepair_types.getPMFLookTable(ind1);

			double xi = particle_state.getParticlePositionX(i);
			double yi = particle_state.getParticlePositionY(i);
			double zi = particle_state.getParticlePositionZ(i);

			// for ions

			int ind1_numtype= particle_types.getNum(ind1);
			for (int ind3 =0; ind3 < ind1_numtype; ind3++){
				int j = particle_state.getParticleIndicesByCount(ind1,ind3);
				if(i!=j) {

					double xj = particle_state.getParticlePositionX(j);
					double	yj = particle_state.getParticlePositionY(j);
					double	zj = particle_state.getParticlePositionZ(j);

					// Find the nearest image
					if (xj < xi - xlen * 0.5) {
						xj += xlen;
					} else if (xj > xi + xlen * 0.5) {
						xj -= xlen;
					}
					if (yj < yi - ylen * 0.5) {
						yj += ylen;
					} else if (yj > yi + ylen * 0.5) {
						yj -= ylen;
					}
					if (zj < zi - zlen * 0.5) {
						zj += zlen;
					} else if (zj > zi + zlen * 0.5) {
						zj -= zlen;
					}


					double xij = xi - xj;
					double yij = yi - yj;
					double zij = zi - zj;

					double rij = sqrt(xij * xij + yij * yij + zij * zij);

					double zvi =particle_state.getParticleCharge(i);
					double	zvj = particle_state.getParticleCharge(j);

					if(clmb_potential==1 && lj_potential ==1){
						// Assumption: pairs of similar particles are listed in the same order as the particle types

						double ljpot = particlepair_types.getAij(ind1)/pow(rij,12) -
								particlepair_types.getBij(ind1)/pow(rij,6);

						double clmbpot = 332.06364 * zvi * zvj / (solvdiel * rij); // units in kcal/mol
						//energy = energy + ljpot + clmbpot;	// units in kcal/mol
						energy += ljpot;
						energy += clmbpot;

					}else if(clmb_potential==1){
						double clmbpot = 332.06364 * zvi * zvj / (solvdiel * rij); // units in kcal/mol
						//energy = energy + clmbpot;	// units in kcal/mol
						energy += clmbpot;
						//	std::cout << "[energy.cpp(calcPPEnergy_SameType)]: one-to-many, coulomb potential, energy]"
						//								<< energy << ", between pair types " << particlepair_types[ind1].label12 << std::endl;

						/*std::cout << "[energy.cpp(calcPPEnergy_SameType)]: one-to-many, coulomb potential, energy] rij = " <<
						rij << ", coulomb energy (kcal/mol) = " << clmbpot << ", pair label " << particlepair_types[ind1].label12 <<
						" between particle indices "  << i << " and " << j << std::endl;*/

					}
					else if(lj_potential ==1){

						double ljpot = particlepair_types.getAij(ind1)/pow(rij,12) -
								particlepair_types.getBij(ind1)/pow(rij,6);

						energy += ljpot;  // units in kcal/mol
					} else if(pmf_lookup_table==1){
						int n_lookup=	particlepair_types.getNLookup(ind1);
						if (rij < particlepair_types.getRLookup(ind1,0)) {
							// give a very large value, if the radial distance between the two ions is
							// less than the lowest value found in the lookup table.
							energy += 10000.0;
							return (energy);
						} else if (rij > particlepair_types.getRLookup(ind1,n_lookup-1)) {

							// interaction energy for particle pairs at distances greater than the maximum distance recorded in the lookup table
							// are set to zero.
							energy += 0;
						} else {


							int indm = int(
									floor(
											(rij - particlepair_types.getRLookup(ind1,0))
											/ (particlepair_types.getRLookup(ind1,1)
													- particlepair_types.getRLookup(ind1,0))));
							if (indm != n_lookup - 1) {
								double x = (rij - particlepair_types.getRLookup(ind1,indm))
																																																																																																							/ (particlepair_types.getRLookup(ind1,1)
																																																																																																									- particlepair_types.getRLookup(ind1,0));
								energy += cspline_interp(
										particlepair_types.getSplineCoeff(ind1,indm).cubic, x);

							} else if (indm == n_lookup - 1) {

								energy += particlepair_types.getPMFLookup(ind1,n_lookup - 1);
							}
						}
					}
				}
			}

		}
		//std::cout << "[energy.cpp(calcPPEnergy_SameType)]: total energy(kcal/mol) = " << energy << std::endl;
		//break;
	}
	//default:
else{	std::cout <<"[energy.cpp (calcPPEnergy_SameType)]: Wrong value set for state_change. " << std::endl;

		//break;

	}
	return (energy);
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * @brief Calculates the total interaction energy between particle of different types
 *
 * @param[in]	state_change(bool)				True or False option to indicate whether the state of the system has changed or not
 * @param[in]	use_spm(bool)					True or False option to indicate whether to use SPM or not
 * @param[in]	solvent_ion_att(bool)			True or False option to indicate whether there is an attractive square well potential between ion and solvent
 * @param[in]	particle_types(CParticleType_t)	Class instance containing parameters of each particle type
 * @param[in]	particle_state(CParticleState_t)	Class instance containing particle state data
 * @param[in]	particlepair_types(CParticlePairType_t)	Class instance containing parameters of each particle pair type
 * @param[in]	cavity_grid(CCavityGrid_t) 		Class instance containing parameters of the cavity grid data
 * @param[in]	box(CBox_t)				Class instance containing parameters of the simulation box
 * @param[in]	solvdiel(double)		Solvent dielectric constant
 * @param[in]	ptclindex(int) 		Vector element index of particle state vector
 * @param[in]	ptcltype(int)		Vector element index of particle type vector
 *
 * @return 	Energy(double) in kcal/mol`
 *
 */
double calcParticlePairEnergy_DifferentType(bool state_change,bool use_spm,bool solvent_ion_att,const
		CParticleType_t &particle_types,const CParticleState_t	&particle_state, const
		CParticlePairType_t &particlepair_types,const CCavityGrid_t &cavity_grid,const CBox_t &box,
		double solvdiel, int ptclindex,int ptcltype){

	double energy = 0.0;

	double xlen = box.getBoxXLen();
	double ylen = box.getBoxYLen();
	double zlen = box.getBoxZLen();

	int nparticle_types = particle_types.getNumParticleTypes();
	int nion_types = nparticle_types;

	int solventtype_index;

	if(use_spm == true){
		--nion_types;  // Note:Assuming that the solvent is the last item in the particle type list.
		solventtype_index = nparticle_types-1;
	}
	//switch (state_change){
	//case false:
	/*Calculate interaction energy between all particles of different types */
	if(state_change==false){
		for (int ind1=0; ind1<nion_types-1; ind1++){
			int ind1_numtype = particle_types.getNum(ind1);
			std::string label1= particle_types.getLabel(ind1);

			for (int ind2 =0; ind2 < ind1_numtype; ind2++){
				int i = particle_state.getParticleIndicesByCount(ind1,ind2);

				double xi = particle_state.getParticlePositionX(i);
				double yi = particle_state.getParticlePositionY(i);
				double zi = particle_state.getParticlePositionZ(i);

				for (int ind3 =ind1+1; ind3 < nion_types; ind3++){
					std::string label2= label1;
					label2.append("_");
					label2.append(particle_types.getLabel(ind3));
					// find index of particle pair
					int pair_index = particlepair_types.getIndexFromLabel(label2);
					bool clmb_potential = particlepair_types.getCoulombPotential(pair_index);
					bool lj_potential = particlepair_types.getLennardJonesPotential(pair_index);
					bool pmf_lookup_table = particlepair_types.getPMFLookTable(ind1);

					int ind3_numtype = particle_types.getNum(ind3);
					for (int ind4 =0; ind4 < ind3_numtype; ind4++){

						int j = particle_state.getParticleIndicesByCount(ind3,ind4);

						double xj = particle_state.getParticlePositionX(j);
						double	yj = particle_state.getParticlePositionY(j);
						double	zj = particle_state.getParticlePositionZ(j);

						// Find the nearest image
						if (xj < xi - xlen * 0.5) {
							xj +=xlen;
						} else if (xj > xi + xlen * 0.5) {
							xj -= xlen;
						}
						if (yj < yi - ylen * 0.5) {
							yj += ylen;
						} else if (yj > yi + ylen * 0.5) {
							yj  -= ylen;
						}
						if (zj < zi - zlen * 0.5) {
							zj += zlen;
						} else if (zj > zi + zlen * 0.5) {
							zj -= zlen;
						}

						double xij = xi - xj;
						double yij = yi - yj;
						double zij = zi - zj;

						double rij = sqrt(xij * xij + yij * yij + zij * zij);
						double zvi =particle_state.getParticleCharge(i);
						double	zvj = particle_state.getParticleCharge(j);


						//	std::cout << "label =" << label2 << ", rij = " << rij << std::endl;

						if(clmb_potential==1 && lj_potential ==1){
							// Assumption: pairs of similar particles are listed in the
							//same order as the particle types

							double ljpot = particlepair_types.getAij(pair_index)/pow(rij,12) -
									particlepair_types.getBij(pair_index)/pow(rij,6);

							double clmbpot = 332.06364 * zvi * zvj / (solvdiel * rij); // units in kcal/mol

							energy += ljpot;
							energy += clmbpot;

						}else if(clmb_potential==1){
							double clmbpot = 332.06364 * zvi * zvj / (solvdiel * rij); // units in kcal/mol
							//	energy = energy + clmbpot;	// units in kcal/mol
							//	std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: all-to-all, coulomb potential, energy]"
							//														<< energy << ", between pair types " << particlepair_types[pair_index].label12 << std::endl;
							/*std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: all-to-all, coulomb potential, energy] rij = " <<
									rij << ", coulomb energy (kcal/mol) = " << clmbpot << ", pair label " << particlepair_types[pair_index].label12
									<< " between particle indices "  << i << " and " << j <<  std::endl;*/
							energy += clmbpot;
						}
						else if(lj_potential ==1){

							double ljpot = particlepair_types.getAij(pair_index)/pow(rij,12) -
									particlepair_types.getBij(pair_index)/pow(rij,6);

							energy += ljpot;
						}else if(pmf_lookup_table==1){
							int n_lookup=	particlepair_types.getNLookup(pair_index);
							if (rij < particlepair_types.getRLookup(pair_index,0)) {
								// give a very large value, if the radial distance between the two ions is
								// less than the lowest value found in the lookup table.
								energy += 10000.0;
								return (energy);
							} else if (rij > particlepair_types.getRLookup(pair_index,n_lookup-1)) {

								// interaction energy for particle pairs at distances greater than the maximum distance recorded in the lookup table
								// are set to zero.
								energy += 0;
							} else {
								int indm = int(
										floor(
												(rij - particlepair_types.getRLookup(pair_index,0))
												/ (particlepair_types.getRLookup(pair_index,1)
														- particlepair_types.getRLookup(pair_index,0))));
								if (indm != n_lookup - 1) {
									double x = (rij - particlepair_types.getRLookup(pair_index,indm))
																																																																																																											/ (particlepair_types.getRLookup(pair_index,1)
																																																																																																													- particlepair_types.getRLookup(pair_index,0));
									energy += cspline_interp(
											particlepair_types.getSplineCoeff(pair_index,indm).cubic, x);

								} else if (indm == n_lookup - 1) {

									energy += particlepair_types.getPMFLookup(pair_index,n_lookup - 1);
								}
							}
						}

					}
				}



			}
		} // ind1

		// Lennard Jones interaction between ion and solvent
		if(use_spm==true){
			for (int ind1=0; ind1<nion_types; ind1++){
				int ind1_numtype = particle_types.getNum(ind1);
				bool lj_potential_w;
				std::string label1w= particle_types.getLabel(ind1);
				//	std::cout << "[energy.cpp(all ion-water pairs)]: selected particle type" << label1w <<
				//		", total number = " << ind1_numtype <<std::endl;
				label1w.append("_Water");

				int ion_water_pair_index;

				int lj_cutoff ;

				// Lennard-Jones interaction between ions and water is true or false

				// find index of particle pair
				ion_water_pair_index = particlepair_types.getIndexFromLabel(label1w);

				lj_potential_w = particlepair_types.getLennardJonesPotential(ion_water_pair_index);

				if(lj_potential_w==true){
					// Lennard-Jones interaction between ion and water
					lj_cutoff = particlepair_types.getLJCutoff(label1w);
					// get cavity grid cell index of the ion
					for (int ind2 =0; ind2 < ind1_numtype; ind2++){
						int i = particle_state.getParticleIndicesByCount(ind1,ind2);

						double xi = particle_state.getParticlePositionX(i);
						double yi = particle_state.getParticlePositionY(i);
						double zi = particle_state.getParticlePositionZ(i);

						int ptcl_cellindex1D = particle_state.getParticleGridCellIndex1D(i);
						// get the LJ cell index of the cavity grid cell index
						int lj_ptcl_cellindex1D = cavity_grid.getLJCellindexByCavityCellIndex(ptcl_cellindex1D);

						int num_lj_neigborcells = cavity_grid.getLJCellindexOffsetNindices();

						//	std::cout << "[energy.cpp]: pair label" << label1w  << ", ptcl_cellindex1D = " << ptcl_cellindex1D << ", LJ ptcl index = " <<
						//											lj_ptcl_cellindex1D << std::endl;


						std::vector<int> lj_neighbor_cellindices = cavity_grid.getLJNeighborCellIndices(lj_ptcl_cellindex1D);
						for (int ljcell=0; ljcell < num_lj_neigborcells; ljcell++){
							// get cavity grid cell indices that are occupied
							int lj_cellindex1D = lj_neighbor_cellindices[ljcell];
							struct_ljcell_cavityindices lj_cavity_cellindices =
									particle_types.getLJCellCavityCellIndices(solventtype_index, lj_cellindex1D);
							int num_cavity_cells = lj_cavity_cellindices.num;
							//	std::cout << "...lj_cellindex1D of other neighboring particle << " << label1w << ", " <<
							//				lj_cellindex1D << ", num cavity cells = " << num_cavity_cells << std::endl;

							if(num_cavity_cells > 0){
								for (int lj_cavitycell=0; lj_cavitycell < num_cavity_cells; lj_cavitycell++){
									int lj_cavity_cellindex = lj_cavity_cellindices.cavity_1Dcellindices[lj_cavitycell];
									int solvent_index = particle_state.getParticleIndicesByCellindex(solventtype_index,lj_cavity_cellindex);

									// get (x,y,z) position

									double xj = particle_state.getParticlePositionX(solvent_index);
									double	yj = particle_state.getParticlePositionY(solvent_index);
									double	zj = particle_state.getParticlePositionZ(solvent_index);
									if(particle_state.getParticle(solvent_index).label!="Water"){
										std::cerr << "[energy.cpp (calcParticlePairEnergy_DifferentType)]: particle type is not Water." << std::endl;
									}

									// Find the nearest image
									if (xj < xi - xlen * 0.5) {
										xj +=xlen;
									} else if (xj > xi + xlen * 0.5) {
										xj -= xlen;
									}
									if (yj < yi - ylen * 0.5) {
										yj += ylen;
									} else if (yj > yi + ylen * 0.5) {
										yj  -= ylen;
									}
									if (zj < zi - zlen * 0.5) {
										zj += zlen;
									} else if (zj > zi + zlen * 0.5) {
										zj -= zlen;
									}

									double xij = xi - xj;
									double yij = yi - yj;
									double zij = zi - zj;

									double rij = sqrt(xij * xij + yij * yij + zij * zij);
									//std::cout << "label =" << label1w << ", rij = " << rij << std::endl;

									if(rij <= lj_cutoff){
										double ljpot = particlepair_types.getAij(ion_water_pair_index)/pow(rij,12) -
												particlepair_types.getBij(ion_water_pair_index)/pow(rij,6);

										//std::cout << "rij = " << rij << ", ljpot = " << ljpot << std::endl;
										energy += ljpot;
									}
								}
							}

						}
					}
				}
			}

		}

		//std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: total energy(kcal/mol) between non-identical ions = " << energy << std::endl
		// between solvent-ion pairs
		if(use_spm==true && solvent_ion_att==true){
			//std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: Calculating interaction energy between ions and solvent." << std::endl;
			std::string label2 = particle_types.getLabel(nparticle_types-1);   // solvent label
			for (int ind1=0; ind1<nparticle_types-1; ind1++){  // only ions; hence, nparticle_types-1
				std::string label12= particle_types.getLabel(ind1);   // ion label
				label12.append("_");
				label12.append(label2);
				// find index of solvent-ion pair
				int pair_index = particlepair_types.getIndexFromLabel(label12);
				if(particlepair_types.getSqrWellPotential(pair_index)==1){
					// find max value of well width

					double hs_cutoff = particlepair_types.getHardSphereCutoff(label12);
					double sqrwell_cutoff = (particlepair_types.getSqrwell_wdth2contact_ratio(pair_index)+1.0)*hs_cutoff;

					int numtype_ind1 = particle_types.getNum(ind1);

					// find grid cells of solvent spheres that are within the square well width for each ion type
					for (int ind2=0;ind2<numtype_ind1;ind2++){
						int i = particle_state.getParticleIndicesByCount(ind1,ind2);

						int lix = particle_state.getParticleGridCellIndexX(i);
						int liy = particle_state.getParticleGridCellIndexY(i);
						int liz = particle_state.getParticleGridCellIndexZ(i);

						struct_cellindex_offset cellindex_offset = particlepair_types.getSqrWellWidthCellIndexOffsetMatrix(label12);
						int ncx = cellindex_offset.ncx;
						int ncy = cellindex_offset.ncy;
						int ncz = cellindex_offset.ncz;
						int Ncx = cellindex_offset.Ncx;
						int Ncy = cellindex_offset.Ncy;
						int Ncz = cellindex_offset.Ncz;

						int ljx,ljy,ljz;

						for (int ind3=0;ind3<2*ncx+1;ind3++){
							ljx = lix + cellindex_offset.delxj[ind3] + cellindex_offset.bcmap_x[lix][ind3];


							for (int ind4 = 0; ind4<2*ncy+1; ind4++){
								ljy = liy + cellindex_offset.delyj[ind4] + cellindex_offset.bcmap_y[liy][ind4];

								for (int ind5 = 0; ind5 < 2*ncz+1; ind5++){
									ljz = liz + cellindex_offset.delzj[ind5] + cellindex_offset.bcmap_z[liz][ind5];

									int cellindex_1D =  cavity_grid.getCellIndex1D(ljx,ljy,ljz);
									int solvent_index = particle_state.getParticleIndicesByCellindex(nparticle_types-1,cellindex_1D);

									if (solvent_index!=-1){

										double xi = particle_state.getParticlePositionX(i);
										double yi = particle_state.getParticlePositionY(i);
										double zi = particle_state.getParticlePositionZ(i);

										double xj = particle_state.getParticlePositionX(solvent_index);
										double yj = particle_state.getParticlePositionY(solvent_index);
										double zj = particle_state.getParticlePositionZ(solvent_index);


										// Find the nearest image
										if (xj < xi - xlen * 0.5) {
											xj +=xlen;
										} else if (xj > xi + xlen * 0.5) {
											xj -= xlen;
										}
										if (yj < yi - ylen * 0.5) {
											yj += ylen;
										} else if (yj > yi + ylen * 0.5) {
											yj  -= ylen;
										}
										if (zj < zi - zlen * 0.5) {
											zj += zlen;
										} else if (zj > zi + zlen * 0.5) {
											zj -= zlen;
										}

										double xij = xi - xj;
										double yij = yi - yj;
										double zij = zi - zj;

										double rij = sqrt(xij * xij + yij * yij + zij * zij);


										if(rij <= sqrwell_cutoff && rij >= hs_cutoff){

											energy += particlepair_types.getSqrWellDepth(pair_index);
										}
									}
								}
							}
						}

					}
				}

			}
		}

		//std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: total energy(kcal/mol) = " << energy << std::endl;
		//break;
	} else if(state_change==true){
		/*Calculate interaction energy between 1 particle and all dissimilar particles */
	//case true:

	//{
		// between ion pairs
		int ind1, i;
		ind1 = ptcltype-1;   // get particle type index
		i = ptclindex;   // get particle index in particles vector

		double xi = particle_state.getParticlePositionX(i);
		double yi = particle_state.getParticlePositionY(i);
		double zi = particle_state.getParticlePositionZ(i);


		if(ptcltype <=nion_types) { // the particle is an ion (not a solvent)

			for (int ind2 =0; ind2 < nion_types; ind2++){
				if(ind2!=ind1){
					std::string label12= particle_types.getLabel(ind1);
					label12.append("_");
					label12.append(particle_types.getLabel(ind2));
					// find index of particle pair


					int pair_index = particlepair_types.getIndexFromLabel(label12);
					int ind2_numtype = particle_types.getNum(ind2);

					bool clmb_potential = particlepair_types.getCoulombPotential(pair_index);
					bool lj_potential = particlepair_types.getLennardJonesPotential(pair_index);
					bool pmf_lookup_table = particlepair_types.getPMFLookTable(ind1);


					for (int ind3 =0; ind3 <ind2_numtype ; ind3++){

						int j = particle_state.getParticleIndicesByCount(ind2,ind3);

						double xj = particle_state.getParticlePositionX(j);
						double	yj = particle_state.getParticlePositionY(j);
						double	zj = particle_state.getParticlePositionZ(j);

						// Find the nearest image
						if (xj < xi - xlen * 0.5) {
							xj += xlen;
						} else if (xj > xi + xlen * 0.5) {
							xj -= xlen;
						}
						if (yj < yi - ylen * 0.5) {
							yj += ylen;
						} else if (yj > yi + ylen * 0.5) {
							yj -= ylen;
						}
						if (zj < zi - zlen * 0.5) {
							zj += zlen;
						} else if (zj > zi + zlen * 0.5) {
							zj -= zlen;
						}

						double xij = xi - xj;
						double yij = yi - yj;
						double zij = zi - zj;

						double rij = sqrt(xij * xij + yij * yij + zij * zij);

						double zvi =particle_state.getParticleCharge(i);
						double	zvj = particle_state.getParticleCharge(j);

						//std::cout << "label =" << label12 << ", rij = " << rij << std::endl;

						if(clmb_potential==1 && lj_potential ==1){
							// Assumption: pairs of similar particles are listed in the same order as the particle types

							double ljpot = particlepair_types.getAij(pair_index)/pow(rij,12) -
									particlepair_types.getBij(pair_index)/pow(rij,6);

							double clmbpot = 332.06364 * zvi *	 zvj / (solvdiel * rij); // units in kcal/mol
							//energy = energy + ljpot + clmbpot;	// units in kcal/mol
							energy += ljpot;
							energy += clmbpot;

						}else if(clmb_potential==1){
							double clmbpot = 332.06364 * zvi * zvj / (solvdiel * rij); // units in kcal/mol
							//energy = energy + clmbpot;	// units in kcal/mol
							energy += clmbpot;
							//	std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: 1-to-many, coulomb potential, energy]"
							//		<< energy << ", between pair types " << particlepair_types[pair_index].label12 << std::endl;

							/*	std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: one-to-many, coulomb potential, energy] rij = " <<
																rij << ", coulomb energy (kcal/mol) = " << clmbpot << ", pair label "
																<< particlepair_types[pair_index].label12 <<
																" between particle indices "  << i << " and " << j << std::endl;
							 */
						}
						else if(lj_potential ==1){

							double ljpot = particlepair_types.getAij(pair_index)/pow(rij,12) -
									particlepair_types.getBij(pair_index)/pow(rij,6);

							//energy = energy + ljpot;	// units in kcal/mol
							energy += ljpot;
						}else if(pmf_lookup_table==1){
							int n_lookup=	particlepair_types.getNLookup(pair_index);
							if (rij < particlepair_types.getRLookup(pair_index,0)) {
								// give a very large value, if the radial distance between the two ions is
								// less than the lowest value found in the lookup table.
								energy += 10000.0;
								return (energy);
							} else if (rij > particlepair_types.getRLookup(pair_index,n_lookup-1)) {

								// interaction energy for particle pairs at distances greater than the maximum distance recorded in the lookup table
								// are set to zero.
								energy += 0;
							} else {


								int indm = int(
										floor(
												(rij - particlepair_types.getRLookup(pair_index,0))
												/ (particlepair_types.getRLookup(pair_index,1)
														- particlepair_types.getRLookup(pair_index,0))));
								if (indm != n_lookup - 1) {
									double x = (rij - particlepair_types.getRLookup(pair_index,indm))
																																																																																																											/ (particlepair_types.getRLookup(pair_index,1)
																																																																																																													- particlepair_types.getRLookup(pair_index,0));
									energy += cspline_interp(
											particlepair_types.getSplineCoeff(pair_index,indm).cubic, x);


									//	energy =  energy + cspline_interp(ionpair_pmf[indk].spline_coeff[indm].cubic,x);
								} else if (indm == n_lookup - 1) {
									energy += particlepair_types.getPMFLookup(pair_index,n_lookup - 1);
								}
							}
						}

					}
				}
			}

			if(use_spm==true){
				// Lennard-Jones interaction between ions and water is true or false

				std::string label1w= particle_types.getLabel(ind1);
				//	std::cout << "[energy.cpp (ion vs all neighbor water]: selected ion type = " << label1w << std::endl;

				label1w.append("_Water");

				// find index of particle pair
				double ion_water_pair_index = particlepair_types.getIndexFromLabel(label1w);
				double lj_potential_w = particlepair_types.getLennardJonesPotential(ion_water_pair_index);

				if(lj_potential_w==true){
					// Lennard-Jones interaction between ion and water
					double lj_cutoff = particlepair_types.getLJCutoff(label1w);

					// get cavity grid cell index of the ion
					int ptcl_cellindex1D = particle_state.getParticleGridCellIndex1D(i);
					// get the LJ cell index of the cavity grid cell index
					int lj_ptcl_cellindex1D = cavity_grid.getLJCellindexByCavityCellIndex(ptcl_cellindex1D);
					int num_lj_neighborcells = cavity_grid.getLJCellindexOffsetNindices();

					//	std::cout << "[energy.cpp]: pair label" << label1w  << ", ptcl_cellindex1D = " << ptcl_cellindex1D << ", LJ ptcl index = " <<
					//															lj_ptcl_cellindex1D << std::endl;

					std::vector<int> lj_neighbor_cellindices = cavity_grid.getLJNeighborCellIndices(lj_ptcl_cellindex1D);
					for (int ljcell=0; ljcell < num_lj_neighborcells; ljcell++){
						// get cavity grid cell indices that are occupied
						int lj_cellindex1D = lj_neighbor_cellindices[ljcell];
						struct_ljcell_cavityindices lj_cavity_cellindices =
								particle_types.getLJCellCavityCellIndices(solventtype_index, lj_cellindex1D);
						int num_cavity_cells = lj_cavity_cellindices.num;
						//	std::cout << "...lj_cellindex1D of other neighboring particle << " << label1w << ", " <<
						//			 lj_cellindex1D << ", num cavity cells = " << num_cavity_cells << std::endl;

						if(num_cavity_cells > 0){
							for (int lj_cavitycell=0; lj_cavitycell < num_cavity_cells; lj_cavitycell++){
								int lj_cavity_cellindex = lj_cavity_cellindices.cavity_1Dcellindices[lj_cavitycell];
								int solvent_index = particle_state.getParticleIndicesByCellindex(solventtype_index,lj_cavity_cellindex);

								// get (x,y,z) position

								double xj = particle_state.getParticlePositionX(solvent_index);
								double	yj = particle_state.getParticlePositionY(solvent_index);
								double	zj = particle_state.getParticlePositionZ(solvent_index);
								if(particle_state.getParticle(solvent_index).label!="Water"){
									std::cerr << "[energy.cpp (calcParticlePairEnergy_DifferentType)]: particle type is not Water." << std::endl;
								}

								// Find the nearest image
								if (xj < xi - xlen * 0.5) {
									xj +=xlen;
								} else if (xj > xi + xlen * 0.5) {
									xj -= xlen;
								}
								if (yj < yi - ylen * 0.5) {
									yj += ylen;
								} else if (yj > yi + ylen * 0.5) {
									yj  -= ylen;
								}
								if (zj < zi - zlen * 0.5) {
									zj += zlen;
								} else if (zj > zi + zlen * 0.5) {
									zj -= zlen;
								}

								double xij = xi - xj;
								double yij = yi - yj;
								double zij = zi - zj;

								double rij = sqrt(xij * xij + yij * yij + zij * zij);

								//std::cout << "label =" << label1w << ", rij = " << rij << std::endl;

								if(rij <= lj_cutoff){
									double ljpot = particlepair_types.getAij(ion_water_pair_index)/pow(rij,12) -
											particlepair_types.getBij(ion_water_pair_index)/pow(rij,6);

									//	std::cout << "rij = " << rij << ", ljpot = " << ljpot << std::endl;
									energy += ljpot;
								}
							}
						}

					}

				}
			}
			// attractive square well potentials between solvent-ion pairs
			if(use_spm==true && solvent_ion_att==true){
				std::string label12= particle_types.getLabel(ind1);
				label12.append("_");
				label12.append(particle_types.getLabel(nparticle_types-1));

				int pair_index = particlepair_types.getIndexFromLabel(label12);
				if(particlepair_types.getSqrWellPotential(pair_index) ==1){
					// find max value of well width

					//		double hs_cutoff = particlepair_types.getRadius1(pair_index) + particlepair_types.getRadius2(pair_index);
					double hs_cutoff = particlepair_types.getHardSphereCutoff(label12);
					double sqrwell_cutoff = (particlepair_types.getSqrwell_wdth2contact_ratio(pair_index)+1.0)*hs_cutoff;

					// find grid cells of solvent spheres that are within the square well width

					int lix = particle_state.getParticleGridCellIndexX(i);
					int liy = particle_state.getParticleGridCellIndexY(i);
					int liz = particle_state.getParticleGridCellIndexZ(i);

					struct_cellindex_offset cellindex_offset = particlepair_types.getSqrWellWidthCellIndexOffsetMatrix(label12);

					int ncx = cellindex_offset.ncx;
					int ncy = cellindex_offset.ncy;
					int ncz = cellindex_offset.ncz;
					int Ncx = cellindex_offset.Ncx;
					int Ncy = cellindex_offset.Ncy;
					int Ncz = cellindex_offset.Ncz;

					int ljx,ljy,ljz;

					for (int ind3=0;ind3<2*ncx+1;ind3++){
						ljx = lix + cellindex_offset.delxj[ind3] + cellindex_offset.bcmap_x[lix][ind3];

						for (int ind4 = 0; ind4<2*ncy+1; ind4++){
							ljy = liy + cellindex_offset.delyj[ind4] + cellindex_offset.bcmap_y[liy][ind4];

							for (int ind5 = 0; ind5 < 2*ncz+1; ind5++){
								ljz = liz + cellindex_offset.delzj[ind5] + cellindex_offset.bcmap_z[liz][ind5];

								int cellindex_1D = cavity_grid.getCellIndex1D(ljx,ljy,ljz);
								int solvent_index = particle_state.getParticleIndicesByCellindex(nparticle_types-1,cellindex_1D);

								if (solvent_index!=-1){
									double xi = particle_state.getParticlePositionX(i);
									double yi = particle_state.getParticlePositionY(i);
									double zi = particle_state.getParticlePositionZ(i);

									double xj = particle_state.getParticlePositionX(solvent_index);
									double yj = particle_state.getParticlePositionY(solvent_index);
									double zj = particle_state.getParticlePositionZ(solvent_index);


									// Find the nearest image
									if (xj < xi - xlen * 0.5) {
										xj +=xlen;
									} else if (xj > xi + xlen * 0.5) {
										xj -= xlen;
									}
									if (yj < yi - ylen * 0.5) {
										yj += ylen;
									} else if (yj > yi + ylen * 0.5) {
										yj  -= ylen;
									}
									if (zj < zi - zlen * 0.5) {
										zj += zlen;
									} else if (zj > zi + zlen * 0.5) {
										zj -= zlen;
									}

									double xij = xi - xj;
									double yij = yi - yj;
									double zij = zi - zj;

									double rij = sqrt(xij * xij + yij * yij + zij * zij);

									if(rij <= sqrwell_cutoff && rij >= hs_cutoff){
										//energy = energy + particlepair_types[pair_index].sqrwell_depth;

										energy += particlepair_types.getSqrWellDepth(pair_index);

										//			std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: ion-solvent, square well potential, energy]"
										//					<< energy << ", between pair types " << particlepair_types[pair_index].label12 << std::endl;

									}
								}
							}
						}
					}

				}
			}


		}
		if(ptcltype == nion_types + 1 && use_spm==true){ // the particle is a solvent

			// Lennard-Jones interaction between water and ions is true or false

			int num_lj_cells = cavity_grid.getLJGridNumCells();
			//std::cout << "selected particle is Water" << ", num_lj_cells = " << num_lj_cells << std::endl;

			if (num_lj_cells > 0){
				for (int ind2=0; ind2 <nion_types;ind2++){
					std::string label1w= particle_types.getLabel(ind2);
					label1w.append("_Water");
					// find index of particle pair
					double ion_water_pair_index = particlepair_types.getIndexFromLabel(label1w);
					double lj_potential_w = particlepair_types.getLennardJonesPotential(ion_water_pair_index);


					if(lj_potential_w==true){
						// Lennard-Jones interaction between ion and water
						double lj_cutoff = particlepair_types.getLJCutoff(label1w);

						// get cavity grid cell index of the ion
						int ptcl_cellindex1D = particle_state.getParticleGridCellIndex1D(i);
						// get the LJ cell index of the cavity grid cell index
						int lj_ptcl_cellindex1D = cavity_grid.getLJCellindexByCavityCellIndex(ptcl_cellindex1D);
						int num_lj_neighborcells = cavity_grid.getLJCellindexOffsetNindices();

						//	std::cout << "[energy.cpp]: pair label" << label1w  << ", ptcl_cellindex1D = " << ptcl_cellindex1D << ", LJ ptcl index = " <<
						//															lj_ptcl_cellindex1D << std::endl;


						std::vector<int> lj_neighbor_cellindices = cavity_grid.getLJNeighborCellIndices(lj_ptcl_cellindex1D);
						for (int ljcell=0; ljcell < num_lj_neighborcells; ljcell++){
							// get cavity grid cell indices that are occupied
							int lj_cellindex1D = lj_neighbor_cellindices[ljcell];
							struct_ljcell_cavityindices lj_cavity_cellindices =
									particle_types.getLJCellCavityCellIndices(ind2, lj_cellindex1D);
							int num_cavity_cells = lj_cavity_cellindices.num;
							//			std::cout << "...lj_cellindex1D of other neighboring particle << " << label1w <<", " <<
							//					 lj_cellindex1D << ", num cavity cells = " << num_cavity_cells << std::endl;

							if (num_cavity_cells > 0){
								for (int lj_cavitycell=0; lj_cavitycell < num_cavity_cells; lj_cavitycell++){
									int lj_cavity_cellindex = lj_cavity_cellindices.cavity_1Dcellindices[lj_cavitycell];
									int ion_index = particle_state.getParticleIndicesByCellindex(ind2,lj_cavity_cellindex);

									// get (x,y,z) position

									double xj = particle_state.getParticlePositionX(ion_index);
									double	yj = particle_state.getParticlePositionY(ion_index);
									double	zj = particle_state.getParticlePositionZ(ion_index);
									if(particle_state.getParticle(ion_index).label!=particle_types.getLabel(ind2)){
										std::cerr << "[energy.cpp (calcParticlePairEnergy_DifferentType)]: ion types do not match." << std::endl;
									}

									// Find the nearest image
									if (xj < xi - xlen * 0.5) {
										xj +=xlen;
									} else if (xj > xi + xlen * 0.5) {
										xj -= xlen;
									}
									if (yj < yi - ylen * 0.5) {
										yj += ylen;
									} else if (yj > yi + ylen * 0.5) {
										yj  -= ylen;
									}
									if (zj < zi - zlen * 0.5) {
										zj += zlen;
									} else if (zj > zi + zlen * 0.5) {
										zj -= zlen;
									}

									double xij = xi - xj;
									double yij = yi - yj;
									double zij = zi - zj;

									double rij = sqrt(xij * xij + yij * yij + zij * zij);


									if(rij <= lj_cutoff){
										double ljpot = particlepair_types.getAij(ion_water_pair_index)/pow(rij,12) -
												particlepair_types.getBij(ion_water_pair_index)/pow(rij,6);

										//			std::cout << "rij = " << rij << ", ljpot = " << ljpot << std::endl;

										energy += ljpot;
									}
								}
							}
						}

					}
				}
			}

			// SOLVENT ION SQUARE WELL POTENTIAL
			if(solvent_ion_att==true){ // solvent ion square well attractive potential
				for (int ind2 =0; ind2 < nion_types; ind2++){

					std::string label12= particle_types.getLabel(ind1);
					label12.append("_");
					label12.append(particle_types.getLabel(ind2));
					// find index of particle pair

					int pair_index = particlepair_types.getIndexFromLabel(label12);
					if(particlepair_types.getSqrWellPotential(pair_index) ==1){
						// find max value of well width

						//double hs_cutoff = particlepair_types.getRadius1(pair_index) + particlepair_types.getRadius2(pair_index);
						double hs_cutoff = particlepair_types.getHardSphereCutoff(label12);

						double sqrwell_cutoff = (particlepair_types.getSqrwell_wdth2contact_ratio(pair_index)+1.0)*hs_cutoff;
						// find grid cells of solvent spheres that are within the square well width


						int lix = particle_state.getParticleGridCellIndexX(i);
						int liy = particle_state.getParticleGridCellIndexY(i);
						int liz = particle_state.getParticleGridCellIndexZ(i);

						struct_cellindex_offset cellindex_offset = particlepair_types.getSqrWellWidthCellIndexOffsetMatrix(label12);
						int ncx = cellindex_offset.ncx;
						int ncy = cellindex_offset.ncy;
						int ncz = cellindex_offset.ncz;
						//	int Ncx = cellindex_offset.Ncx;
						//	int Ncy = cellindex_offset.Ncy;
						//	int Ncz = cellindex_offset.Ncz;

						int ljx,ljy,ljz;

						for (int ind3=0;ind3<2*ncx+1;ind3++){
							ljx = lix + cellindex_offset.delxj[ind3] + cellindex_offset.bcmap_x[lix][ind3];


							for (int ind4 = 0; ind4<2*ncy+1; ind4++){
								ljy = liy + cellindex_offset.delyj[ind4] + cellindex_offset.bcmap_y[liy][ind4];

								for (int ind5 = 0; ind5 < 2*ncz+1; ind5++){
									ljz = liz + cellindex_offset.delzj[ind5] + cellindex_offset.bcmap_z[liz][ind5];

									int cellindex_1D =  cavity_grid.getCellIndex1D(ljx,ljy,ljz);
									int ion_index = particle_state.getParticleIndicesByCellindex(ind2,cellindex_1D);

									if (ion_index!=-1){


										double xi = particle_state.getParticlePositionX(i);
										double yi = particle_state.getParticlePositionY(i);
										double zi = particle_state.getParticlePositionZ(i);



										double xj = particle_state.getParticlePositionX(ion_index);
										double	yj = particle_state.getParticlePositionY(ion_index);
										double	zj = particle_state.getParticlePositionZ(ion_index);


										// Find the nearest image
										if (xj < xi - xlen * 0.5) {
											xj +=xlen;
										} else if (xj > xi + xlen * 0.5) {
											xj -= xlen;
										}
										if (yj < yi - ylen * 0.5) {
											yj += ylen;
										} else if (yj > yi + ylen * 0.5) {
											yj  -= ylen;
										}
										if (zj < zi - zlen * 0.5) {
											zj += zlen;
										} else if (zj > zi + zlen * 0.5) {
											zj -= zlen;
										}

										double xij = xi - xj;
										double yij = yi - yj;
										double zij = zi - zj;

										double rij = sqrt(xij * xij + yij * yij + zij * zij);

										if(rij <= sqrwell_cutoff && rij >= hs_cutoff){
											//energy = energy + particlepair_types[pair_index].sqrwell_depth;
											energy += particlepair_types.getSqrWellDepth(pair_index);
											//			std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: solvent-ion, square well potential, energy]"
											//					<< energy << ", between pair types " << particlepair_types[pair_index].label12 << std::endl;

										}
									}
								}
							}
						}

					}

				}

			}

		}
		//	std::cout << "[energy.cpp(calcPPEnergy_DifferentType)]: total energy(kcal/mol) = " << energy << std::endl;
		//break;
	}else{

		std::cout <<"[energy.cpp (calcPPEnergy_SameType)]: Wrong value set for state_change. Should be true or false. " << std::endl;

	}

	return (energy);
}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Calculates total interaction energy between particles
 *
 * @param[in]	state_change(bool)			True or False option to indicate whether the state of the system has changed or not
 * @param[in]	use_spm(bool)				True or False option to indicate whether to use SPM or not
 * @param[in]	use_mwater(bool)			True or False option to indicate whether to use MWater model
 * @param[in]	solvent_ion_att(bool)		True or False option to indicate whether there is an attractive square well potential between ion and solvent
 * @param[in]	solvdiel(double)			Solvent dielectric constant
 * @param[in]	particle_state(CParticleState_t) 	Class instance containing particle state data
 * @param[in]	particle_types(CParticleType_t)	Class instance containing parameters of each particle type
 * @param[in]	mwater(CMWater)				Class instance containing parameters of MWater model
 * @param[in]	particlepair_types(CParticlePairType_t)	Class instance containing parameters of each particle pair type
 * @param[in]	cavity_grid(CCavityGrid_t)		Class instance containing parameters of the cavity grid data
 * @param[in]	box(CBox_t)			Class instance containing parameters of the simulation box
 * @param[in]	ptclindex(int)		Vector element index of particle state vector
 * @param[in]	ptcltype(int)		Vector element index of particle type vector
 *
 *
 * @return 	Energy(double) in kcal/mol
 *
 */
double calcParticlePairEnergy(bool state_change,bool use_spm,bool use_mwater,bool solvent_ion_att,
		double solvdiel,const CParticleState_t &particle_state,const CParticlePairType_t &particlepair_types,
		const CParticleType_t &particle_types, CMWater &mwater, const CCavityGrid_t & cavity_grid,
		const CBox_t &box,int ptclindex,int ptcltype){


	double energy = 0.0;

	// calculate energy between particles of the same type.
	energy  = calcParticlePairEnergy_SameType(state_change,use_spm,use_mwater,particle_types,particle_state, particlepair_types,
			mwater,cavity_grid, box,solvdiel,ptclindex, ptcltype);

	energy += calcParticlePairEnergy_DifferentType(state_change,use_spm,solvent_ion_att,particle_types,
			particle_state,particlepair_types,cavity_grid, box,solvdiel,ptclindex,ptcltype);

	return (energy);
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/**
 * @brief Calculate interaction energy between solute and particles
 *
 * @param[in]	state_change(bool)			True or False option to indicate whether the state of the system has changed or not
 * @param[in]	use_spm(bool)				True or False option to indicate whether to use SPM or not
 * @param[in]	particle_state(CParticleState_t)	Class instance containing particle state data
 * @param[in]	particle_types(CParticleType_t)	Class instance containing parameters of each particle type
 * @param[in]	solute(CSoluteModel_t)			Class instance containing parameters of the solute model
 * @param[in]	box(CBox_t)					Class instance containing parameters of the simulation box
 * @param[in]	ptclindex(int)				Vector element index of particle state vector
 *
 *
 * @return 	Energy(double) in kcal/mol
 */
double calcSoluteParticleEnergy(bool state_change,bool use_spm,const CParticleState_t &particle_state,
		const CParticleType_t &particle_types, const CSoluteModel_t &solute,const CBox_t &box,int ptclindex){
	double energy = 0.0;

	int nion_types = particle_types.getNumParticleTypes();
	if(use_spm == true){
		--nion_types;  // Note:Assuming that the solvent is the last item in the particle type list.
	}

	//switch (state_change){

	//case false:
	if (state_change==false){
		for (int ind1=0; ind1<nion_types; ind1++){
			int ind1_numtype = particle_types.getNum(ind1);
			for (int ind2 =0; ind2 < ind1_numtype; ind2++){

				double x = particle_state.getParticlePositionX(ind1);

				double y =	particle_state.getParticlePositionY(ind1);
				double z = particle_state.getParticlePositionZ(ind1);
				double potential = solute.getSoluteAllAtom_Potential(x,y,z,box) ;

				energy += potential*particle_state.getParticleCharge(ind1);
			}

		}
	//	break;

	//case true:
	}else if(state_change==true){
		double x = particle_state.getParticlePositionX(ptclindex);

		double y =	particle_state.getParticlePositionY(ptclindex);
		double z = particle_state.getParticlePositionZ(ptclindex);
		double potential = solute.getSoluteAllAtom_Potential(x,y,z,box) ;

		energy += potential*particle_state.getParticleCharge(ptclindex);


		//break;
	}

	return energy;
}
