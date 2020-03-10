/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_patchyparticlerect.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

//NR
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"
#include <iostream>
#include <iomanip>
#define pi_cst 3.14159265359

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairPatchyparticlerect::PairPatchyparticlerect(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairPatchyparticlerect::~PairPatchyparticlerect()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    
    // NR
    int n = atom->ntypes;
    for (int i = 1; i <= n; i++)
    {
      for (int j = 1; j <= n; j++)
      {
	if(patchInter[i][j].typeP1 != NULL)
	{
	  delete[] patchInter[i][j].typeP1;
	  patchInter[i][j].typeP1 = NULL;
	}
	if(patchInter[i][j].typeP2 != NULL)
	{
	  delete[] patchInter[i][j].typeP2;
	  patchInter[i][j].typeP2 = NULL;
	}
      }
    }
    memory->destroy(patchInter);
    delete[] pMorse;
    delete[] pGeometry;
    
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairPatchyparticlerect::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  //NR
  double **torque = atom->torque;
  double delti[3],deltj[3]; // torque increment;;
  double *qi;
  double ni[3];
  double *qj;
  double nj[3];
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double rij[3];
  double rnorm;
  double expTerm;
  double morseFactor;
  double morseDerivation;
  double angleTerm;
  double v1[3];
  double v2[3];
  double v3[3];
  double v4[3];
  double fi[3];


  ///// RECT /////
  double ai[3];
  double bi[3];
  double rji_a[3];
  double rji_a_norm(0.0);
  double rji_a_norm_cube(0.0);
  double rji_b[3];
  double rji_b_norm(0.0);
  double rji_b_norm_cube(0.0);
  double rji_sca_bi(0.0);
  double rji_sca_ai(0.0);
  double costi_a(0.0);
  double costi_b(0.0);
  double ti_a(0.0);
  double ti_b(0.0);
  double cosFactori_a(0.0);
  double cosFactori_b(0.0);
  double dti_dcosti_a(0.0);
  double dti_dcosti_b(0.0);
  double ni_sca_rij(0.0);
  double ni_sca_rji_a(0.0);
  double ni_sca_rji_b(0.0);

  double aj[3];
  double bj[3];
  double rij_a[3];
  double rij_a_norm(0.0);
  double rij_a_norm_cube(0.0);
  double rij_b[3];
  double rij_b_norm(0.0);
  double rij_b_norm_cube(0.0);
  double rij_sca_bj(0.0);
  double rij_sca_aj(0.0);
  double costj_a(0.0);
  double costj_b(0.0);
  double tj_a(0.0);
  double tj_b(0.0);
  double cosFactorj_a(0.0);
  double cosFactorj_b(0.0);
  double dtj_dcostj_a(0.0);
  double dtj_dcostj_b(0.0);
  double nj_sca_rij(0.0);
  double nj_sca_rij_a(0.0);
  double nj_sca_rij_b(0.0);

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      fi[0] = 0.0;
      fi[1] = 0.0;
      fi[2] = 0.0;
      delti[0] = 0.0;
      delti[1] = 0.0;
      delti[2] = 0.0;
      deltj[0] = 0.0;
      deltj[1] = 0.0;
      deltj[2] = 0.0;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      rij[0] = delx;
      rij[1] = dely;
      rij[2] = delz;
      rnorm = sqrt(rsq);

      if (rsq < cutsq[itype][jtype] && factor_lj != 0) // NR try add not calculate if factor_lj = 0
      {
        // Get particle orientation
        qi = bonus[i].quat;
        qj = bonus[j].quat;

        // Init patch energy
        double energyPatch(0.0);

        // Create container to remember patch j information
        int patch_j_done[patchInter[itype][jtype].nbP2];
        double costj_a_vect[patchInter[itype][jtype].nbP2];
        double costj_b_vect[patchInter[itype][jtype].nbP2];
        int patch_j_interact[patchInter[itype][jtype].nbP2];
	double nj_vect_0[patchInter[itype][jtype].nbP2];
	double nj_vect_1[patchInter[itype][jtype].nbP2];
	double nj_vect_2[patchInter[itype][jtype].nbP2];
	double aj_vect_0[patchInter[itype][jtype].nbP2]; 
	double aj_vect_1[patchInter[itype][jtype].nbP2];
	double aj_vect_2[patchInter[itype][jtype].nbP2];
	double bj_vect_0[patchInter[itype][jtype].nbP2];
	double bj_vect_1[patchInter[itype][jtype].nbP2];
	double bj_vect_2[patchInter[itype][jtype].nbP2];
	double rij_sca_bj_vect[patchInter[itype][jtype].nbP2];
	double rij_a_vect_0[patchInter[itype][jtype].nbP2];
	double rij_a_vect_1[patchInter[itype][jtype].nbP2];
	double rij_a_vect_2[patchInter[itype][jtype].nbP2];
	double rij_a_norm_vect[patchInter[itype][jtype].nbP2];
	double rij_sca_aj_vect[patchInter[itype][jtype].nbP2];
	double rij_b_vect_0[patchInter[itype][jtype].nbP2];
	double rij_b_vect_1[patchInter[itype][jtype].nbP2];
	double rij_b_vect_2[patchInter[itype][jtype].nbP2];
	double rij_b_norm_vect[patchInter[itype][jtype].nbP2];

	// NR: Loop over all patches
	for(int k(0) ; k < patchInter[itype][jtype].nbP1 ; k++)
	{
      // Calculate particle i patch orientation in space
      q_n_body_to_space(qi, pGeometry[patchInter[itype][jtype].typeP1[k]].nP, ni);
      q_n_body_to_space(qi, pGeometry[patchInter[itype][jtype].typeP1[k]].aP, ai);
      q_n_body_to_space(qi, pGeometry[patchInter[itype][jtype].typeP1[k]].bP, bi);
      MathExtra::norm3(ni);
      MathExtra::norm3(ai);
      MathExtra::norm3(bi);

      // Check if particle j is positioned such that it can interact with this patch on particle i
      rji_sca_bi = -rij[0]*bi[0] - rij[1]*bi[1] - rij[2]*bi[2];
      rji_a[0] = -rij[0] - rji_sca_bi*bi[0];
      rji_a[1] = -rij[1] - rji_sca_bi*bi[1];
      rji_a[2] = -rij[2] - rji_sca_bi*bi[2];
      rji_a_norm = sqrt(rji_a[0]*rji_a[0] + rji_a[1]*rji_a[1] + rji_a[2]*rji_a[2]);
      costi_a = MathExtra::dot3(ni,rji_a) / rji_a_norm;
      if(costi_a > 1){costi_a = 1.0;}if(costi_a < -1){costi_a = -1.0;}
     


      rji_sca_ai = -rij[0]*ai[0] - rij[1]*ai[1] - rij[2]*ai[2];
      rji_b[0] = -rij[0] - rji_sca_ai*ai[0];
      rji_b[1] = -rij[1] - rji_sca_ai*ai[1];
      rji_b[2] = -rij[2] - rji_sca_ai*ai[2];
      rji_b_norm = sqrt(rji_b[0]*rji_b[0] + rji_b[1]*rji_b[1] + rji_b[2]*rji_b[2]);
      costi_b = MathExtra::dot3(ni,rji_b) / rji_b_norm;
      if(costi_b > 1){costi_b = 1.0;}if(costi_b < -1){costi_b = -1.0;}

      if(costi_a >= pGeometry[patchInter[itype][jtype].typeP1[k]].cosThetaMax_a && costi_b >= pGeometry[patchInter[itype][jtype].typeP1[k]].cosThetaMax_b)
      {
          for(int l(0) ; l < patchInter[itype][jtype].nbP2 ; l++)
          {
            int idxPM(pGeometry[patchInter[itype][jtype].typeP1[k]].typeP*nbPatchType + pGeometry[patchInter[itype][jtype].typeP2[l]].typeP);
            if (rsq < pMorse[idxPM].cutoff*pMorse[idxPM].cutoff)
            {

              bool p_interact(false);

              //if(patch_j_done[l] != 666)
              if(true)
	      {
                  // Calculate particle j patch orientation in space
                  q_n_body_to_space(qj, pGeometry[patchInter[itype][jtype].typeP2[l]].nP, nj);
                  q_n_body_to_space(qj, pGeometry[patchInter[itype][jtype].typeP2[l]].aP, aj);
                  q_n_body_to_space(qj, pGeometry[patchInter[itype][jtype].typeP2[l]].bP, bj);
                  MathExtra::norm3(nj);
                  MathExtra::norm3(aj);
                  MathExtra::norm3(bj);

                  // Check if particle i is positioned such that it can interact with this patch on particle j
                  rij_sca_bj = rij[0]*bj[0] + rij[1]*bj[1] + rij[2]*bj[2];
                  rij_a[0] = rij[0] - rij_sca_bj*bj[0];
                  rij_a[1] = rij[1] - rij_sca_bj*bj[1];
                  rij_a[2] = rij[2] - rij_sca_bj*bj[2];
                  rij_a_norm = sqrt(rij_a[0]*rij_a[0] + rij_a[1]*rij_a[1] + rij_a[2]*rij_a[2]);
                  costj_a = MathExtra::dot3(nj,rij_a) / rij_a_norm;
		  if(costj_a > 1){costj_a = 1.0;}if(costj_a < -1){costj_a = -1.0;}

                  rij_sca_aj = rij[0]*aj[0] + rij[1]*aj[1] + rij[2]*aj[2];
                  rij_b[0] = rij[0] - rij_sca_aj*aj[0];
                  rij_b[1] = rij[1] - rij_sca_aj*aj[1];
                  rij_b[2] = rij[2] - rij_sca_aj*aj[2];
                  rij_b_norm = sqrt(rij_b[0]*rij_b[0] + rij_b[1]*rij_b[1] + rij_b[2]*rij_b[2]);
                  costj_b = MathExtra::dot3(nj,rij_b) / rij_b_norm;
		  if(costj_b > 1){costj_b = 1.0;}if(costj_b < -1){costj_b = -1.0;}
		 
                  costj_a_vect[l] = costj_a;
                  costj_b_vect[l] = costj_b;
		  
		  nj_vect_0[l] = nj[0]; nj_vect_1[l] = nj[1]; nj_vect_2[l] = nj[2];
		  aj_vect_0[l] = aj[0]; aj_vect_1[l] = aj[1]; aj_vect_2[l] = aj[2];
		  bj_vect_0[l] = bj[0]; bj_vect_1[l] = bj[1]; bj_vect_2[l] = bj[2];
		  rij_sca_bj_vect[l] = rij_sca_bj;
		  rij_a_vect_0[l] = rij_a[0]; rij_a_vect_1[l] = rij_a[1]; rij_a_vect_2[l] = rij_a[2];
		  rij_a_norm_vect[l] = rij_a_norm;
		  rij_sca_aj_vect[l] = rij_sca_aj;
		  rij_b_vect_0[l] = rij_b[0]; rij_b_vect_1[l] = rij_b[1]; rij_b_vect_2[l] = rij_b[2];
		  rij_b_norm_vect[l] = rij_b_norm;
                  
                  p_interact = (costj_a >= pGeometry[patchInter[itype][jtype].typeP2[l]].cosThetaMax_a) && (costj_b >= pGeometry[patchInter[itype][jtype].typeP2[l]].cosThetaMax_b);
                  patch_j_interact[l] = p_interact;
                  patch_j_done[l] = 666;
              }
              else
              {
                  p_interact = patch_j_interact[l];
              }

              if(p_interact)
              {
                    // Compute Morse and Morse derivation
                    double rA(exp( -pMorse[idxPM].b * ( rnorm - pMorse[idxPM].Re ) ));
                    expTerm  = 1.0 - rA ;
                    morseFactor = pMorse[idxPM].De * ( expTerm*expTerm - 1.0 );
                    morseDerivation = 2 * pMorse[idxPM].De * pMorse[idxPM].b * rA * (1.0 - rA);

                    ti_a = acos(costi_a) * 180.0 / pi_cst;
                    cosFactori_a = cos(pi_cst * ti_a / (2.0 * pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax_a));
                    ti_b = acos(costi_b) * 180.0 / pi_cst;
                    cosFactori_b = cos(pi_cst * ti_b / (2.0 * pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax_b));
                    ni_sca_rij = MathExtra::dot3(ni,rij);
                    ni_sca_rji_a = MathExtra::dot3(ni,rji_a);
                    ni_sca_rji_b = MathExtra::dot3(ni,rji_b);

                    costj_a = costj_a_vect[l];
                    costj_b = costj_b_vect[l];
		    nj[0] = nj_vect_0[l]; nj[1] = nj_vect_1[l]; nj[2] = nj_vect_2[l];
		    aj[0] = aj_vect_0[l]; aj[1] = aj_vect_1[l]; aj[2] = aj_vect_2[l];
		    bj[0] = bj_vect_0[l]; bj[1] = bj_vect_1[l]; bj[2] = bj_vect_2[l];
		    
		    rij_sca_bj = rij_sca_bj_vect[l];
		    rij_a[0] = rij_a_vect_0[l]; rij_a[1] = rij_a_vect_1[l]; rij_a[2] = rij_a_vect_2[l];
		    rij_a_norm = rij_a_norm_vect[l];
		    rij_sca_aj = rij_sca_aj_vect[l];
		    rij_b[0] = rij_b_vect_0[l]; rij_b[1] = rij_b_vect_1[l]; rij_b[2] = rij_b_vect_2[l];
		    rij_b_norm = rij_b_norm_vect[l];
		    
                    tj_a = acos(costj_a) * 180.0 / pi_cst;
                    tj_b = acos(costj_b) * 180.0 / pi_cst;
                    cosFactorj_a = cos(pi_cst * tj_a / (2.0 * pGeometry[patchInter[itype][jtype].typeP2[l]].thetaMax_a));
                    cosFactorj_b = cos(pi_cst * tj_b / (2.0 * pGeometry[patchInter[itype][jtype].typeP2[l]].thetaMax_b));
                    nj_sca_rij = MathExtra::dot3(nj,rij);
                    nj_sca_rij_a = MathExtra::dot3(nj,rij_a);
                    nj_sca_rij_b = MathExtra::dot3(nj,rij_b);
		    

                    angleTerm = cosFactori_a * cosFactori_b * cosFactorj_a * cosFactorj_b;

                    // Bunch of very useful variables
                    if(costi_a*costi_a == 1.0)
                    {
                      dti_dcosti_a = 0.0;
                    }
                    else
                    {
                      dti_dcosti_a = - 1.0 / sqrt(1.0 - costi_a*costi_a);
                    }
                    if(costi_b*costi_b == 1.0)
                    {
                      dti_dcosti_b = 0.0;
                    }
                    else
                    {
                      dti_dcosti_b = - 1.0 / sqrt(1.0 - costi_b*costi_b);
                    }

                    if(costj_a*costj_a == 1.0)
                    {
                      dtj_dcostj_a = 0.0;
                    }
                    else
                    {
                      dtj_dcostj_a = - 1.0 / sqrt(1.0 - costj_a*costj_a);
                    }
                    if(costj_b*costj_b == 1.0)
                    {
                      dtj_dcostj_b = 0.0;
                    }
                    else
                    {
                      dtj_dcostj_b = - 1.0 / sqrt(1.0 - costj_b*costj_b);
                    }
                    rji_a_norm_cube = rji_a_norm * rji_a_norm * rji_a_norm;
                    rji_b_norm_cube = rji_b_norm * rji_b_norm * rji_b_norm;
                    rij_a_norm_cube = rij_a_norm * rij_a_norm * rij_a_norm;
                    rij_b_norm_cube = rij_b_norm * rij_b_norm * rij_b_norm;

                    // Compute force
                    double prefacti_a = dti_dcosti_a * sin(pi_cst * ti_a / (2.0 * pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax_a)) * pi_cst / (2.0 * pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax_a * pi_cst / 180.0);
                    double prefacti_b = dti_dcosti_b * sin(pi_cst * ti_b / (2.0 * pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax_b)) * pi_cst / (2.0 * pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax_b * pi_cst / 180.0);
                    double prefactj_a = dtj_dcostj_a * sin(pi_cst * tj_a / (2.0 * pGeometry[patchInter[itype][jtype].typeP2[l]].thetaMax_a)) * pi_cst / (2.0 * pGeometry[patchInter[itype][jtype].typeP2[l]].thetaMax_a * pi_cst / 180.0);
                    double prefactj_b = dtj_dcostj_b * sin(pi_cst * tj_b / (2.0 * pGeometry[patchInter[itype][jtype].typeP2[l]].thetaMax_b)) * pi_cst / (2.0 * pGeometry[patchInter[itype][jtype].typeP2[l]].thetaMax_b * pi_cst / 180.0);
                    for(int m(0) ; m < 3 ; m++)
                    {
                        fi[m] += -(angleTerm*morseDerivation/rnorm) * rij[m] + morseFactor *
                        (
                         prefacti_a * ( -ni[m] / rji_a_norm - rji_a[m] * ni_sca_rij / rji_a_norm_cube ) * cosFactori_b * cosFactorj_a * cosFactorj_b +
                         prefacti_b * ( -ni[m] / rji_b_norm - rji_b[m] * ni_sca_rij / rji_b_norm_cube ) * cosFactori_a * cosFactorj_a * cosFactorj_b +
                         prefactj_a * ( nj[m] / rij_a_norm - rij_a[m] * nj_sca_rij / rij_a_norm_cube ) * cosFactori_a * cosFactori_b * cosFactorj_b +
                         prefactj_b * ( nj[m] / rij_b_norm - rij_b[m] * nj_sca_rij / rij_b_norm_cube ) * cosFactori_a * cosFactori_b * cosFactorj_a
                        );
                    }
                    /*
		    std::cout << "F = " << sqrt(fi[0]*fi[0] + fi[1]*fi[1] + fi[2]*fi[2]) << std::endl;
		    std::cout << "rji_a_norm = " << rji_a_norm << std::endl;
		    std::cout << "rji_b_norm = " << rji_b_norm << std::endl;
		    std::cout << "rij_a_norm = " << rij_a_norm << std::endl;
		    std::cout << "rij_b_norm = " << rij_b_norm << std::endl;
		    std::cout << "angleterm = " << angleTerm << std::endl;
		    std::cout << "morseDerivation = " << morseDerivation << std::endl;
		    std::cout << "rnorm = " << rnorm << std::endl;
		    std::cout << "cosFactor = " << cosFactori_a << "	" << cosFactori_b << "	" << cosFactorj_a << "	" << cosFactorj_b << std::endl;
		    std::cout << "see " << std::setprecision(12) << pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax_a << "	" << acos(costi_a)* 180.0 / pi_cst << "	" << std::setprecision(36) <<  costi_a << std::endl;
	            std::cout << "see " << std::setprecision(12) << pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax_a << "	" << acos(costi_b)* 180.0 / pi_cst << "	" <<  std::setprecision(36) << costi_b << std::endl;
		    */
		    /*
		    			std::cout << std::endl;
			std::cout << "rij_b_norm = " << rij_b_norm << std::endl;
			std::cout << rij[0] << "	" << rij[1] << "	" << rij[2] << std::endl;
			std::cout << nj[0] << "	" << nj[1] << "	" << nj[2] << std::endl;
			std::cout << aj[0] << "	" << aj[1] << "	" << aj[2] << std::endl;
			std::cout << bj[0] << "	" << bj[1] << "	" << bj[2] << std::endl;
			std::cout << "rij_sca_aj = " << rij_sca_aj << std::endl;
			std::cout << std::endl;
			*/
		    
                    // Compute torque on i
                    MathExtra::cross3(ni, rji_a, v1);
                    MathExtra::cross3(ni, rji_b, v2);
                    MathExtra::cross3(bi, rij, v3);
                    MathExtra::cross3(ai, rij, v4);
                    for(int m(0) ; m < 3 ; m++)
                    {
                        delti[m] +=  morseFactor *
                        (
                         v1[m] * prefacti_a * cosFactori_b * cosFactorj_a * cosFactorj_b / rji_a_norm +
                         v2[m] * prefacti_b * cosFactori_a * cosFactorj_a * cosFactorj_b / rji_b_norm +
                         v3[m] * prefacti_a * cosFactori_b * cosFactorj_a * cosFactorj_b * ni_sca_rji_a * (rji_sca_bi) / rji_a_norm_cube +
                         v4[m] * prefacti_b * cosFactori_a * cosFactorj_a * cosFactorj_b * ni_sca_rji_b * (rji_sca_ai) / rji_b_norm_cube
                        );
                    }

                    // Compute torque on j
                    MathExtra::cross3(nj, rij_a, v1);
                    MathExtra::cross3(nj, rij_b, v2);
                    MathExtra::cross3(bj, rij, v3);
                    MathExtra::cross3(aj, rij, v4);
                    for(int m(0) ; m < 3 ; m++)
                    {
                        deltj[m] +=  morseFactor *
                        (
                         v1[m] * prefactj_a * cosFactorj_b * cosFactori_a * cosFactori_b / rij_a_norm +
                         v2[m] * prefactj_b * cosFactorj_a * cosFactori_a * cosFactori_b / rij_b_norm +
                         v3[m] * prefactj_a * cosFactorj_b * cosFactori_a * cosFactori_b * nj_sca_rij_a * (-rij_sca_bj) / rij_a_norm_cube +
                         v4[m] * prefactj_b * cosFactorj_a * cosFactori_a * cosFactori_b * nj_sca_rij_b * (-rij_sca_aj) / rij_b_norm_cube
                        );
                    }

                     energyPatch += morseFactor * angleTerm;
		     //std::cout << " k l e " << k << "	" << l << "	" << morseFactor*angleTerm << std::endl;
		    //std::cout << ti_a << "	" << ti_b << "	" << tj_a << "	" << tj_b << std::endl;

              }
            }
          }

	}

	}

	/// NR, Scale force and torque for bond if needed ///
	fi[0] *= factor_lj;
	fi[1] *= factor_lj;
	fi[2] *= factor_lj;
	delti[0] *= factor_lj;
	delti[1] *= factor_lj;
	delti[2] *= factor_lj;

	/// NR, Add force and torque to atom i ///
        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
	torque[i][0] += delti[0];
        torque[i][1] += delti[1];
        torque[i][2] += delti[2];

	/// NR, Add force and torque to atom j if newton_pair ///
        if (newton_pair || j < nlocal) {
          f[j][0] -= fi[0];
          f[j][1] -= fi[1];
          f[j][2] -= fi[2];
	  deltj[0] *= factor_lj;
	  deltj[1] *= factor_lj;
	  deltj[2] *= factor_lj;
	  torque[j][0] += deltj[0];
	  torque[j][1] += deltj[1];
	  torque[j][2] += deltj[2];

	  /// NR, check local angular momentum conservation for debug purpose ///
	  
	  /*double torqueForceCst[3];
	  MathExtra::cross3(rij, fi, v1);
	  torqueForceCst[0] = delti[0] + deltj[0] + v1[0];
	  torqueForceCst[1] = delti[1] + deltj[1] + v1[1];
	  torqueForceCst[2] = delti[2] + deltj[2] + v1[2];
	  double tfCst(sqrt(torqueForceCst[0]*torqueForceCst[0] + torqueForceCst[1]*torqueForceCst[1] + torqueForceCst[2]*torqueForceCst[2]));
	  if( tfCst > 1e-10)
	  {
	    std::cout << "WARNING: FORCE-TORQUE EQUILIBRIUM WRONG !!!" << std::endl;
	    std::cout << tfCst << std::endl;
	  }*/
	  
	  
        }

	/// NR, Calculate energy ///
        if (eflag)
	{
          evdwl = factor_lj * energyPatch;
        }

        /// NR, Calculate stuff ///
        if (evflag)
	{
	  ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,fi[0],fi[1],fi[2],delx,dely,delz);
	}
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */



/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPatchyparticlerect::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  
  //NR
  memory->create(patchInter,n+1,n+1,"pair:patchinter");
  
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPatchyparticlerect::settings(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  
  
  //NR : global patch settings for the system
  nbPatchType = force->numeric(FLERR,arg[1]);
  if (narg < 3 + 4 * nbPatchType * nbPatchType) error->all(FLERR,"Illegal pair_style command");
  pMorse = new PatchMorseRect[nbPatchType*nbPatchType];
  int count(2);
  for(int i(0) ; i < nbPatchType*nbPatchType ; i++)
  {
    pMorse[i].cutoff = force->numeric(FLERR,arg[count]); count ++;
    pMorse[i].De = force->numeric(FLERR,arg[count]); count ++;
    pMorse[i].b = force->numeric(FLERR,arg[count]); count ++;
    pMorse[i].Re = force->numeric(FLERR,arg[count]); count ++;
  }
  nbPatchSystem = force->numeric(FLERR,arg[count]); count ++;
  if (narg < 3 + 4 * nbPatchType * nbPatchType + nbPatchSystem * 12) error->all(FLERR,"Illegal pair_style command");
  pGeometry = new PatchGeometryRect[nbPatchSystem];
  for(int i(0) ; i < nbPatchSystem ; i++)
  {
    pGeometry[i].typeP = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].nP[0] = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].nP[1] = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].nP[2] = force->numeric(FLERR,arg[count]); count ++;
    double scaleNP(sqrt(pGeometry[i].nP[0]*pGeometry[i].nP[0]+pGeometry[i].nP[1]*pGeometry[i].nP[1]+pGeometry[i].nP[2]*pGeometry[i].nP[2]));
    pGeometry[i].nP[0] = pGeometry[i].nP[0] / scaleNP;
    pGeometry[i].nP[1] = pGeometry[i].nP[1] / scaleNP;
    pGeometry[i].nP[2] = pGeometry[i].nP[2] / scaleNP;
    
    pGeometry[i].aP[0] = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].aP[1] = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].aP[2] = force->numeric(FLERR,arg[count]); count ++;
    double scaleAP(sqrt(pGeometry[i].aP[0]*pGeometry[i].aP[0]+pGeometry[i].aP[1]*pGeometry[i].aP[1]+pGeometry[i].aP[2]*pGeometry[i].aP[2]));
    pGeometry[i].aP[0] = pGeometry[i].aP[0] / scaleAP;
    pGeometry[i].aP[1] = pGeometry[i].aP[1] / scaleAP;
    pGeometry[i].aP[2] = pGeometry[i].aP[2] / scaleAP;
    
    pGeometry[i].bP[0] = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].bP[1] = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].bP[2] = force->numeric(FLERR,arg[count]); count ++;
    double scaleBP(sqrt(pGeometry[i].bP[0]*pGeometry[i].bP[0]+pGeometry[i].bP[1]*pGeometry[i].bP[1]+pGeometry[i].bP[2]*pGeometry[i].bP[2]));
    pGeometry[i].bP[0] = pGeometry[i].bP[0] / scaleBP;
    pGeometry[i].bP[1] = pGeometry[i].bP[1] / scaleBP;
    pGeometry[i].bP[2] = pGeometry[i].bP[2] / scaleBP;
    
    pGeometry[i].thetaMax_a = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].cosThetaMax_a = cos(pGeometry[i].thetaMax_a * pi_cst / 180.0);
    
    pGeometry[i].thetaMax_b = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].cosThetaMax_b = cos(pGeometry[i].thetaMax_b * pi_cst / 180.0);
  }
  
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPatchyparticlerect::coeff(int narg, char **arg)
{
  if (narg < 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  
  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int nbp1 = force->numeric(FLERR,arg[2]);
  int nbp2 = force->numeric(FLERR,arg[3]);
  
  int nbPatches(nbp1+nbp2);
  if ( narg < 4 + nbPatches || narg > 5 + nbPatches)
    error->all(FLERR,"Incorrect args for pair coefficients");
  
  if (!allocated) allocate();

  double cut_one = cut_global;
  if (narg == 5 + nbPatches ) cut_one = force->numeric(FLERR,arg[4+ nbPatches]);
  
  // Read patches
  PatchInteractionRect patchInterTemp;
  patchInterTemp.nbP1 = nbp1;
  patchInterTemp.nbP2 = nbp2;
  patchInterTemp.typeP1 = new int[nbp1];
  patchInterTemp.typeP2 = new int[nbp2];
  int count(0);
  for(int i(0) ; i < nbp1 ; i++)
  {
    patchInterTemp.typeP1[i] = force->numeric(FLERR,arg[4 + count]); count ++;
  }
  for(int i(0) ; i < nbp2 ; i++)
  {
    patchInterTemp.typeP2[i] = force->numeric(FLERR,arg[4 + count]); count ++;
  }
  

  count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      
      patchInter[i][j] = patchInterTemp;
      
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPatchyparticlerect::init_style()
{  
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR,"Pair gayberne requires atom style ellipsoid");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPatchyparticlerect::init_one(int i, int j)
{
  //std::cout << "HERE " << i << " " << j << std::endl;
  if (setflag[i][j] == 0) {
    patchInter[i][j].nbP1 = 0;
    patchInter[i][j].typeP1 = NULL;
    patchInter[i][j].nbP2 = 0;
    patchInter[i][j].typeP2 = NULL;
    cut[i][j] = 0.0;
  }

  if (offset_flag && (cut[i][j] > 0.0)) {
    /*double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));*/
    offset[i][j] = 0.0; ///// TO DO: THINK ABOUT THAT
  } else offset[i][j] = 0.0;

  // NR
  patchInter[j][i].nbP1 = patchInter[i][j].nbP2;
  patchInter[j][i].nbP2 = patchInter[i][j].nbP1;
  

  //if(patchInter[j][i].typeP1 == NULL)
  //{
    for(int k(0) ; k < patchInter[j][i].nbP1 ; k++)
    {
      if(k == 0) { patchInter[j][i].typeP1 = new int[patchInter[j][i].nbP1]; }
      patchInter[j][i].typeP1[k] = patchInter[i][j].typeP2[k];
    }
  //}
  //if(patchInter[j][i].typeP2 == NULL)
  //{
    for(int k(0) ; k < patchInter[j][i].nbP2 ; k++)
    {
      if(k == 0) { patchInter[j][i].typeP2 = new int[patchInter[j][i].nbP2]; }
      patchInter[j][i].typeP2[k] = patchInter[i][j].typeP1[k];
    }
  //}

  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPatchyparticlerect::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	
	// NR
        fwrite(&patchInter[i][j].nbP1,sizeof(int),1,fp);
	for(int k(0) ; k < patchInter[i][j].nbP1 ; k++)
	{
	  fwrite(&patchInter[i][j].typeP1[k],sizeof(int),1,fp);
	}
	fwrite(&patchInter[i][j].nbP2,sizeof(int),1,fp);
	for(int k(0) ; k < patchInter[i][j].nbP2 ; k++)
	{
	  fwrite(&patchInter[i][j].typeP2[k],sizeof(int),1,fp);
	}
        
	
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPatchyparticlerect::read_restart(FILE *fp)
{
  
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
	  
	// NR
        fread(&patchInter[i][j].nbP1,sizeof(int),1,fp);
	patchInter[i][j].typeP1 = new int[patchInter[i][j].nbP1];
	for(int k(0) ; k < patchInter[i][j].nbP1 ; k++)
	{
	  fread(&patchInter[i][j].typeP1[k],sizeof(int),1,fp);
	}
	fread(&patchInter[i][j].nbP2,sizeof(int),1,fp);
	patchInter[i][j].typeP2 = new int[patchInter[i][j].nbP2];
	for(int k(0) ; k < patchInter[i][j].nbP2 ; k++)
	{
	  fread(&patchInter[i][j].typeP2[k],sizeof(int),1,fp);
	}
	  
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        //NR
        MPI_Bcast(&patchInter[i][j].nbP1,1,MPI_INT,0,world);
	for(int k(0) ; k < patchInter[i][j].nbP1 ; k++)
	{
	  MPI_Bcast(&patchInter[i][j].typeP1[k],1,MPI_INT,0,world);
	}
        MPI_Bcast(&patchInter[i][j].nbP2,1,MPI_INT,0,world);
	for(int k(0) ; k < patchInter[i][j].nbP2 ; k++)
	{
	  MPI_Bcast(&patchInter[i][j].typeP2[k],1,MPI_INT,0,world);
	}
	
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPatchyparticlerect::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
  
  //NR
  fwrite(&nbPatchType,sizeof(int),1,fp);
  for(int k(0) ; k < nbPatchType*nbPatchType ; k++)
  {
    //fwrite(&pMorse[k],sizeof(PatchMorseRect),1,fp);
    fwrite(&pMorse[k].cutoff,sizeof(double),1,fp);
    fwrite(&pMorse[k].De,sizeof(double),1,fp);
    fwrite(&pMorse[k].b,sizeof(double),1,fp);
    fwrite(&pMorse[k].Re,sizeof(double),1,fp);
  }
  fwrite(&nbPatchSystem,sizeof(int),1,fp);
  for(int k(0) ; k < nbPatchSystem ; k++)
  {
    //fwrite(&pGeometry[k],sizeof(PatchGeometryRect),1,fp);
    fwrite(&pGeometry[k].typeP,sizeof(int),1,fp);
    fwrite(&pGeometry[k].nP[0],sizeof(double),1,fp);
    fwrite(&pGeometry[k].nP[1],sizeof(double),1,fp);
    fwrite(&pGeometry[k].nP[2],sizeof(double),1,fp);
    
    fwrite(&pGeometry[k].aP[0],sizeof(double),1,fp);
    fwrite(&pGeometry[k].aP[1],sizeof(double),1,fp);
    fwrite(&pGeometry[k].aP[2],sizeof(double),1,fp);
    
    fwrite(&pGeometry[k].bP[0],sizeof(double),1,fp);
    fwrite(&pGeometry[k].bP[1],sizeof(double),1,fp);
    fwrite(&pGeometry[k].bP[2],sizeof(double),1,fp);
    
    fwrite(&pGeometry[k].thetaMax_a,sizeof(double),1,fp);
    fwrite(&pGeometry[k].cosThetaMax_a,sizeof(double),1,fp);
    
    fwrite(&pGeometry[k].thetaMax_b,sizeof(double),1,fp);
    fwrite(&pGeometry[k].cosThetaMax_b,sizeof(double),1,fp);
  }  
  
  
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPatchyparticlerect::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
    
  //NR
  fread(&nbPatchType,sizeof(int),1,fp);
  pMorse = new PatchMorseRect[nbPatchType*nbPatchType];
  for(int k(0) ; k < nbPatchType ; k++)
  {
    //fread(&pMorse[k],sizeof(PatchMorseRect),1,fp);
    fread(&pMorse[k].cutoff,sizeof(double),1,fp);
    fread(&pMorse[k].De,sizeof(double),1,fp);
    fread(&pMorse[k].b,sizeof(double),1,fp);
    fread(&pMorse[k].Re,sizeof(double),1,fp);
  }
  fread(&nbPatchSystem,sizeof(int),1,fp);
  pGeometry = new PatchGeometryRect[nbPatchSystem];
  for(int k(0) ; k < nbPatchSystem ; k++)
  {
    //fread(&pGeometry[k],sizeof(PatchGeometryRect),1,fp);
    fread(&pGeometry[k].typeP,sizeof(int),1,fp);
    fread(&pGeometry[k].nP[0],sizeof(double),1,fp);
    fread(&pGeometry[k].nP[1],sizeof(double),1,fp);
    fread(&pGeometry[k].nP[2],sizeof(double),1,fp);
    
    fread(&pGeometry[k].aP[0],sizeof(double),1,fp);
    fread(&pGeometry[k].aP[1],sizeof(double),1,fp);
    fread(&pGeometry[k].aP[2],sizeof(double),1,fp);
    
    fread(&pGeometry[k].bP[0],sizeof(double),1,fp);
    fread(&pGeometry[k].bP[1],sizeof(double),1,fp);
    fread(&pGeometry[k].bP[2],sizeof(double),1,fp);
    
    fread(&pGeometry[k].thetaMax_a,sizeof(double),1,fp);
    fread(&pGeometry[k].cosThetaMax_a,sizeof(double),1,fp);
    
    fread(&pGeometry[k].thetaMax_b,sizeof(double),1,fp);
    fread(&pGeometry[k].cosThetaMax_b,sizeof(double),1,fp);
  } 
  
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
  
  //NR
  MPI_Bcast(&nbPatchType,1,MPI_INT,0,world);
  for(int k(0) ; k < nbPatchType ; k++)
  {
    //fread(&pMorse[k],sizeof(PatchMorseRect),1,fp);
    MPI_Bcast(&pMorse[k].cutoff,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pMorse[k].De,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pMorse[k].b,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pMorse[k].Re,1,MPI_DOUBLE,0,world);
    
  }
  MPI_Bcast(&nbPatchSystem,1,MPI_INT,0,world);
  for(int k(0) ; k < nbPatchSystem ; k++)
  {
    //fread(&pGeometry[k],sizeof(PatchGeometryRect),1,fp);
    MPI_Bcast(&pGeometry[k].typeP,1,MPI_INT,0,world);
    MPI_Bcast(&pGeometry[k].nP[0],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].nP[1],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].nP[2],1,MPI_DOUBLE,0,world);
    
    MPI_Bcast(&pGeometry[k].aP[0],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].aP[1],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].aP[2],1,MPI_DOUBLE,0,world);
    
    MPI_Bcast(&pGeometry[k].bP[0],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].bP[1],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].bP[2],1,MPI_DOUBLE,0,world);
    
    MPI_Bcast(&pGeometry[k].thetaMax_a,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].cosThetaMax_a,1,MPI_DOUBLE,0,world);
    
    MPI_Bcast(&pGeometry[k].thetaMax_b,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].cosThetaMax_b,1,MPI_DOUBLE,0,world);

  } 
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairPatchyparticlerect::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
  {
    fprintf(fp,"%d %d %d ",i, patchInter[i][i].nbP1, patchInter[i][i].nbP2);
    for(int k(0) ; k < patchInter[i][i].nbP1 ; k++)
    {
      fprintf(fp,"%d ",patchInter[i][i].typeP1[k]);
    }
    for(int k(0) ; k < patchInter[i][i].nbP2 ; k++)
    {
      fprintf(fp,"%d ",patchInter[i][i].typeP2[k]);
    }
    fprintf(fp,"\n");
    
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairPatchyparticlerect::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
  {
    for (int j = i; j <= atom->ntypes; j++)
    {
      fprintf(fp,"%d %d %d %d ",i,j, patchInter[i][j].nbP1, patchInter[i][j].nbP2);
      for(int k(0) ; k < patchInter[i][j].nbP1 ; k++)
      {
	fprintf(fp,"%d ",patchInter[i][j].typeP1[k]);
      }
      for(int k(0) ; k < patchInter[i][j].nbP2 ; k++)
      {
	fprintf(fp,"%d ",patchInter[i][j].typeP2[k]);
      }
      fprintf(fp,"%g \n",cut[i][j]);
    }
  }
}

/* ---------------------------------------------------------------------- */


void PairPatchyparticlerect::q_n_body_to_space(double *q, double *nbody, double *nspace)
{
  nspace[0] = (q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]) * nbody[0] + 2.0 * (q[1]*q[2] - q[0]*q[3]) * nbody[1] + 2.0 * (q[1]*q[3] + q[0]*q[2]) * nbody[2];
  nspace[1] = (2.0 * (q[1]*q[2] + q[0]*q[3])) * nbody[0] + (q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]) * nbody[1] + 2.0 * (q[2]*q[3] - q[0]*q[1])* nbody[2];
  nspace[2] = (2.0 * (q[1]*q[3] - q[0]*q[2])) * nbody[0] + (2.0 * (q[2]*q[3] + q[0]*q[1])) * nbody[1] + (q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]) * nbody[2];
}
