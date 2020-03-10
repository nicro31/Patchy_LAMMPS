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
#include "pair_patchyparticle.h"
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

PairPatchyparticle::PairPatchyparticle(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairPatchyparticle::~PairPatchyparticle()
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

void PairPatchyparticle::compute(int eflag, int vflag)
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
  double costi;
  double ti;
  double costj;
  double tj;
  double rji[3];
  double rij[3];
  double rnorm;
  double expTerm;
  double morseFactor;
  double morseDerivation;
  double angleTerm;
  double cosFactori;
  double cosFactorj;
  double dti_dcosti;
  double dtj_dcostj;
  double dti_drij[3];
  double dtj_drij[3];
  double v1[3];
  double v2[3];
  double v3[3];
  double v4[3];
  double v1Coeff;
  double v2Coeff;
  double v3Coeff;
  double v4Coeff;
  double fi[3];
  double fi_temp[3];

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
      
      if (rsq < cutsq[itype][jtype]) 
      {
	
	double energyPatch(0.0);
	
	// NR: Loop over all patches
	for(int k(0) ; k < patchInter[itype][jtype].nbP1 ; k++)
	{
	  for(int l(0) ; l < patchInter[itype][jtype].nbP2 ; l++)
	  {
	    int idxPM(pGeometry[patchInter[itype][jtype].typeP1[k]].typeP*nbPatchType + pGeometry[patchInter[itype][jtype].typeP2[l]].typeP);
	 	    
	    if (rsq < pMorse[idxPM].cutoff*pMorse[idxPM].cutoff)
	    {
	      // Calculate patch orientation in space
	      qi = bonus[i].quat;
	      qj = bonus[j].quat;
	      q_n_body_to_space(qi, pGeometry[patchInter[itype][jtype].typeP1[k]].nP, ni);
	      MathExtra::norm3(ni);
	      q_n_body_to_space(qj, pGeometry[patchInter[itype][jtype].typeP2[l]].nP, nj);
	      MathExtra::norm3(nj);
	      
	      // Check if patches interact
	      costi = - MathExtra::dot3(ni,rij) / rnorm;
	      costj = MathExtra::dot3(nj,rij) / rnorm;
	      
	      if(costi >= pGeometry[patchInter[itype][jtype].typeP1[k]].cosThetaMax && costj >= pGeometry[patchInter[itype][jtype].typeP2[l]].cosThetaMax)
	      {	
		
		// Compute force
		double rA(exp( -pMorse[idxPM].b * ( rnorm - pMorse[idxPM].Re ) ));
		expTerm  = 1.0 - rA ;
		morseFactor = pMorse[idxPM].De * ( expTerm*expTerm - 1.0 );
		ti = acos(costi) * 180.0 / pi_cst;
		tj = acos(costj) * 180.0 / pi_cst;
		cosFactori = cos(pi_cst * ti / (2.0 * pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax));
		cosFactorj = cos(pi_cst * tj / (2.0 * pGeometry[patchInter[itype][jtype].typeP2[l]].thetaMax));
		angleTerm = cosFactori * cosFactorj;
		
		if(costi*costi == 1.0)
		{
		  dti_dcosti = 0.0;
		}
		else
		{
		  dti_dcosti = - 1.0 / sqrt(1.0 - costi*costi);
		}
		if(costj*costj == 1.0)
		{
		  dtj_dcostj = 0.0;
		}
		else
		{
		  dtj_dcostj = - 1.0 / sqrt(1.0 - costj*costj);
		}
		
		MathExtra::copy3(ni, v1);
		MathExtra::scale3(-1.0/rnorm,v1);
		MathExtra::copy3(rij,v2);
		v2Coeff = (MathExtra::dot3(ni,rij)/(rnorm*rnorm*rnorm));
		MathExtra::scale3(v2Coeff,v2);
		MathExtra::add3(v1, v2, dti_drij);
		
		MathExtra::copy3(nj, v1);
		MathExtra::scale3(1.0/rnorm,v1);
		MathExtra::copy3(rij,v2);
		v2Coeff = (-MathExtra::dot3(nj,rij)/(rnorm*rnorm*rnorm));
		MathExtra::scale3(v2Coeff,v2);
		MathExtra::add3(v1, v2, dtj_drij);
		
		morseDerivation = 2 * pMorse[idxPM].De * pMorse[idxPM].b * rA * (1.0 - rA);
		MathExtra::copy3(rij,v4);
		MathExtra::scale3(-angleTerm*morseDerivation/rnorm,v4);
		
		MathExtra::copy3(dti_drij, v1);
		v1Coeff = cosFactorj * dti_dcosti * sin(pi_cst * ti / (2.0 * pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax)) * pi_cst / (2.0 * pGeometry[patchInter[itype][jtype].typeP1[k]].thetaMax * pi_cst / 180.0); 
		MathExtra::copy3(dtj_drij, v2);
		v2Coeff = cosFactori * dtj_dcostj * sin(pi_cst * tj / (2.0 * pGeometry[patchInter[itype][jtype].typeP2[l]].thetaMax)) * pi_cst / (2.0 * pGeometry[patchInter[itype][jtype].typeP2[l]].thetaMax * pi_cst / 180.0);
		MathExtra::scale3(v1Coeff, v1);
		MathExtra::scale3(v2Coeff, v2);
		MathExtra::add3(v1,v2,v3);
		MathExtra::scale3(morseFactor, v3);
		MathExtra::add3(v4, v3, fi_temp);
		
		fi[0] = fi[0] + fi_temp[0];
		fi[1] = fi[1] + fi_temp[1];
		fi[2] = fi[2] + fi_temp[2];
		
		// Compute torques
		// On i
		 MathExtra::cross3(rij, ni, v1);
		 MathExtra::scale3(morseFactor * v1Coeff/rnorm, v1);
		 delti[0] = delti[0] + v1[0];
		 delti[1] = delti[1] + v1[1];
		 delti[2] = delti[2] + v1[2];
		 
		 // On j
		 MathExtra::cross3(rij, nj, v2);
		 MathExtra::scale3(-morseFactor * v2Coeff/rnorm, v2);
		 deltj[0] = deltj[0] + v2[0];
		 deltj[1] = deltj[1] + v2[1];
		 deltj[2] = deltj[2] + v2[2];
		 
		 energyPatch += morseFactor * angleTerm;
		
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
	  /*
	  double torqueForceCst[3];
	  MathExtra::cross3(rij, fi, v1);
	  torqueForceCst[0] = delti[0] + deltj[0] + v1[0];
	  torqueForceCst[1] = delti[1] + deltj[1] + v1[1];
	  torqueForceCst[2] = delti[2] + deltj[2] + v1[2];
	  double tfCst(sqrt(torqueForceCst[0]*torqueForceCst[0] + torqueForceCst[1]*torqueForceCst[1] + torqueForceCst[2]*torqueForceCst[2]));
	  if( tfCst > 1e-10)
	  {
	    std::cout << "WARNING: FORCE-TORQUE EQUILIBRIUM WRONG !!!" << std::endl;
	    std::cout << tfCst << std::endl;
	  }
	  */
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

void PairPatchyparticle::allocate()
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

void PairPatchyparticle::settings(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  
  
  //NR : global patch settings for the system
  nbPatchType = force->numeric(FLERR,arg[1]);
  if (narg < 3 + 4 * nbPatchType * nbPatchType) error->all(FLERR,"Illegal pair_style command");
  pMorse = new PatchMorse[nbPatchType*nbPatchType];
  int count(2);
  for(int i(0) ; i < nbPatchType*nbPatchType ; i++)
  {
    pMorse[i].cutoff = force->numeric(FLERR,arg[count]); count ++;
    pMorse[i].De = force->numeric(FLERR,arg[count]); count ++;
    pMorse[i].b = force->numeric(FLERR,arg[count]); count ++;
    pMorse[i].Re = force->numeric(FLERR,arg[count]); count ++;
  }
  nbPatchSystem = force->numeric(FLERR,arg[count]); count ++;
  if (narg < 3 + 4 * nbPatchType * nbPatchType + nbPatchSystem * 5) error->all(FLERR,"Illegal pair_style command");
  pGeometry = new PatchGeometry[nbPatchSystem];
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
    pGeometry[i].thetaMax = force->numeric(FLERR,arg[count]); count ++;
    pGeometry[i].cosThetaMax = cos(pGeometry[i].thetaMax * pi_cst / 180.0);
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

void PairPatchyparticle::coeff(int narg, char **arg)
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
  PatchInteraction patchInterTemp;
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

void PairPatchyparticle::init_style()
{  
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR,"Pair gayberne requires atom style ellipsoid");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPatchyparticle::init_one(int i, int j)
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

void PairPatchyparticle::write_restart(FILE *fp)
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

void PairPatchyparticle::read_restart(FILE *fp)
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

void PairPatchyparticle::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
  
  //NR
  fwrite(&nbPatchType,sizeof(int),1,fp);
  for(int k(0) ; k < nbPatchType*nbPatchType ; k++)
  {
    //fwrite(&pMorse[k],sizeof(PatchMorse),1,fp);
    fwrite(&pMorse[k].cutoff,sizeof(double),1,fp);
    fwrite(&pMorse[k].De,sizeof(double),1,fp);
    fwrite(&pMorse[k].b,sizeof(double),1,fp);
    fwrite(&pMorse[k].Re,sizeof(double),1,fp);
  }
  fwrite(&nbPatchSystem,sizeof(int),1,fp);
  for(int k(0) ; k < nbPatchSystem ; k++)
  {
    //fwrite(&pGeometry[k],sizeof(PatchGeometry),1,fp);
    fwrite(&pGeometry[k].typeP,sizeof(int),1,fp);
    fwrite(&pGeometry[k].nP[0],sizeof(double),1,fp);
    fwrite(&pGeometry[k].nP[1],sizeof(double),1,fp);
    fwrite(&pGeometry[k].nP[2],sizeof(double),1,fp);
    fwrite(&pGeometry[k].thetaMax,sizeof(double),1,fp);
    fwrite(&pGeometry[k].cosThetaMax,sizeof(double),1,fp);
  }  
  
  
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPatchyparticle::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
    
  //NR
  fread(&nbPatchType,sizeof(int),1,fp);
  pMorse = new PatchMorse[nbPatchType*nbPatchType];
  for(int k(0) ; k < nbPatchType ; k++)
  {
    //fread(&pMorse[k],sizeof(PatchMorse),1,fp);
    fread(&pMorse[k].cutoff,sizeof(double),1,fp);
    fread(&pMorse[k].De,sizeof(double),1,fp);
    fread(&pMorse[k].b,sizeof(double),1,fp);
    fread(&pMorse[k].Re,sizeof(double),1,fp);
  }
  fread(&nbPatchSystem,sizeof(int),1,fp);
  pGeometry = new PatchGeometry[nbPatchSystem];
  for(int k(0) ; k < nbPatchSystem ; k++)
  {
    //fread(&pGeometry[k],sizeof(PatchGeometry),1,fp);
    fread(&pGeometry[k].typeP,sizeof(int),1,fp);
    fread(&pGeometry[k].nP[0],sizeof(double),1,fp);
    fread(&pGeometry[k].nP[1],sizeof(double),1,fp);
    fread(&pGeometry[k].nP[2],sizeof(double),1,fp);
    fread(&pGeometry[k].thetaMax,sizeof(double),1,fp);
    fread(&pGeometry[k].cosThetaMax,sizeof(double),1,fp);
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
    //fread(&pMorse[k],sizeof(PatchMorse),1,fp);
    MPI_Bcast(&pMorse[k].cutoff,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pMorse[k].De,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pMorse[k].b,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pMorse[k].Re,1,MPI_DOUBLE,0,world);
    
  }
  MPI_Bcast(&nbPatchSystem,1,MPI_INT,0,world);
  for(int k(0) ; k < nbPatchSystem ; k++)
  {
    //fread(&pGeometry[k],sizeof(PatchGeometry),1,fp);
    MPI_Bcast(&pGeometry[k].typeP,1,MPI_INT,0,world);
    MPI_Bcast(&pGeometry[k].nP[0],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].nP[1],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].nP[2],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].thetaMax,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&pGeometry[k].cosThetaMax,1,MPI_DOUBLE,0,world);

  } 
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairPatchyparticle::write_data(FILE *fp)
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

void PairPatchyparticle::write_data_all(FILE *fp)
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


void PairPatchyparticle::q_n_body_to_space(double *q, double *nbody, double *nspace)
{
  nspace[0] = (q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]) * nbody[0] + 2.0 * (q[1]*q[2] - q[0]*q[3]) * nbody[1] + 2.0 * (q[1]*q[3] + q[0]*q[2]) * nbody[2];
  nspace[1] = (2.0 * (q[1]*q[2] + q[0]*q[3])) * nbody[0] + (q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]) * nbody[1] + 2.0 * (q[2]*q[3] - q[0]*q[1])* nbody[2];
  nspace[2] = (2.0 * (q[1]*q[3] - q[0]*q[2])) * nbody[0] + (2.0 * (q[2]*q[3] + q[0]*q[1])) * nbody[1] + (q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]) * nbody[2];
}
