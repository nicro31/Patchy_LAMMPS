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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "bond_cgcell_harmonic.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

//NR
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"
#include <iostream>
#include <iomanip>
#define pi_cst 3.14159265359

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondCgcellHarmonic::BondCgcellHarmonic(LAMMPS *lmp) : Bond(lmp)
{
  reinitflag = 1;
}

/* ---------------------------------------------------------------------- */

BondCgcellHarmonic::~BondCgcellHarmonic()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
    memory->destroy(kA);
    memory->destroy(gammaA0);
    memory->destroy(kD);
    memory->destroy(nD0);
    memory->destroy(phiD0);
  }
}

/* ---------------------------------------------------------------------- */

void BondCgcellHarmonic::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,rk,rmu;
  double cosGammaa,deltaGammaa,kdga;
  double cosGammab,deltaGammab,kdgb;
  double eangle,tanglea,tangleb;
  double fi1[3],fi2[3];
  double fx, fy, fz, fmod, fmod_sqrtff;

  ebond = 0.0;
  eangle = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  
  //NR
  double **torque = atom->torque;
  double delta[3],deltb[3]; // torque increment;;
  double *qa,ax[3],ay[3],az[3];
  double *qb,bx[3],by[3],bz[3];
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  
  //NR, new for bond
  double bondx, bondy,bondz; 
  double fbondx, fbondy, fbondz;
  double tbondx1, tbondy1, tbondz1;
  double tbondx2, tbondy2, tbondz2;
  
  //NR, new for angle
  double anglecst; 
  double tanglex1, tangley1, tanglez1;
  

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    
    //NR
    qa=bonus[i1].quat;
    MathExtra::q_to_exyz(qa,ax,ay,az);
    qb=bonus[i2].quat;
    MathExtra::q_to_exyz(qb,bx,by,bz);
    //INVERT bx for second particle !!!
    bx[0] = -bx[0];
    bx[1] = -bx[1];
    bx[2] = -bx[2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    // iRef .... iDip
    
    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - r0[type];
    rk = k[type] * dr;
    rmu = r * 1.0;
    
    // NR, calculate bond potential
    bondx = -delx + r0[type]*(bx[0]-ax[0]);
    bondy = -dely + r0[type]*(bx[1]-ax[1]);
    bondz = -delz + r0[type]*(bx[2]-ax[2]);
    // Bond force on i1
    fbondx = 2*k[type]*bondx;
    fbondy = 2*k[type]*bondy;
    fbondz = 2*k[type]*bondz;
    // Torque on i1
    tbondx1 = ax[1]*fbondz*r0[type] - ax[2]*fbondy*r0[type];
    tbondy1 = ax[2]*fbondx*r0[type] - ax[0]*fbondz*r0[type];
    tbondz1 = ax[0]*fbondy*r0[type] - ax[1]*fbondx*r0[type];
    // Torque on i2
    tbondx2 = - bx[1]*fbondz*r0[type] + bx[2]*fbondy*r0[type];
    tbondy2 = - bx[2]*fbondx*r0[type] + bx[0]*fbondz*r0[type];
    tbondz2 = - bx[0]*fbondy*r0[type] + bx[1]*fbondx*r0[type];  
    
    // NR, calculate angle potential
    anglecst = ax[0]*bx[0] + ax[1]*bx[1] + ax[2]*bx[2] - cos(gammaA0[type]*pi_cst/180.0);
    // Torque on i1 (Note there is no force and torque on i2 is opposite to torque on i1)
    tanglex1 = -2*kA[type]*anglecst * (ax[1]*bx[2] - ax[2]*bx[1]);
    tangley1 = -2*kA[type]*anglecst * (ax[2]*bx[0] - ax[0]*bx[2]);
    tanglez1 = -2*kA[type]*anglecst * (ax[0]*bx[1] - ax[1]*bx[0]);
    

    ////////////////////////////////////////////////////////////////////
    
    
    // NR, calculate forces and torques from dihedral potential
    long double gammaD;
    long double phiD;
    double r_by[3];
    double r_ay[3];
    double rAB[3];
    long double n_r_by;
    long double n_r_ay;
    long double r_by__r_ay;
    long double r__by;
    long double r__ay;
    long double lambdaD;
    long double muD;
    rAB[0] = delx; 
    rAB[1] = dely;
    rAB[2] = delz;
    double skew_r_T[3][3];
    double skew_ay_T[3][3];
    double skew_by_T[3][3];
    double mu_tempA[3];
    double mu_tempB[3];
    double muA[3];
    double muB[3];
    double torqueAD[3];
    double torqueBD[3];
    double fAD[3];
    double fBD[3];
    
    skew_r_T[0][0] = 0; skew_r_T[0][1] = rAB[2]; skew_r_T[0][2] = -rAB[1];
    skew_r_T[1][0] = -rAB[2]; skew_r_T[1][1] = 0; skew_r_T[1][2] = rAB[0];
    skew_r_T[2][0] = rAB[1]; skew_r_T[2][1] = -rAB[0]; skew_r_T[2][2] = 0;
    
    skew_ay_T[0][0] = 0; skew_ay_T[0][1] = ay[2]; skew_ay_T[0][2] = -ay[1];
    skew_ay_T[1][0] = -ay[2]; skew_ay_T[1][1] = 0; skew_ay_T[1][2] = ay[0];
    skew_ay_T[2][0] = ay[1]; skew_ay_T[2][1] = -ay[0]; skew_ay_T[2][2] = 0;
    
    skew_by_T[0][0] = 0; skew_by_T[0][1] = by[2]; skew_by_T[0][2] = -by[1];
    skew_by_T[1][0] = -by[2]; skew_by_T[1][1] = 0; skew_by_T[1][2] = by[0];
    skew_by_T[2][0] = by[1]; skew_by_T[2][1] = -by[0]; skew_by_T[2][2] = 0;
    
    MathExtra::cross3(rAB,by,r_by);
    n_r_by = MathExtra::len3(r_by);
    MathExtra::cross3(rAB,ay,r_ay);
    n_r_ay = MathExtra::len3(r_ay);
    r_by__r_ay = MathExtra::dot3(r_by,r_ay);
    r__by = MathExtra::dot3(rAB,by);
    r__ay = MathExtra::dot3(rAB,ay);
    long double gSign(MathExtra::dot3(r_by,ay));
    gammaD = - ( (gSign > 0) ? 1 : -1 );
    lambdaD = r_by__r_ay / (n_r_by*n_r_ay);
    if(lambdaD > 1.0)
      lambdaD = 1.0;
    if(lambdaD < -1.0)
      lambdaD = -1.0;
    phiD = gammaD * acos(lambdaD);

    // Dihedral torque for A
    double v1[3];
    double v2[3];
    double v3[3];
    double v1_coeff(rsq/(n_r_ay*n_r_by));
    v1[0] = by[0]; v1[1] = by[1]; v1[2] = by[2];
    MathExtra::scale3(v1_coeff,v1);
    double v2_coeff(r__by/(n_r_ay*n_r_by));
    v2[0] = rAB[0]; v2[1] = rAB[1]; v2[2] = rAB[2];
    MathExtra::scale3(v2_coeff,v2);
    MathExtra::matvec(skew_r_T,r_ay,v3);
    double v3_coeff(r_by__r_ay/(n_r_ay*n_r_ay*n_r_ay*n_r_by));
    MathExtra::scale3(v3_coeff,v3);
    MathExtra::sub3(v1,v2,mu_tempA);
    MathExtra::sub3(mu_tempA,v3,muA);
    MathExtra::cross3(ay,muA,torqueAD);
    double torqueAD_coeff(-gammaD*nD0[type]*kD[type]*sin(nD0[type]*phiD-phiD0[type]*pi_cst/180.0) / sqrt(1-(lambdaD*lambdaD)));
    MathExtra::scale3(torqueAD_coeff,torqueAD);
    
    // Dihedral torque for B
    v1[0] = ay[0]; v1[1] = ay[1]; v1[2] = ay[2];
    MathExtra::scale3(v1_coeff,v1);
    v2_coeff = r__ay/(n_r_ay*n_r_by);
    v2[0] = rAB[0]; v2[1] = rAB[1]; v2[2] = rAB[2];
    MathExtra::scale3(v2_coeff,v2);
    MathExtra::matvec(skew_r_T,r_by,v3);
    v3_coeff = r_by__r_ay/(n_r_by*n_r_by*n_r_by*n_r_ay);
    MathExtra::scale3(v3_coeff,v3);
    MathExtra::sub3(v1,v2,mu_tempB);
    MathExtra::sub3(mu_tempB,v3,muB);
    MathExtra::cross3(by,muB,torqueBD);
    double torqueBD_coeff(torqueAD_coeff);
    MathExtra::scale3(torqueBD_coeff,torqueBD);  
    
    // Dihedral force for A and B
    double fD_coeff(torqueAD_coeff);
    double M[3];
    double N[3];
    v1[0] = rAB[0]; v1[1] = rAB[1]; v1[2] = rAB[2];
    double by__ay;
    by__ay = MathExtra::dot3(by,ay);
    v1_coeff = 2*by__ay;
    MathExtra::scale3(v1_coeff,v1);
    v2[0] = ay[0]; v2[1] = ay[1]; v2[2] = ay[2];
    v2_coeff = -r__by;
    MathExtra::scale3(v2_coeff, v2);
    v3[0] = by[0]; v3[1] = by[1]; v3[2] = by[2];
    v3_coeff = -r__ay;
    MathExtra::scale3(v3_coeff, v3);
    M[0] = v1[0] + v2[0] + v3[0];
    M[1] = v1[1] + v2[1] + v3[1];
    M[2] = v1[2] + v2[2] + v3[2];
    double M_coeff(1.0/(n_r_by*n_r_ay));
    MathExtra::scale3(M_coeff, M);
    double N1V[3];
    MathExtra::matvec(skew_ay_T,r_ay,N1V);
    double N2V[3];
    MathExtra::matvec(skew_by_T,r_by,N2V);
    double N1_coeff(1.0/(n_r_by*n_r_ay*n_r_ay*n_r_ay));
    double N2_coeff(1.0/(n_r_ay*n_r_by*n_r_by*n_r_by));
    MathExtra::scale3(N1_coeff,N1V);
    MathExtra::scale3(N2_coeff,N2V);
    MathExtra::add3(N1V,N2V,N);
    MathExtra::scale3(r_by__r_ay,N);
    double mu[3];
    MathExtra::add3(M,N,mu);
    MathExtra::scale3(fD_coeff, mu);
    fAD[0] = mu[0]; fAD[1] = mu[1]; fAD[2] = mu[2];
    fBD[0] = -mu[0]; fBD[1] = -mu[1]; fBD[2] = -mu[2];
    double edihedral(0.0);
    
    if(lambdaD == 1.0 || lambdaD == -1.0)
    {
      double torque_coeff_lambdaCase(kD[type]*nD0[type]*sin(phiD0[type]*pi_cst/180.0)/r);
      fAD[0] = 0;fAD[1] = 0;fAD[2] = 0;
      fBD[0] = 0;fBD[1] = 0;fBD[2] = 0;
      torqueAD[0] = -torque_coeff_lambdaCase*rAB[0];
      torqueAD[1] = -torque_coeff_lambdaCase*rAB[1];
      torqueAD[2] = -torque_coeff_lambdaCase*rAB[2];
      torqueBD[0] = torque_coeff_lambdaCase*rAB[0];
      torqueBD[1] = torque_coeff_lambdaCase*rAB[1];
      torqueBD[2] = torque_coeff_lambdaCase*rAB[2];
    }
    

    
    // force & energy

    if (r > 0.0) fbond = -2.0*rk/r;
    else fbond = 0.0;
    if (eflag)
    {
      ebond = k[type] * (bondx*bondx+bondy*bondy+bondz*bondz);
      eangle = kA[type] * anglecst * anglecst; 
      edihedral = kD[type] * (1.0 + cos(nD0[type]*phiD-phiD0[type]*pi_cst/180.0));
    }
    
    //double fangle(sqrt(fi1[0]*fi1[0] + fi1[1]*fi1[1] + fi1[2]*fi1[2]));
    
    double fdihedral(sqrt(fAD[0]*fAD[0] + fAD[1]*fAD[1] + fAD[2]*fAD[2]));
  



    // apply force and torque to each of 2 atoms
    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fbondx + fAD[0];
      f[i1][1] += fbondy + fAD[1];
      f[i1][2] += fbondz + fAD[2];

      // NR
      torque[i1][0] += tbondx1 + tanglex1 + torqueAD[0];
      torque[i1][1] += tbondy1 + tangley1 + torqueAD[1];
      torque[i1][2] += tbondz1 + tanglez1 + torqueAD[2];
      
      
    }
    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= fbondx - fBD[0];
      f[i2][1] -= fbondy - fBD[1];
      f[i2][2] -= fbondz - fBD[2];
      
      //NR
      torque[i2][0] += tbondx2 - tanglex1 + torqueBD[0];
      torque[i2][1] += tbondy2 - tangley1 + torqueBD[1];
      torque[i2][2] += tbondz2 - tanglez1 + torqueBD[2];
    }
    
    /// TO CORRECT ///
    //if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond + eangle + edihedral,fbond+fangle+fdihedral,delx,dely,delz);
    if (evflag) ev_tally_xyz(i1,i2,nlocal,newton_bond,ebond + eangle + edihedral,fbondx + fAD[0],fbondy + fAD[1],fbondz + fAD[2],delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondCgcellHarmonic::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(r0,n+1,"bond:r0");
  memory->create(kA,n+1,"bond:kA");
  memory->create(gammaA0,n+1,"bond:gammaA0");
  memory->create(kD,n+1,"bond:kD");
  memory->create(nD0,n+1,"bond:nD0");
  memory->create(phiD0,n+1,"bond:phiD0");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondCgcellHarmonic::coeff(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double r0_one = force->numeric(FLERR,arg[2]);
  double kA_one = force->numeric(FLERR,arg[3]);
  double gammaA0_one = force->numeric(FLERR,arg[4]);
  double kD_one = force->numeric(FLERR,arg[5]);
  double nD0_one = force->numeric(FLERR,arg[6]);
  double phiD0_one = force->numeric(FLERR,arg[7]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    kA[i] = kA_one;
    gammaA0[i] = gammaA0_one;
    kD[i] = kD_one;
    nD0[i] = nD0_one;
    phiD0[i] = phiD0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondCgcellHarmonic::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondCgcellHarmonic::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&kA[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&gammaA0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&kD[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&nD0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&phiD0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondCgcellHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&kA[1],sizeof(double),atom->nbondtypes,fp);
    fread(&gammaA0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&kD[1],sizeof(double),atom->nbondtypes,fp);
    fread(&nD0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&phiD0[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kA[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gammaA0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kD[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&nD0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&phiD0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondCgcellHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g\n",i,k[i],r0[i],kA[i],gammaA0[i],kD[i],nD0[i],phiD0[i]);
}

/* ---------------------------------------------------------------------- */

double BondCgcellHarmonic::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  int i1,i2,n;
  double delx,dely,delz,ebond,fbond;
  double r,dr,rk,rmu;
  double cosGammaa,deltaGammaa,kdga;
  double cosGammab,deltaGammab,kdgb;
  double eangle,tanglea,tangleb;
  double fi1[3],fi2[3];
  double fx, fy, fz, fmod, fmod_sqrtff;

  ebond = 0.0;
  eangle = 0.0;
  //ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  
  //NR
  double **torque = atom->torque;
  double delta[3],deltb[3]; // torque increment;;
  double *qa,ax[3],ay[3],az[3];
  double *qb,bx[3],by[3],bz[3];
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  
  //NR, new for bond
  double bondx, bondy,bondz; 
  double fbondx, fbondy, fbondz;
  double tbondx1, tbondy1, tbondz1;
  double tbondx2, tbondy2, tbondz2;
  
  //NR, new for angle
  double anglecst; 
  double tanglex1, tangley1, tanglez1;
  
  
    i1 = i;
    i2 = j;
    
    //NR
    qa=bonus[i1].quat;
    MathExtra::q_to_exyz(qa,ax,ay,az);
    qb=bonus[i2].quat;
    MathExtra::q_to_exyz(qb,bx,by,bz);
    //INVERT bx for second particle !!!
    bx[0] = -bx[0];
    bx[1] = -bx[1];
    bx[2] = -bx[2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    // iRef .... iDip
    
    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - r0[type];
    rk = k[type] * dr;
    
    //NR, calculate angle potential
    rmu = r * 1.0;
    
    // NR, calculate bond potential
    bondx = -delx + r0[type]*(bx[0]-ax[0]);
    bondy = -dely + r0[type]*(bx[1]-ax[1]);
    bondz = -delz + r0[type]*(bx[2]-ax[2]);
    // Bond force on i1
    fbondx = 2*k[type]*bondx;
    fbondy = 2*k[type]*bondy;
    fbondz = 2*k[type]*bondz;
    // Torque on i1
    tbondx1 = ax[1]*fbondz*r0[type] - ax[2]*fbondy*r0[type];
    tbondy1 = ax[2]*fbondx*r0[type] - ax[0]*fbondz*r0[type];
    tbondz1 = ax[0]*fbondy*r0[type] - ax[1]*fbondx*r0[type];
    // Torque on i2
    tbondx2 = - bx[1]*fbondz*r0[type] + bx[2]*fbondy*r0[type];
    tbondy2 = - bx[2]*fbondx*r0[type] + bx[0]*fbondz*r0[type];
    tbondz2 = - bx[0]*fbondy*r0[type] + bx[1]*fbondx*r0[type];  
    
    // NR, calculate angle potential
    anglecst = ax[0]*bx[0] + ax[1]*bx[1] + ax[2]*bx[2] - cos(gammaA0[type]*pi_cst/180.0);
    // Torque on i1 (Note there is no force and torque on i2 is opposite to torque on i1)
    tanglex1 = -2*kA[type]*anglecst * (ax[1]*bx[2] - ax[2]*bx[1]);
    tangley1 = -2*kA[type]*anglecst * (ax[2]*bx[0] - ax[0]*bx[2]);
    tanglez1 = -2*kA[type]*anglecst * (ax[0]*bx[1] - ax[1]*bx[0]);

    ////////////////////////////////////////////////////////////////////
    
    
    // NR, calculate forces and torques from dihedral potential
    long double gammaD;
    long double phiD;
    double r_by[3];
    double r_ay[3];
    double rAB[3];
    long double n_r_by;
    long double n_r_ay;
    long double r_by__r_ay;
    long double r__by;
    long double r__ay;
    long double lambdaD;
    long double muD;
    rAB[0] = delx; 
    rAB[1] = dely;
    rAB[2] = delz;
    double skew_r_T[3][3];
    double skew_ay_T[3][3];
    double skew_by_T[3][3];
    double mu_tempA[3];
    double mu_tempB[3];
    double muA[3];
    double muB[3];
    double torqueAD[3];
    double torqueBD[3];
    double fAD[3];
    double fBD[3];
    
    skew_r_T[0][0] = 0; skew_r_T[0][1] = rAB[2]; skew_r_T[0][2] = -rAB[1];
    skew_r_T[1][0] = -rAB[2]; skew_r_T[1][1] = 0; skew_r_T[1][2] = rAB[0];
    skew_r_T[2][0] = rAB[1]; skew_r_T[2][1] = -rAB[0]; skew_r_T[2][2] = 0;
    
    skew_ay_T[0][0] = 0; skew_ay_T[0][1] = ay[2]; skew_ay_T[0][2] = -ay[1];
    skew_ay_T[1][0] = -ay[2]; skew_ay_T[1][1] = 0; skew_ay_T[1][2] = ay[0];
    skew_ay_T[2][0] = ay[1]; skew_ay_T[2][1] = -ay[0]; skew_ay_T[2][2] = 0;
    
    skew_by_T[0][0] = 0; skew_by_T[0][1] = by[2]; skew_by_T[0][2] = -by[1];
    skew_by_T[1][0] = -by[2]; skew_by_T[1][1] = 0; skew_by_T[1][2] = by[0];
    skew_by_T[2][0] = by[1]; skew_by_T[2][1] = -by[0]; skew_by_T[2][2] = 0;
    
    MathExtra::cross3(rAB,by,r_by);
    n_r_by = MathExtra::len3(r_by);
    MathExtra::cross3(rAB,ay,r_ay);
    n_r_ay = MathExtra::len3(r_ay);
    r_by__r_ay = MathExtra::dot3(r_by,r_ay);
    r__by = MathExtra::dot3(rAB,by);
    r__ay = MathExtra::dot3(rAB,ay);
    long double gSign(MathExtra::dot3(r_by,ay));
    gammaD = - ( (gSign > 0) ? 1 : -1 );
    lambdaD = r_by__r_ay / (n_r_by*n_r_ay);
    if(lambdaD > 1.0)
      lambdaD = 1.0;
    if(lambdaD < -1.0)
      lambdaD = -1.0;
    phiD = gammaD * acos(lambdaD);

    // Dihedral torque for A
    double v1[3];
    double v2[3];
    double v3[3];
    double v1_coeff(rsq/(n_r_ay*n_r_by));
    v1[0] = by[0]; v1[1] = by[1]; v1[2] = by[2];
    MathExtra::scale3(v1_coeff,v1);
    double v2_coeff(r__by/(n_r_ay*n_r_by));
    v2[0] = rAB[0]; v2[1] = rAB[1]; v2[2] = rAB[2];
    MathExtra::scale3(v2_coeff,v2);
    MathExtra::matvec(skew_r_T,r_ay,v3);
    double v3_coeff(r_by__r_ay/(n_r_ay*n_r_ay*n_r_ay*n_r_by));
    MathExtra::scale3(v3_coeff,v3);
    MathExtra::sub3(v1,v2,mu_tempA);
    MathExtra::sub3(mu_tempA,v3,muA);
    MathExtra::cross3(ay,muA,torqueAD);
    double torqueAD_coeff(-gammaD*nD0[type]*kD[type]*sin(nD0[type]*phiD-phiD0[type]*pi_cst/180.0) / sqrt(1-(lambdaD*lambdaD)));
    MathExtra::scale3(torqueAD_coeff,torqueAD);
    
    // Dihedral torque for B
    v1[0] = ay[0]; v1[1] = ay[1]; v1[2] = ay[2];
    MathExtra::scale3(v1_coeff,v1);
    v2_coeff = r__ay/(n_r_ay*n_r_by);
    v2[0] = rAB[0]; v2[1] = rAB[1]; v2[2] = rAB[2];
    MathExtra::scale3(v2_coeff,v2);
    MathExtra::matvec(skew_r_T,r_by,v3);
    v3_coeff = r_by__r_ay/(n_r_by*n_r_by*n_r_by*n_r_ay);
    MathExtra::scale3(v3_coeff,v3);
    MathExtra::sub3(v1,v2,mu_tempB);
    MathExtra::sub3(mu_tempB,v3,muB);
    MathExtra::cross3(by,muB,torqueBD);
    double torqueBD_coeff(torqueAD_coeff);
    MathExtra::scale3(torqueBD_coeff,torqueBD);  
    
    // Dihedral force for A and B
    double fD_coeff(torqueAD_coeff);
    double M[3];
    double N[3];
    v1[0] = rAB[0]; v1[1] = rAB[1]; v1[2] = rAB[2];
    double by__ay;
    by__ay = MathExtra::dot3(by,ay);
    v1_coeff = 2*by__ay;
    MathExtra::scale3(v1_coeff,v1);
    v2[0] = ay[0]; v2[1] = ay[1]; v2[2] = ay[2];
    v2_coeff = -r__by;
    MathExtra::scale3(v2_coeff, v2);
    v3[0] = by[0]; v3[1] = by[1]; v3[2] = by[2];
    v3_coeff = -r__ay;
    MathExtra::scale3(v3_coeff, v3);
    M[0] = v1[0] + v2[0] + v3[0];
    M[1] = v1[1] + v2[1] + v3[1];
    M[2] = v1[2] + v2[2] + v3[2];
    double M_coeff(1.0/(n_r_by*n_r_ay));
    MathExtra::scale3(M_coeff, M);
    double N1V[3];
    MathExtra::matvec(skew_ay_T,r_ay,N1V);
    double N2V[3];
    MathExtra::matvec(skew_by_T,r_by,N2V);
    double N1_coeff(1.0/(n_r_by*n_r_ay*n_r_ay*n_r_ay));
    double N2_coeff(1.0/(n_r_ay*n_r_by*n_r_by*n_r_by));
    MathExtra::scale3(N1_coeff,N1V);
    MathExtra::scale3(N2_coeff,N2V);
    MathExtra::add3(N1V,N2V,N);
    MathExtra::scale3(r_by__r_ay,N);
    double mu[3];
    MathExtra::add3(M,N,mu);
    MathExtra::scale3(fD_coeff, mu);
    fAD[0] = mu[0]; fAD[1] = mu[1]; fAD[2] = mu[2];
    fBD[0] = -mu[0]; fBD[1] = -mu[1]; fBD[2] = -mu[2];
    double edihedral(0.0);
    
    if(lambdaD == 1.0 || lambdaD == -1.0)
    {
      double torque_coeff_lambdaCase(kD[type]*nD0[type]*sin(phiD0[type]*pi_cst/180.0)/r);
      fAD[0] = 0;fAD[1] = 0;fAD[2] = 0;
      fBD[0] = 0;fBD[1] = 0;fBD[2] = 0;
      torqueAD[0] = -torque_coeff_lambdaCase*rAB[0];
      torqueAD[1] = -torque_coeff_lambdaCase*rAB[1];
      torqueAD[2] = -torque_coeff_lambdaCase*rAB[2];
      torqueBD[0] = torque_coeff_lambdaCase*rAB[0];
      torqueBD[1] = torque_coeff_lambdaCase*rAB[1];
      torqueBD[2] = torque_coeff_lambdaCase*rAB[2];
    }
    

    
    // force & energy

    if (r > 0.0) fbond = sqrt(fbondx*fbondx + fbondy*fbondy+fbondz*fbondz);
    else fbond = 0.0;

    ebond = k[type] * (bondx*bondx+bondy*bondy+bondz*bondz);
    eangle = kA[type] * anglecst * anglecst; 
    edihedral = kD[type] * (1.0 + cos(nD0[type]*phiD-phiD0[type]*pi_cst/180.0));
  
    double fangle(0);
    double fdihedral(sqrt(fAD[0]*fAD[0] + fAD[1]*fAD[1] + fAD[2]*fAD[2]));
  
    
    
  fforce = 0;
  if (r > 0.0) fforce = fbond + fangle/r + fdihedral/r;
  return ebond + eangle + edihedral;
}

/* ----------------------------------------------------------------------
    Return ptr to internal members upon request.
------------------------------------------------------------------------ */
void *BondCgcellHarmonic::extract( char *str, int &dim )
{
  dim = 1;
  if (strcmp(str,"kappa")==0) return (void*) k;
  if (strcmp(str,"r0")==0) return (void*) r0;
  if (strcmp(str,"kA")==0) return (void*) kA;
  if (strcmp(str,"gammaA0")==0) return (void*) gammaA0;
  if (strcmp(str,"kD")==0) return (void*) kD;
  if (strcmp(str,"nD0")==0) return (void*) nD0;
  if (strcmp(str,"phiD0")==0) return (void*) phiD0;
  return NULL;
}


