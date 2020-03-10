/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(patchyparticle,PairPatchyparticle)

#else

#ifndef LMP_PAIR_PATCHYPARTICLE_H
#define LMP_PAIR_PATCHYPARTICLE_H

#include "pair.h"

struct PatchInteraction
{
  // Patches on particle 1
  int nbP1;
  int *typeP1;
  
  // Patches on particle 2
  int nbP2;
  int *typeP2;
};

struct PatchMorse 
{
  // Cutoff for interaction
  double cutoff;
  
  // Morse potential parameters
  double De;
  double b;
  double Re;
};

struct PatchGeometry
{
  // Patch type
  int typeP;
  
  // Angular parameters
  double nP[3];
  double thetaMax;  
  double cosThetaMax;
};

namespace LAMMPS_NS {

class PairPatchyparticle : public Pair {
 public:
  PairPatchyparticle(class LAMMPS *);
  virtual ~PairPatchyparticle();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);


 protected:
  double cut_global;
  double **cut;
  double **offset;
  
  //NR : global patch settings for the system
  int nbPatchType;
  PatchMorse *pMorse;
  int nbPatchSystem;
  PatchGeometry *pGeometry;
  
  //NR : patch settings per particle pair
  PatchInteraction **patchInter;

  virtual void allocate();
  
  // Add function to calculate transformation of vector from body-frame to space-frame
  void q_n_body_to_space(double *q, double *nbody, double *nspace);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/
