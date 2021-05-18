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

PairStyle(POL,PairPOL)

#else

#ifndef LMP_PAIR_POL_H
#define LMP_PAIR_POL_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPOL : public Pair {
 public:
  PairPOL(class LAMMPS *);
  virtual ~PairPOL();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_global;
  double *rad;
  double **cut,**fmag,**offset, **k_at, **f_rand, sigma2;
  int chain_type, Npol, **dotsign, *pollist, *NpolNeigh, *Jmin, **PolNeigh, boundlow, boundhigh;

  class RanMars *random;
  int seed;
  virtual void allocate();
};

}

#endif
#endif

/*
	Example calls/formatting:
	pair_coeff i j fmag chain_type k_one Npol boundlow boundhigh sig_rep f_rand seed cut
  pair_coeff i j fmag chain_type k_one Npol boundlow boundhigh sig_rep f_rand seed
 
  sig_rep: 

*/

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
