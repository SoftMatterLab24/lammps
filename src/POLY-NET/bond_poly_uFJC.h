/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(poly/uFJC,BondPolyuFJC);
// clang-format on
#else

#ifndef LMP_BOND_POLY_uFJC_H
#define LMP_BOND_POLY_uFJC_H

#include "bond_poly.h"

namespace LAMMPS_NS {

class BondPolyuFJC : public BondPoly {
 public:
  BondPolyuFJC(class LAMMPS *);
  ~BondPolyuFJC() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void settings(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, double, int, int, double &) override;
  //int pack_forward_comm(int, int *, double *, int, int *) override;
  //void unpack_forward_comm(int, int, double *) override;
  //int pack_reverse_comm(int, int, double *) override;
  //void unpack_reverse_comm(int, int *, double *) override;

 protected:
  double *k0, *kappa, *fcrit, *gamma, *lamc;
  int smooth_flag, normalize_flag, stretch_flag, flag_distribution;

  int nmax;
  char *id_fix_property_bond;

  struct Table {
   int ninput;
   int *iatomfile, *jatomfile;
   double r0;
   double *Nfile, *bfile;
   double *N, *b;
  };

  // Pointers for random numbers
  int seed;
  class RanMars *random;

  // parameters for pdf styles
  double lo, hi, mu, sigma, lambda, l_min, l_max, alpha, beta;
  double b_dist;
  char *dist_type = nullptr;

  int tabstyle, tablength, ntables, *tabindex;
  Table *tables;

  void allocate();
  void store_data();
  double store_bond(int, int, int);

  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, char *, char *);
  void bcast_table(Table *);
  
  void param_extract(Table *, char *);
  void update_table(int);
  void bond_lookup(int, int, int, double &, double &);
  double sample_cdf(long int, int);
  double bond_uniform(long int, int);
  double barrier(double, double, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif