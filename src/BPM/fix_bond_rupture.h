/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(bond/rupture,FixBondRupture);
// clang-format on
#else

#ifndef LMP_FIX_BOND_RUPTURE_H
#define LMP_FIX_BOND_RUPTURE_H

#include "fix.h"
namespace LAMMPS_NS {

class FixBondRupture : public Fix {
 public:
  FixBondRupture(class LAMMPS *, int, char **);
  ~FixBondRupture() override;
  int setmask() override;
  void init() override;
  void post_integrate() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_reverse_comm_size(int, int) override;

 protected:
  int me, nprocs, ndata, update_flag;

  int commflag;
  int nmax, maxbond;
  double probability;
  double **bprob;

  // Partner tracking for cross-processor bonds
  tagint *break_partner = nullptr;

  // pointer to shared bond history fix (if created)
  class FixBondHistory *fix_bond_history = nullptr;
  char *id_fix_bond_history_rupture = nullptr;

  int n_histories;
  std::vector<Fix *> histories;

  // internal bondstore bookkeeping (delegated to FixBondHistory)
  int updated_bond_flag = 0;
  int *setflag = nullptr;

  // Pointers for random numbers
  class RanMars *random;

  // parameters for pdf styles
  double *lo, *hi, *mu, *sigma, *lambda, *alpha, *beta;
  int *use_dist;
  double *dist_data;
  char *dist_type;

  // Default arguments
  int btype, seed;
  double rcrit, rcritsq, rcritsq_g, lamc, p_fraction, zeta, kappa, omega, dt_eq, k0, f0, ks0, kc0, fs0, fc0;

  // Flags for styles
  int flag_dist, flag_stretch, flag_fraction, flag_slip, flag_slip_catch, flag_rate;

  // Flags for keywords
  int flag_table, flag_distribution, flag_crit, flag_prob, flag_tilt;

  // Internal methods/functions
  void store_data();
  double store_bond(int, int, int);
  double sample_cdf(long int, int);
  double bond_uniform(long int, int);
  double barrier(double, double, double, double);

  void process_broken(int, int);
  void update_topology();

  // Create an array to store bonds broken this timestep
  // and since the last neighbor list build
  std::vector<std::pair<tagint, tagint>> new_broken_pairs;
  std::vector<std::pair<tagint, tagint>> new_created_pairs;
  
  // Structs and methods for reading/processing table files
  // Define struct for bond table
  struct Table {
   int ninput;
   int *iatomfile, *jatomfile;
   double *datafile;
   double *data;
  };

  int tabstyle, tablength, ntables, *tabindex;
  Table *tables;

  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, char *, char *);
  void bcast_table(Table *);
  void param_extract(Table *, char *);

};

}    // namespace LAMMPS_NS

#endif
#endif

