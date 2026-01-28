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
BondStyle(bpm/gkv,BondBPMGKV);
// clang-format on
#else

#ifndef LMP_BOND_BPM_GKV_H
#define LMP_BOND_BPM_GKV_H

#include "bond_bpm.h"

namespace LAMMPS_NS {

class BondBPMGKV : public BondBPM {
 public:
  BondBPMGKV(class LAMMPS *);
  ~BondBPMGKV() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void settings(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, double, int, int, double &) override;
  void *extract(const char *, int &) override;
  //int pack_forward_comm(int, int *, double *, int, int *) override;
  //void unpack_forward_comm(int, int, double *) override;
  //int pack_reverse_comm(int, int, double *) override;
  //void unpack_reverse_comm(int, int *, double *) override;

 protected:
  double *Ks, *Kj, *rcrit, *gamma, *zeta, *aT, *lamc;
  int smooth_flag, normalize_flag, temperature_flag, stretch_flag;

  int nmax;
  char *id_fix_property_bond;
  double *aT_temp;
  double dt_temp;

  struct Table {
   int ninput;
   double r0;
   int *iatomfile, *jatomfile, *Nfile, *N;
   double *bfile, *b;
  };

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
  void bond_lookup(int, int, int, int &, double &);
  void direct_solve(double, int, int, double, double &);
  void iter_solve(double, int, int, double, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif