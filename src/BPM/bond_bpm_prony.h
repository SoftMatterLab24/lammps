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
BondStyle(bpm/prony,BondBPMProny);
// clang-format on
#else

#ifndef LMP_BOND_BPM_PRONY_H
#define LMP_BOND_BPM_PRONY_H

#include "bond_bpm.h"

namespace LAMMPS_NS {

class BondBPMProny : public BondBPM {
 public:
  BondBPMProny(class LAMMPS *);
  ~BondBPMProny() override;
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
  double *k0, *ecrit, *gamma, *lamc, *eplastic, *aT;
  int smooth_flag, normalize_flag, nonlinear_flag, plastic_flag, temperature_flag;

  int index_vol, index_vol0, nmax;
  char *id_fix_property_bond;
  double *len_current, **H;
  double *aT_temp;
  double dt_temp;

  struct Table {
   int ninput, nninput, fpflag;
   double fplo, fphi, r0;
   double lo, hi;
   double *kfile, *etafile, *expfile, *Nfile, *bfile;
   double *k, *eta, *expj, *N, *b;
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
  void nonlinear_param_extract(Table *, char *);
  void update_table(int);
};

}    // namespace LAMMPS_NS

#endif
#endif