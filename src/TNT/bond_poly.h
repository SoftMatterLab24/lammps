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
BondStyle(poly, BondPOLY);
// clang-format on
#else

#ifndef LMP_BOND_POLY_H
#define LMP_BOND_POLY_H

#include "bond.h"

#include <vector>

namespace LAMMPS_NS {

class Fix;

class BondPOLY : public Bond {
    public:
     BondPOLY(class LAMMPS *);
     ~BondPOLY() override;
     void compute(int, int) override;
     void init_style() override;
     void settings(int, char **) override;
     void coeff(int, char **) override;
     double equilibrium_distance(int) override;
     void write_restart(FILE *) override;
     void read_restart(FILE *) override;
     void write_restart_settings(FILE *) override;
     void read_restart_settings(FILE *) override;
     int pack_forward_comm(int, int *, double *, int, int *) override;
     void unpack_forward_comm(int, int, double *) override;
     int pack_reverse_comm(int, int, double *) override;
     void unpack_reverse_comm(int, int *, double *) override;
     double single(int, double, int, int, double &) override;

    protected:

     double *r0, *k, *fc;
     int store_local_freq, nhistory, update_flag, hybrid_flag;

     int  nmax;
     char *id_fix_property_bond;

     struct Table {
      int ninput, fpflag;
      double fplo, fphi, r0;
      double lo, hi;
      double *nfile, *bfile;
      double *n, *b;
     };

     int tabstyle, tablength;
     
  
     int ntables;
     Table *tables;
     int *tabindex;
    
     void allocate();
     void store_data();
     double store_bond(int, int, int);

     void null_table(Table *);
     void free_table(Table *);
     void read_table(Table *, char *, char *);
     void bcast_table(Table *);
  
     void param_extract(Table *, char *);
     void bond_lookup(int, int, double &, double &);

     std::vector<int> leftover_iarg;

     char *id_fix_dummy, *id_fix_dummy2;
     char *id_fix_update, *id_fix_bond_history;
     char *id_fix_store_local, *id_fix_prop_atom;
     class FixStoreLocal *fix_store_local;
     class FixBondHistory *fix_bond_history;
     class FixUpdateSpecialBonds *fix_update_special_bonds;

     void process_broken(int, int);
     typedef void (BondPOLY::*FnPtrPack)(int, int, int);
     FnPtrPack *pack_choice;    // ptrs to pack functions
     double *output_data;

     int prop_atom_flag, nvalues, overlay_flag, break_flag, ignore_special_flag;
     int index_x_ref, index_y_ref, index_z_ref;

     int n_histories;
     std::vector<Fix *> histories;

     void pack_id1(int, int, int);
     void pack_id2(int, int, int);
     void pack_time(int, int, int);
     void pack_x(int, int, int);
     void pack_y(int, int, int);
     void pack_z(int, int, int);
     void pack_x_ref(int, int, int);
     void pack_y_ref(int, int, int);
     void pack_z_ref(int, int, int);
     
   };
   
   }    // namespace LAMMPS_NS
   
   #endif
   #endif