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

namespace LAMMPS_NS {

class BondPOLY : public Bond {
    public:
     BondPOLY(class LAMMPS *_lmp) : Bond(_lmp) {}
     ~BondPOLY() override;
     void compute(int, int) override;
     void settings(int, char **) override;
     void coeff(int, char **) override;
     double equilibrium_distance(int) override;
     void write_restart(FILE *) override;
     void read_restart(FILE *) override;
     void write_restart_settings(FILE *) override;
     void read_restart_settings(FILE *) override;
     //void write_data(FILE *) override;
     double single(int, double, int, int, double &) override;
     
     //void *extract(const char *, int &) override;
   
    protected:
     int tabstyle, tablength;
     double *r0;
  
     struct Table {
      int ninput, fpflag;
      double fplo, fphi, r0;
      double lo, hi;
      double *nfile, *bfile;
      double *n, *b;
     };
  
     int ntables;
     Table *tables;
     int *tabindex;
  
     void allocate();
     void null_table(Table *);
     void free_table(Table *);
     void read_table(Table *, char *, char *);
     void bcast_table(Table *);
  
     void param_extract(Table *, char *);
  
     void bond_lookup(int, int, double &, double &);
   
     
   };
   
   }    // namespace LAMMPS_NS
   
   #endif
   #endif