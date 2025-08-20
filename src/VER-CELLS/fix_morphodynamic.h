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
FixStyle(morphodynamic,FixMorphoDynamic);
// clang-format on
#else

//PRB{These Statements ensure that this class is only inlcuded once in the project}
#ifndef LMP_FIX_MORPHODYNAMIC_H
#define LMP_FIX_MORPHODYNAMIC_H

#include "fix.h"
#include <vector>

using namespace std;

namespace LAMMPS_NS {

class FixMorphoDynamic : public Fix {
 public:
  FixMorphoDynamic(class LAMMPS *, int, char **);
  ~FixMorphoDynamic() override;
  void post_constructor();
  int setmask() override;
  void init() override;
  void setup(int) override;
  //void init_list(int, class NeighList *) override;
  // void min_setup(int) override;

  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  void grow_arrays(int) override;
  void set_arrays(int) override;
  double memory_usage() override;

 protected:
  int me, nprocs;
  int nmax;

 private:
  FILE *fp1;    //write neighs list to txt file
  FILE *fp2;    //write vertices list to txt file
  FILE *fp3;    //write num of T1 transitions to txt file
  std::string file1, file2, file3;
  //class NeighList *list;

  //Default values (change if required)
  int per_atom_out = 10;    //Number of columns in per atom data output
  int neighs_MAX = 30;      //Maximum number of neighbours for each atom

  //Voronoi Force Data members
  class Compute *vcompute;     // To read data from compute voronoi (only for initialization)
  char *id_compute_voronoi;    // A Pointer variable

  /*CONSTANT STYLE VARIABLES*/

  //simulation params
  double rcut, ksoft;    //soft repsulsion betweeen cell center
  double c1, c2;         //stopping cell centers from crossing bonds
  int num_MC;            //Number of MC iterations

  /*EQUAL STYLE VARIABLES - evaluated at every time step*/
  char *nu0_str, *nevery_output_str;
  int nu0_style, nevery_output_style;
  int nu0_var, nevery_output_var;
  int nevery_output;     //Output data frequency
  double nu0;           //friction coefficient, number of MC iterations

  //fix berendsen
  char *berflagx_str, *berflagy_str, *berdamp_str;    // later read as input
  int  berflagx_style, berflagy_style, berdamp_style;
  int  berflagx_var, berflagy_var, berdamp_var;
  int  berflagx, berflagy;
  double berdamp;

  //Tissue properties
  char *kpstr, *p0str, *KTstr, *neveryMCstr, *Jsstr, *Jnstr, *fastr;
  int kpstyle, p0style, KTstyle, neveryMCstyle, Jsstyle, Jnstyle, fastyle;
  int kpvar, p0var, KTvar, neveryMCvar, Jsvar, Jnvar, favar;

  double kp, p0, KT, Js, Jn, fa;
  int neveryMC;

  //flag variables
  int flag_init;

  // Declare per atom arrays
  // Note that they do persist accorss time steps
  // But, we might not be looking at the same atom at the next time step
  // because an atom might have moved to a different processor
  // Ideally use them if on single processor
  // OR if you don't need to access per atom data from the previous time step.

  double **cell_shape;    //Note that these are dynamic arrays
  double **vertices;

  // Declare local array
  int **bond_life;    // Again this is a dynamic array
  int num_bonds;

  class RanMars *wgn;
  int seed_wgn;

  // For invoking fix property/atom
  char *new_fix_id;    //To store atom property: neighs_list
  int index;

  // Common (Dyn Tri force) Data members

  char *idregion;
  class Region *region;    // A pointer variable (but we don't delete it!)

  int ilevel_respa;

  int commflag;    // for communicating data of ghost atoms

  //Declare functions

  // For cyclically arranging a set of points
  void arrange_cyclic(tagint *, int, int);

  // For storing cell info: area, peri and energy
  void get_cell_data1(double *, double *, tagint *, int, int);
  void get_cell_data2(int, double *, int, double **, tagint *);
  void get_trans_state_geom(double *, double[][2], int);
  double get_trans_state_energy(int, tagint *, int, double[][4], double **, tagint[4]);

  void sortverticesCCW(double[][2], int, double, double);

  // For obtaining common neighs for a pair of atoms
  int get_comm_neigh(tagint *, tagint *, tagint *, int, int);

  // See if atoms are already bonded or not

  int unbonded(tagint *, int, tagint);

  //adding and deleting neighs
  //void add_neigh(int, tagint, tagint *, int);
  void remove_neigh(int, tagint, tagint *, int, int);

  void print_neighs_list(tagint *, int, int);

  bool isConcave(double *, double *, double *, double *);
  double crossProduct(double *, double *, double *);
  void getCP(double *, double *, double *);
  void normalize(double *);

  //void arrange_cyclic_new(tagint *, int, int, double **);
  bool is_cyclic_perm(tagint *, tagint *, int);
  bool is_center_outside(double[], int, int, double, double);

  // to check self intersecting polygons
  bool onSegment(double, double, double, double, double, double);
  int orientation(double, double, double, double, double, double);
  bool doIntersect(double, double, double, double, double, double, double, double);
  bool isSelfIntersecting(double[], int);
  void arrange_cyclic_moved(tagint *, int, int, double[][2]);
  void get_cell_config_moved(double *, tagint *, int, int, double[][2]);
  int nearest_image(int, int, double[][2]);

  // Implementing Fix press/berendsen
  std::vector<Fix *> rfix;    // indices of rigid fixes
  double press_current[2];
  double press_start[2];
  double press_stop[2] = {0.0};    // This models traction free boundary conditions
  double press_target[2];
  double dilation_MDN[2];

  void resize_box();
};

}    // namespace LAMMPS_NS

#endif
#endif