/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#include "fix_morphodynamic.h"

#include "arg_info.h"
#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
// #include "cell.hh"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h" // Include the header file for the Error class
#include "fix_deform.h"
#include "group.h"
#include "input.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "random_mars.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <algorithm>
#include <bits/stdc++.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <random>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <voro++.hh>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;
using namespace std::chrono;
using namespace voro;

enum { NONE, CONSTANT, EQUAL, ATOM };

/* ---------------------------------------------------------------------- */

FixMorphoDynamic::FixMorphoDynamic(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_compute_voronoi(nullptr), cell_shape(nullptr), vertices(nullptr),
    bond_life(nullptr), wgn(nullptr), idregion(nullptr), region(nullptr), kpstr(nullptr), 
    p0str(nullptr), KTstr(nullptr), neveryMCstr(nullptr), Jsstr(nullptr), Jnstr(nullptr),
    fastr(nullptr), berflagx_str(nullptr), berflagy_str(nullptr),
    nevery_output_str(nullptr), nu0_str(nullptr), berdamp_str(nullptr)
{
  if (narg < 21) error->all(FLERR, "Illegal fix MorphoDynamic command: not sufficient args");

  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);

  // initialize Marsaglia RNG with processor-unique seed (for adding noise in active motion)
  // seed_wgn = 1;
  // wgn = new RanMars(lmp, seed_wgn + comm->me);
  // srand(10);

  dynamic_group_allow = 1;
  energy_peratom_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  thermo_energy = thermo_virial = 1;
  evflag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  /*Read in the simulation input paramters*/

  // for initial setup
  id_compute_voronoi = utils::strdup(arg[3]);

  rcut = utils::numeric(FLERR, arg[4], false, lmp);
  ksoft = utils::numeric(FLERR, arg[5], false, lmp);
  c1 = utils::numeric(FLERR, arg[6], false, lmp);
  c2 = utils::numeric(FLERR, arg[7], false, lmp);
  num_MC = utils::numeric(FLERR, arg[8], false, lmp);

  if (utils::strmatch(arg[9], "^v_")) {
    nevery_output_str = utils::strdup(arg[9] + 2);
  } else {
    nevery_output = utils::numeric(FLERR, arg[9], false, lmp);
    nevery_output_style = CONSTANT;
  }

  if (utils::strmatch(arg[10], "^v_")) {
    nu0_str = utils::strdup(arg[10] + 2);
  } else {
    nu0 = utils::numeric(FLERR, arg[10], false, lmp);
    nu0_style = CONSTANT;
  }
  
  //Berendsen flag variables
  if (utils::strmatch(arg[11], "^v_")) {
    berdamp_str = utils::strdup(arg[11] + 2);
  } else {
    berdamp = utils::numeric(FLERR, arg[11], false, lmp);
    berdamp_style = CONSTANT;
  }
  if (utils::strmatch(arg[12], "^v_")) {
    berflagx_str = utils::strdup(arg[12] + 2);
  } else {
    berflagx = utils::numeric(FLERR, arg[12], false, lmp);
    berflagx_style = CONSTANT;
  }
  if (utils::strmatch(arg[13], "^v_")) {
    berflagy_str = utils::strdup(arg[13] + 2);
  } else {
    berflagy = utils::numeric(FLERR, arg[13], false, lmp);
    berflagy_style = CONSTANT;
  }

  //Biophysical params

  if (utils::strmatch(arg[14], "^v_")) {
    kpstr = utils::strdup(arg[14] + 2);
  } else {
    kp = utils::numeric(FLERR, arg[14], false, lmp);
    kpstyle = CONSTANT;
  }
  if (utils::strmatch(arg[15], "^v_")) {
    p0str = utils::strdup(arg[15] + 2);
  } else {
    p0 = utils::numeric(FLERR, arg[15], false, lmp);
    p0style = CONSTANT;
  }
  if (utils::strmatch(arg[16], "^v_")) {
    KTstr = utils::strdup(arg[16] + 2);
  } else {
    KT = utils::numeric(FLERR, arg[16], false, lmp);
    KTstyle = CONSTANT;
  }
  if (utils::strmatch(arg[17], "^v_")) {
    neveryMCstr = utils::strdup(arg[17] + 2);
  } else {
    neveryMC = utils::numeric(FLERR, arg[17], false, lmp);
    neveryMCstyle = CONSTANT;
  }
  if (utils::strmatch(arg[18], "^v_")) {
    Jsstr = utils::strdup(arg[18] + 2);
  } else {
    Js = utils::numeric(FLERR, arg[18], false, lmp);
    Jsstyle = CONSTANT;
  }
  if (utils::strmatch(arg[19], "^v_")) {
    Jnstr = utils::strdup(arg[19] + 2);
  } else {
    Jn = utils::numeric(FLERR, arg[19], false, lmp);
    Jnstyle = CONSTANT;
  }
  if (utils::strmatch(arg[20], "^v_")) {
    fastr = utils::strdup(arg[20] + 2);
  } else {
    fa = utils::numeric(FLERR, arg[20], false, lmp);
    fastyle = CONSTANT;
  }

  /*This fix takes in input a per-atom array
  produced by compute voronoi*/

  vcompute = modify->get_compute_by_id(id_compute_voronoi);
  if (!vcompute)
    error->all(FLERR, "Could not find compute ID {} for voronoi compute", id_compute_voronoi);

  //parse values for optional arguments
  nevery = 1;    // Using default value for now

  if (narg > 21) {
    idregion = utils::strdup(arg[21]);
    region = domain->get_region_by_id(idregion);
  }

  // Initialize nmax and virial pointer
  nmax = atom->nmax;
  cell_shape = nullptr;
  vertices = nullptr;
  bond_life = nullptr;

  // Specify attributes for dumping connectivity (neighs_array)
  // This fix generates a per-atom array with specified columns as output,
  // containing information for owned atoms (nlocal on each processor) (accessed from dump file)

  peratom_flag = 1;
  size_peratom_cols = per_atom_out + 2;    // 2 extra columns for x and y
  peratom_freq = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class
  if (peratom_flag) {
    FixMorphoDynamic::grow_arrays(atom->nmax);
    atom->add_callback(Atom::GROW);
  }

  //Define file names for outputting
  file1 = "neighs_list.txt";
  file2 = "vertices_list.txt";
  file3 = "num_T1.txt";        //change to num_T1 if intereted in number of T1 transitions

  flag_init = 0;
}

/* ---------------------------------------------------------------------- */

FixMorphoDynamic::~FixMorphoDynamic()
{ 
  delete[] id_compute_voronoi;
  delete[] nevery_output_str;
  delete[] nu0_str;
  delete[] berdamp_str;
  delete[] berflagx_str;
  delete[] berflagy_str;
  delete[] kpstr;
  delete[] p0str;
  delete[] KTstr;
  delete[] neveryMCstr;
  delete[] Jsstr;
  delete[] Jnstr;
  delete[] fastr;

  delete[] idregion;
  delete wgn;

  memory->destroy(cell_shape);
  memory->destroy(vertices);
  memory->destroy(bond_life);

  // unregister callbacks to this fix from atom class
  if (peratom_flag) { atom->delete_callback(id, Atom::GROW); }

  if (new_fix_id && modify->nfix) modify->delete_fix(new_fix_id);
  delete[] new_fix_id;
}

/* ---------------------------------------------------------------------- */
// returntype classname :: functidentifier(args)

int FixMorphoDynamic::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMorphoDynamic::post_constructor()
{
  // Create call to fix property/atom for storing neighs IDs

  new_fix_id = utils::strdup(
      id + std::string("_FIX_PA"));    // This is the name of the new fix property/atom
  modify->add_fix(fmt::format("{} {} property/atom i2_neighs {} ghost yes", new_fix_id,
                              group->names[igroup], std::to_string(neighs_MAX)));
  // the modify command is creating the fix proptery/atom call

  int tmp1, tmp2;    //these are flag variables that will be returned by find_custom
  index = atom->find_custom("neighs", tmp1, tmp2);
}

/* ---------------------------------------------------------------------- */

void FixMorphoDynamic::init()
{ 
  //check EQUAL style variables
  if (nevery_output_str) {
    nevery_output_var = input->variable->find(nevery_output_str);
    if (nevery_output_var < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", nevery_output_str);
    if (input->variable->equalstyle(nevery_output_var))
      nevery_output_style = EQUAL;
  }
  if(nu0_str){
    nu0_var = input->variable->find(nu0_str);
    if(nu0_var < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", nu0_str);
    if(input->variable->equalstyle(nu0_var))
      nu0_style = EQUAL;
  }
  if(berdamp_str){
    berdamp_var = input->variable->find(berdamp_str);
    if(berdamp_var < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", berdamp_str);
    if(input->variable->equalstyle(berdamp_var))
      berdamp_style = EQUAL;
  }
  if(berflagx_str){
    berflagx_var = input->variable->find(berflagx_str);
    if(berflagx_var < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", berflagx_str);
    if(input->variable->equalstyle(berflagx_var))
      berflagx_style = EQUAL;
  }
  if(berflagy_str){
    berflagy_var = input->variable->find(berflagy_str);
    if(berflagy_var < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", berflagy_str);
    if(input->variable->equalstyle(berflagy_var))
      berflagy_style = EQUAL;
  }
  if (kpstr) {
    kpvar = input->variable->find(kpstr);
    if (kpvar < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", kpstr);
    if (input->variable->equalstyle(kpvar))
      kpstyle = EQUAL;
  }
  if (p0str) {
    p0var = input->variable->find(p0str);
    if (p0var < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", p0str);
    if (input->variable->equalstyle(p0var))
      p0style = EQUAL;
  }
  if (KTstr) {
    KTvar = input->variable->find(KTstr);
    if (KTvar < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", KTstr);
    if (input->variable->equalstyle(KTvar))
      KTstyle = EQUAL;
  }
  if (neveryMCstr) {
    neveryMCvar = input->variable->find(neveryMCstr);
    if (neveryMCvar < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", neveryMCstr);
    if (input->variable->equalstyle(neveryMCvar))
    neveryMCstyle = EQUAL;
  }
  if (Jsstr) {
    Jsvar = input->variable->find(Jsstr);
    if (Jsvar < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", Jsstr);
    if (input->variable->equalstyle(Jsvar))
      Jsstyle = EQUAL;
  }
  if (Jnstr) {
    Jnvar = input->variable->find(Jnstr);
    if (Jnvar < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", Jnstr);
    if (input->variable->equalstyle(Jnvar))
      Jnstyle = EQUAL;
  }
  if (fastr) {
    favar = input->variable->find(fastr);
    if (favar < 0) error->all(FLERR, "Variable {} for fix morphodynamic does not exist", fastr);
    if (input->variable->equalstyle(favar))
      fastyle = EQUAL;
  }

  // set index and check validity of region
  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix morphodynamic does not exist", idregion);
  }

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }

  /*From Fix berendsen*/

  rfix.clear();
  for (auto &ifix : modify->get_fix_list())
    if (ifix->rigid_flag) rfix.push_back(ifix);

  /*compute voronoi was giving weird results (not handling periodicity), hence moved to setup*/

  //We need a full neighbor list, built every Nevery steps

  // neighbor->add_request(this, NeighConst::REQ_FULL);
}

/*-----------------------------------------------------------------------------------------*/

// void FixMorphoDynamic::init_list(int /*id*/, NeighList *ptr)
// {
//   list = ptr;
// }

/*-----------------------------------------------------------------------------------------*/

void FixMorphoDynamic::setup(int vflag)
{
  //Read Equal style variables
  if (nevery_output_style == EQUAL){
    nevery_output = input->variable->compute_equal(nevery_output_var);
  }
  if (nu0_style == EQUAL){
    nu0 = input->variable->compute_equal(nu0_var);
  }
  if (berdamp_style == EQUAL){
    berdamp = input->variable->compute_equal(berdamp_var);
  }
  if (berflagx_style == EQUAL){
    berflagx = input->variable->compute_equal(berflagx_var);
  }
  if (berflagy_style == EQUAL){
    berflagy = input->variable->compute_equal(berflagy_var);
  }
  if (kpstyle == EQUAL){
    kp = input->variable->compute_equal(kpvar);
  }
  if (p0style == EQUAL){
    p0 = input->variable->compute_equal(p0var);
  }
  if (KTstyle == EQUAL){
    KT = input->variable->compute_equal(KTvar);
  }
  if (neveryMCstyle == EQUAL){
    neveryMC = input->variable->compute_equal(neveryMCvar);
  }
  if (Jsstyle == EQUAL){
    Js = input->variable->compute_equal(Jsvar);
  }
  if (Jnstyle == EQUAL){
    Jn = input->variable->compute_equal(Jnvar);
  }
  if (fastyle == EQUAL){
    fa = input->variable->compute_equal(favar);
  }

  if (flag_init == 0) {
    flag_init = 1;
    // Proc info
    double **x = atom->x;    //This is x_n
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    int nall = nlocal + nghost;
    tagint *tag = atom->tag;

    // Invoke compute
    modify->clearstep_compute();
    vcompute = modify->get_compute_by_id(id_compute_voronoi);

    //  int lf = vcompute->local_flag;
    //  int pf = vcompute->peratom_flag;

    /*For peratom array*/

    //   if (!(vcompute->invoked_flag & Compute::INVOKED_PERATOM)) {
    //     vcompute->compute_peratom();
    //     vcompute->invoked_flag |= Compute::INVOKED_PERATOM;
    //   }

    /*We only need local array for our purpose*/

    if (!(vcompute->invoked_flag & Compute::INVOKED_LOCAL)) {
      vcompute->compute_local();
      vcompute->invoked_flag |= Compute::INVOKED_LOCAL;
    }

    int num_rows_loc = vcompute->size_local_rows;

    /*Create neighs array (local) from vcompute->local array*/

    //pointer to neighs
    tagint **neighs = atom->iarray[index];

    // Initialize entries to 0 (space allocation for the only propert/atom array in this code)
    for (int i = 0; i < nall; i++) {
      for (int j = 0; j < neighs_MAX; j++) {
        if (mask[i] & groupbit) { neighs[i][j] = 0; }
      }
    }

    //Allocate sufficient memory as we don't know how many bonds in the system
    memory->create(bond_life, nlocal * 3, 3, "fix_morphodynamic:bond_life");
    num_bonds = 0;
    vector<vector<int>> bond_life_vec(nlocal * 3);

    // Populate neighs array with global ids of the voronoi neighs (initial configuration)
    // Also populate bond_life array (and vector)

    for (int n = 0; n < num_rows_loc; n++) {
      // skip rows with neighbor ID 0 as they denote z surfaces:
      if (int(vcompute->array_local[n][1]) == 0) { continue; }

      int cell_1, cell_2;

      if (int(vcompute->array_local[n][0]) < int(vcompute->array_local[n][1])) {
        cell_1 = int(vcompute->array_local[n][0]);
        cell_2 = int(vcompute->array_local[n][1]);
      } else if (int(vcompute->array_local[n][1]) < int(vcompute->array_local[n][0])) {
        cell_1 = int(vcompute->array_local[n][1]);
        cell_2 = int(vcompute->array_local[n][0]);
      }
      vector<int> target = {cell_1, cell_2, 0};
      if (find(bond_life_vec.begin(), bond_life_vec.end(), target) == bond_life_vec.end()) {
        bond_life_vec.push_back(target);
        bond_life[num_bonds][0] = cell_1;
        bond_life[num_bonds][1] = cell_2;
        bond_life[num_bonds][2] = 0;
        num_bonds += 1;
      }

      // get the local ID of atom
      int i = atom->map(int(vcompute->array_local[n][0]));
      if (i < 0) {
        error->one(FLERR, "Didn't find the atom");
        // Since the array is local, this should not get invoked
      }
      // Get the global ID of the neighbor of cell i
      int neigh_id = int(vcompute->array_local[n][1]);
      int m = 0;
      // Find the first empty space in the neighs array for atom i
      while (neighs[i][m] != 0) { m = m + 1; }
      if (mask[i] & groupbit) { neighs[i][m] = neigh_id; }
    }

    int num_neighs[nall] = {neighs_MAX};

    // Find num_neighs for nall cells
    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < neighs_MAX; j++) {
        if (neighs[i][j] == 0) {
          num_neighs[i] = j;
          break;
        }
      }
      int imj_atom = atom->sametag[i];
      while (imj_atom != -1) {
        num_neighs[imj_atom] = num_neighs[i];
        imj_atom = atom->sametag[imj_atom];
      }
    }

    // Arrange neighs in cyclic manner
    // While arranging this will tell you whether a neighbor has moved out the ghost cutoff or not
    // If thrown an exception, increase the cutoff
    for (int i = 0; i < nlocal; i++) { arrange_cyclic(neighs[i], num_neighs[i], i); }

    // Communicate neighs (to have ordered neighs list of ghost atoms as well)
    commflag = 2;
    comm->forward_comm(this, neighs_MAX);

    /*Create memeory allocations for other per atom info*/
    //allocate memory at once instead of every time step

    nmax = atom->nmax;
    memory->create(cell_shape, nmax, per_atom_out, "fix_morphodynamic:cell_shape");
    memory->create(vertices, nmax, neighs_MAX * 2, "fix_morphodynamic:vertices");

    // Initialize cell-shape arrays to zero
    for (int i = 0; i < nall; i++) {
      for (int j = 0; j < per_atom_out; j++) { cell_shape[i][j] = 0.0; }
    }

    // Initialize vertices info as well
    if (peratom_flag) {
      for (int i = 0; i < nall; i++) {
        for (int j = 0; j < neighs_MAX * 2; j++) { vertices[i][j] = 0.0; }
      }
    }

    //Calculate cell geometry data and populate vertices array (t = 0)
    for (int i = 0; i < nall; i++) {
      get_cell_data1(vertices[i], cell_shape[i], neighs[i], num_neighs[i], i);
    }

    //Calculate cell mechanics data (t = 0)
    for (int i = 0; i < nlocal; i++) {
      get_cell_data2(i, vertices[i], num_neighs[i], cell_shape, neighs[i]);
    }

    // Write 'files' at ts = 0

    fp1 = fopen(file1.c_str(), "w");
    fp2 = fopen(file2.c_str(), "w");
    fp3 = fopen(file3.c_str(), "w");
    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < num_neighs[i]; j++) {
        fprintf(fp1, "%d  ", neighs[i][j]);
        fprintf(fp2, "%f  %f  ", vertices[i][2 * j], vertices[i][2 * j + 1]);
      }
      fprintf(fp1, "\n");
      fprintf(fp2, "\n");
    }
    fprintf(fp1, "\n\n");
    fclose(fp1);
    fprintf(fp2, "\n\n");
    fclose(fp2);

    /*Don't need this verlet or respa stuff (probably!) as long as things are going well*/

    /*Initializing the neighs_array before post force method is called inside setup*/

    // if (utils::strmatch(update->integrate_style, "^verlet")){
    //   post_force(vflag);
    // }
    // else {
    //   (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    //   post_force_respa(vflag, ilevel_respa, 0);
    //   (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
    // }

    // run for zeroth step
    post_force(vflag);
  }
}

/* ---------------------------------------------------------------------- */
/*Not needed this as we are not doing any minimization and only doing dynamic runs*/

// void FixMorphoDynamic::min_setup(int vflag)
// {
//   post_force(vflag);
// }

/* --------------------------------------------------------------------------------- */

void FixMorphoDynamic::post_force(int vflag)
{
  //Read EQUAL style input variables
  if (nevery_output_style == EQUAL){
    nevery_output = input->variable->compute_equal(nevery_output_var);
  }
  if (nu0_style == EQUAL){
    nu0 = input->variable->compute_equal(nu0_var);
  }
  if (berdamp_style == EQUAL){
    berdamp = input->variable->compute_equal(berdamp_var);
  }
  if (berflagx_style == EQUAL){
    berflagx = input->variable->compute_equal(berflagx_var);
  }
  if (berflagy_style == EQUAL){
    berflagy = input->variable->compute_equal(berflagy_var);
  }
  if (kpstyle == EQUAL){
    kp = input->variable->compute_equal(kpvar);
  }
  if (p0style == EQUAL){
    p0 = input->variable->compute_equal(p0var);
  }
  if (KTstyle == EQUAL){
    KT = input->variable->compute_equal(KTvar);
  }
  if (neveryMCstyle == EQUAL){
    neveryMC = input->variable->compute_equal(neveryMCvar);
  }
  if (Jsstyle == EQUAL){
    Js = input->variable->compute_equal(Jsvar);
  }
  if (Jnstyle == EQUAL){
    Jn = input->variable->compute_equal(Jnvar);
  }
  if (fastyle == EQUAL){
    fa = input->variable->compute_equal(favar);
  }

  /*Read in current data*/

  double **x = atom->x;    //This is x_n
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  imageint *image = atom->image;
  tagint *tag = atom->tag;
  double dt = update->dt;
  int *sametag = atom->sametag;    //this returns the next local id of the atom having the same tag

  int natoms = atom->natoms;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  double *cut = comm->cutghost;

  if (update->ntimestep % nevery) return;

  // virial setup
  v_init(vflag);
  // update region if necessary
  if (region) region->prematch();

  int me = comm->me;    //current rank value
  int nprocs = comm->nprocs;

  // Possibly resize arrays
  if (atom->nmax > nmax) {
    memory->destroy(cell_shape);
    memory->destroy(vertices);
    nmax = atom->nmax;
    memory->create(cell_shape, nmax, per_atom_out, "fix_morphodynamic:cell_shape");
    memory->create(vertices, nmax, neighs_MAX * 2, "fix_morphodynamic:vertices");
  }

  tagint **neighs = atom->iarray[index];    //This is basically topology T^{n-1}

  int num_neighs[nall];
  for (int i = 0; i < nall; i++) {
    for (int j = 0; j < neighs_MAX; j++) {
      if (neighs[i][j] == 0) {
        num_neighs[i] = j;
        break;
      }
    }
  }

  /************************************************************************************************************* */
  /*  STEP 1 : First find updated/cyclically rarranged neighs list based on changes in position of cell centers*/
  /************************************************************************************************************* */

  tagint neighs_rearr[nlocal][neighs_MAX];

  if (update->ntimestep > 0) {

    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < neighs_MAX; j++) { neighs_rearr[i][j] = neighs[i][j]; }
      //arrange the neighs cyclically based on x*
      arrange_cyclic(neighs_rearr[i], num_neighs[i], i);
    }

    // Compare the original and rearranged neighs list to see if it is a cyclic permutation or not

    for (int i = 0; i < nlocal; i++) {
      if (is_cyclic_perm(neighs[i], neighs_rearr[i], num_neighs[i])) {
        continue;
      } else {
        error->one(FLERR,
                   "\n !!!! The new neighs list (after cells positions have been updated) is not a "
                   "cyclic permutation of the previous one !!!! \n");
      }
    }

    // Update neighs list (now cyclically permuted based on x_updated)

    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < neighs_MAX; j++) { neighs[i][j] = neighs_rearr[i][j]; }
    }

    // Communicate neighs (to have ordered neighs list of ghost atoms as well)
    commflag = 2;
    comm->forward_comm(this, neighs_MAX);

    // Now that cells have moved, tiling has changed. Update cell-shape data based on x_updated
    // Num of neighs has not changed

    for (int i = 0; i < nall; i++) {
      get_cell_data1(vertices[i], cell_shape[i], neighs[i], num_neighs[i], i);
    }

    // Update cell stress
    for (int i = 0; i < nlocal; i++) {
      get_cell_data2(i, vertices[i], num_neighs[i], cell_shape, neighs[i]);
    }

    //CHECK for self intersecting polygons

    for (int i = 0; i < nlocal; i++) {
      if (isSelfIntersecting(vertices[i], 2 * num_neighs[i])) {
        error->one(FLERR, " !!!!!! Self intersecting polygon %d still detected !!!!! \n", tag[i]);
      }
    }
  }

  //DEBUGGER
  double etotal[num_MC] = {0.0};
  double etotal_0 = 0.0;
  for (int i = 0; i < nlocal; i++) { etotal_0 += cell_shape[i][5]; }
  //DEBUGGER
  /***********************************************************************************************************************/
  /*    STEP 2 : Perform cell rearrangements (T1 transitions) based on Metropolis algorithm (topology Tn = Tn(Tn-1, xn)*/
  /***********************************************************************************************************************/

  int num_cellrearr = 0;
  //Reinitalize age of all bonds (currently in the network) to 0 before MC
  for (int i = 0; i < num_bonds; i++) { bond_life[i][2] = 0; }

  //Start Monte-Carlo
  //random seed for each time. You can use srand(update->ntimestep) time(NULL) + me for a more 'stochastic' like behavior
  //srand(update->ntimestep);
  srand(time(NULL) + me); 

  if ((update->ntimestep > 0) && (update->ntimestep % neveryMC == 0)) {
    for (int n = 0; n < num_MC; n++) {

      //Randomly pick a bond
      int bond_id = rand() % num_bonds;    //returns a random integer between [0, num_bonds-1]

      //allocate memory for storing tags of cells in the Quad in this mc step and some other info
      tagint cell_tags[4] = {0};

      int cell_i = atom->map(bond_life[bond_id][0]);    //local ids of cells forming bond
      int cell_j = atom->map(bond_life[bond_id][1]);

      //error check for num_neighs_celli-->behaving a bit wierdly! when full
      if (num_neighs[cell_i] == 0) { error->one(FLERR, "no of neighs for atom is 0"); }

      // //CHECK 1 --> Bond must be an old one and not a new one
      // if (bond_life[bond_id][2] == 1) {
      //   //OUTPUTTING
      //   if (n == 0) {
      //     etotal[n] = etotal_0;
      //   } else {
      //     etotal[n] = etotal[n - 1];
      //   }
      //   //OUTPUTTING
      //   continue;
      // }

      //CHECK 2 ---> No of neighs for a cell should not fall below 4
      if (num_neighs[cell_i] <= 4 || num_neighs[cell_j] <= 4) {    //OUTPUTTING
        if (n == 0) {
          etotal[n] = etotal_0;
        } else {
          etotal[n] = etotal[n - 1];
        }
        //OUTPUTTING
        continue;
      }

      tagint comm_neigh_ij[2] = {-1};    // array to store common pair atoms
      int num_comm_neighs_ij = get_comm_neigh(comm_neigh_ij, neighs[cell_i], neighs[cell_j],
                                              num_neighs[cell_i], num_neighs[cell_j]);

      //CHECK 3 --> Two bonded atoms (i and j) should have exactly two distinct common neighs
      if (num_comm_neighs_ij != 2) {
        //OUTPUTTING
        if (n == 0){
          etotal[n] = etotal_0; 
        }
        else{
          etotal[n] = etotal[n - 1]; 
        }        
        //OUTPUTTING
        continue;
        // error->one(FLERR,
        //            " \n No of common neighs between bonded atoms (i and j) not exactly 2 --> "
        //            "probably some "
        //            "problem in swapping , This could also happen if there is a triangular cell "
        //            "involved\n");
      }

      if (comm_neigh_ij[0] == comm_neigh_ij[1]) {
        error->one(
            FLERR,
            " \n Both the comm neighs of i and j are same --> might be due to some atom added "
            "twice \n");
      }

      //If common neighs for cell i and j are exactly 2 and distinct then continue

      //this gives the local ids of the common neighs of cell i and j
      int cell_k = atom->map(comm_neigh_ij[0]);
      int cell_l = atom->map(comm_neigh_ij[1]);

      //all cell_i,j,k,l are 'local' ids at this point

      cell_tags[0] = tag[cell_i];
      cell_tags[1] = tag[cell_j];
      cell_tags[2] = tag[cell_k];
      cell_tags[3] = tag[cell_l];

      //Find the nearest images of cells from cell_i (now based on updated positions)

      int cell_j_nearest = domain->closest_image(cell_i, cell_j);
      int cell_k_nearest = domain->closest_image(cell_i, cell_k);
      int cell_l_nearest = domain->closest_image(cell_i, cell_l);

      //CHECK 4 --> Get the correct images of cells and their neighbors

      int flag_valid_imj = 0;

      int closest_imj_of_l_from_k = domain->closest_image(cell_k, cell_l);
      int closest_imj_of_k_from_l = domain->closest_image(cell_l, cell_k);

      if (cell_k_nearest == cell_k &&
          cell_l_nearest ==
              cell_l) {    // if k and l are owned then the nearest images must coincide
        if (cell_k_nearest == closest_imj_of_k_from_l &&
            cell_l_nearest == closest_imj_of_l_from_k) {
          flag_valid_imj = 1;    //swapping allowed as nearest images coincide
        }
      } else if (cell_k_nearest != cell_k &&
                 cell_l_nearest !=
                     cell_l) {    //if both k and l are ghost then nearest images must be different
        if (cell_k_nearest != closest_imj_of_k_from_l &&
            cell_l_nearest != closest_imj_of_l_from_k) {
          flag_valid_imj = 1;    //swapping allowed as nearest images don't coincide
        }
      } else if (cell_k_nearest == cell_k &&
                 cell_l_nearest != cell_l) {    //if k is owned and l is ghost
        if (cell_l_nearest == closest_imj_of_l_from_k &&
            cell_k_nearest != closest_imj_of_k_from_l) {
          flag_valid_imj = 1;
        }
      } else if (cell_k_nearest != cell_k && cell_l_nearest == cell_l) {
        if (cell_k_nearest == closest_imj_of_k_from_l &&
            cell_l_nearest != closest_imj_of_l_from_k) {
          flag_valid_imj = 1;
        }
      }

      if (flag_valid_imj == 0) {
        //OUTPUTTING
        if (n == 0) {
          etotal[n] = etotal_0;
        } else {
          etotal[n] = etotal[n - 1];
        }
        //OUTPUTTING
        continue;    //not a valid swap
      }

      //CHECK 5 --> see if cell_k and cell_l are non-bonded

      int flag_unbonded = unbonded(neighs[cell_k], num_neighs[cell_k], tag[cell_l]);

      if (flag_unbonded == 0) {
        error->one(FLERR, "cells k and l are already bonded but should not be!!!!");
      }

      tagint comm_neigh_kl[3] = {-1};    // array to store common pair atoms
      int num_comm_neighs_kl = get_comm_neigh(comm_neigh_kl, neighs[cell_k], neighs[cell_l],
                                              num_neighs[cell_k], num_neighs[cell_l]);

      if (num_comm_neighs_kl != 2) {
        //They should not form a bond as they have more than two common neighs
        //This is a perfectly valid case (don't worry)
        //OUTPUTTING
        if (n == 0) {
          etotal[n] = etotal_0;
        } else {
          etotal[n] = etotal[n - 1];
        }
        //OUTPUTTING
        continue;
      }

      //CHECK 6 --> See if the quadrialteral formed is convex or concave:
      //i->k->j->l is the order of points of quadralteral in a cyclic manner

      double p1[2] = {x[cell_i][0], x[cell_i][1]};
      double p2[2] = {x[cell_k_nearest][0], x[cell_k_nearest][1]};
      double p3[2] = {x[cell_j_nearest][0], x[cell_j_nearest][1]};
      double p4[2] = {x[cell_l_nearest][0], x[cell_l_nearest][1]};

      //If quadrilateral formed is concave then resultant tiling not valid
      if (isConcave(p1, p2, p3, p4)) {
        //OUTPUTTING
        if (n == 0) {
          etotal[n] = etotal_0;
        } else {
          etotal[n] = etotal[n - 1];
        }
        //OUTPUTTING
        continue;
      }

      /*********** All check completed ---> Attempt T1 swapping *************/

      /*STEP 1 ---> Calculate energy of the transition state*/

      //Find the transition state vertex

      double vertex_ijl[2], vertex_ijk[2], vertex_ikl[2], vertex_jkl[2];

      vertex_ijl[0] = (x[cell_i][0] + x[cell_j_nearest][0] + x[cell_l_nearest][0]) / 3.0;
      vertex_ijl[1] = (x[cell_i][1] + x[cell_j_nearest][1] + x[cell_l_nearest][1]) / 3.0;

      vertex_ijk[0] = (x[cell_i][0] + x[cell_j_nearest][0] + x[cell_k_nearest][0]) / 3.0;
      vertex_ijk[1] = (x[cell_i][1] + x[cell_j_nearest][1] + x[cell_k_nearest][1]) / 3.0;

      vertex_ikl[0] = (x[cell_i][0] + x[cell_k_nearest][0] + x[cell_l_nearest][0]) / 3.0;
      vertex_ikl[1] = (x[cell_i][1] + x[cell_k_nearest][1] + x[cell_l_nearest][1]) / 3.0;

      vertex_jkl[0] = (x[cell_j_nearest][0] + x[cell_k_nearest][0] + x[cell_l_nearest][0]) / 3.0;
      vertex_jkl[1] = (x[cell_j_nearest][1] + x[cell_k_nearest][1] + x[cell_l_nearest][1]) / 3.0;

      double vertex_trans[2];
      //Take the average of vertices of all four triangles involved
      vertex_trans[0] = (vertex_ijl[0] + vertex_ijk[0] + vertex_ikl[0] + vertex_jkl[0]) / 4.0;
      vertex_trans[1] = (vertex_ijl[1] + vertex_ijk[1] + vertex_ikl[1] + vertex_jkl[1]) / 4.0;

      // Create the deformed polygonal shapes for cells i, j, k and l

      double dummy_vertices_i[num_neighs[cell_i] - 1][2];
      double dummy_vertices_j[num_neighs[cell_j_nearest] - 1][2];
      double dummy_vertices_k[num_neighs[cell_k_nearest]][2];
      double dummy_vertices_l[num_neighs[cell_l_nearest]][2];

      //Loop through the vertices of each cell and identify the vertices that needs to be changed/swapped

      double trans_state_geom[4][4] = {0.0};
      double epsilon = 1e-5;

      // For cell i

      int num_dummyvertices = 0;

      for (int i = 0; i < num_neighs[cell_i]; i++) {
        if (abs(vertices[cell_i][2 * i] - vertex_ijk[0]) > epsilon ||
            abs(vertices[cell_i][2 * i + 1] - vertex_ijk[1]) > epsilon) {
          if (abs(vertices[cell_i][2 * i] - vertex_ijl[0]) > epsilon ||
              abs(vertices[cell_i][2 * i + 1] - vertex_ijl[1]) > epsilon) {
            dummy_vertices_i[num_dummyvertices][0] = vertices[cell_i][2 * i];
            dummy_vertices_i[num_dummyvertices][1] = vertices[cell_i][2 * i + 1];
            num_dummyvertices += 1;
          }
        }
      }

      if (num_dummyvertices != num_neighs[cell_i] - 2) {
        error->one(FLERR, "vertices not identified correctly for i!");
      }

      dummy_vertices_i[num_dummyvertices][0] = vertex_trans[0];
      dummy_vertices_i[num_dummyvertices][1] = vertex_trans[1];

      sortverticesCCW(dummy_vertices_i, num_neighs[cell_i] - 1, x[cell_i][0], x[cell_i][1]);

      get_trans_state_geom(trans_state_geom[0], dummy_vertices_i, num_neighs[cell_i] - 1);

      // For cell j

      num_dummyvertices = 0;

      for (int i = 0; i < num_neighs[cell_j_nearest]; i++) {
        if (abs(vertices[cell_j_nearest][2 * i] - vertex_ijk[0]) > epsilon ||
            abs(vertices[cell_j_nearest][2 * i + 1] - vertex_ijk[1]) > epsilon) {
          if (abs(vertices[cell_j_nearest][2 * i] - vertex_ijl[0]) > epsilon ||
              abs(vertices[cell_j_nearest][2 * i + 1] - vertex_ijl[1]) > epsilon) {
            dummy_vertices_j[num_dummyvertices][0] = vertices[cell_j_nearest][2 * i];
            dummy_vertices_j[num_dummyvertices][1] = vertices[cell_j_nearest][2 * i + 1];
            num_dummyvertices += 1;
          }
        }
      }

      if (num_dummyvertices != num_neighs[cell_j_nearest] - 2) {
        error->one(FLERR, "vertices not identified correctly for j!");
      }

      dummy_vertices_j[num_dummyvertices][0] = vertex_trans[0];
      dummy_vertices_j[num_dummyvertices][1] = vertex_trans[1];

      sortverticesCCW(dummy_vertices_j, num_neighs[cell_j_nearest] - 1, x[cell_j_nearest][0],
                      x[cell_j_nearest][1]);

      get_trans_state_geom(trans_state_geom[1], dummy_vertices_j, num_neighs[cell_j_nearest] - 1);

      // For cell k

      num_dummyvertices = 0;

      for (int i = 0; i < num_neighs[cell_k_nearest]; i++) {
        if (abs(vertices[cell_k_nearest][2 * i] - vertex_ijk[0]) > epsilon ||
            abs(vertices[cell_k_nearest][2 * i + 1] - vertex_ijk[1]) > epsilon) {
          dummy_vertices_k[num_dummyvertices][0] = vertices[cell_k_nearest][2 * i];
          dummy_vertices_k[num_dummyvertices][1] = vertices[cell_k_nearest][2 * i + 1];
          num_dummyvertices += 1;
        }
      }

      if (num_dummyvertices != num_neighs[cell_k_nearest] - 1) {
        error->one(FLERR, "vertices not identified correctly for k!");
      }

      dummy_vertices_k[num_dummyvertices][0] = vertex_trans[0];
      dummy_vertices_k[num_dummyvertices][1] = vertex_trans[1];

      sortverticesCCW(dummy_vertices_k, num_neighs[cell_k_nearest], x[cell_k_nearest][0],
                      x[cell_k_nearest][1]);

      get_trans_state_geom(trans_state_geom[2], dummy_vertices_k, num_neighs[cell_k_nearest]);

      // For cell l

      num_dummyvertices = 0;

      for (int i = 0; i < num_neighs[cell_l_nearest]; i++) {
        if (abs(vertices[cell_l_nearest][2 * i] - vertex_ijl[0]) > epsilon ||
            abs(vertices[cell_l_nearest][2 * i + 1] - vertex_ijl[1]) > epsilon) {
          dummy_vertices_l[num_dummyvertices][0] = vertices[cell_l_nearest][2 * i];
          dummy_vertices_l[num_dummyvertices][1] = vertices[cell_l_nearest][2 * i + 1];
          num_dummyvertices += 1;
        }
      }

      if (num_dummyvertices != num_neighs[cell_l_nearest] - 1) {
        error->one(FLERR, "vertices not identified correctly for l!");
      }

      dummy_vertices_l[num_dummyvertices][0] = vertex_trans[0];
      dummy_vertices_l[num_dummyvertices][1] = vertex_trans[1];

      sortverticesCCW(dummy_vertices_l, num_neighs[cell_l_nearest], x[cell_l_nearest][0],
                      x[cell_l_nearest][1]);

      get_trans_state_geom(trans_state_geom[3], dummy_vertices_l, num_neighs[cell_l_nearest]);

      // Get transition state (TS) energy
      double TS_energy[4] = {0.0};

      // for cell i
      TS_energy[0] = get_trans_state_energy(0, neighs[cell_i], num_neighs[cell_i], trans_state_geom,
                                            cell_shape, cell_tags);
      TS_energy[1] = get_trans_state_energy(1, neighs[cell_j], num_neighs[cell_j], trans_state_geom,
                                            cell_shape, cell_tags);
      TS_energy[2] = get_trans_state_energy(2, neighs[cell_k], num_neighs[cell_k], trans_state_geom,
                                            cell_shape, cell_tags);
      TS_energy[3] = get_trans_state_energy(3, neighs[cell_l], num_neighs[cell_l], trans_state_geom,
                                            cell_shape, cell_tags);

      double E_TS = TS_energy[0] + TS_energy[1] + TS_energy[2] + TS_energy[3];

      //store energy before swapping attempt

      double E_before = cell_shape[cell_i][5] + cell_shape[cell_j][5] + cell_shape[cell_k][5] +
          cell_shape[cell_l][5];

      /*STEP 2 ---> Attempt a bond swap*/

      //Update neighs information first for all owned atoms (and their ghost images) involved

      //Delete bond
      for (int m = 0; m < 2; m++) {
        int owned_atom = atom->map(cell_tags[m]);
        tagint neigh_tag;
        if (m == 0) {
          neigh_tag = cell_tags[1];
        } else if (m == 1) {
          neigh_tag = cell_tags[0];
        }
        //First update info on owned atom
        remove_neigh(owned_atom, neigh_tag, neighs[owned_atom], num_neighs[owned_atom], neighs_MAX);
        num_neighs[owned_atom] -= 1;
        if (num_neighs[owned_atom] == 3) { error->one(FLERR, "number of neighs fell below 4"); }
        arrange_cyclic(neighs[owned_atom], num_neighs[owned_atom], owned_atom);
        //Now update info on its ghost images
        int imj_atom = sametag[owned_atom];
        while (imj_atom != -1) {
          for (int ii = 0; ii < neighs_MAX; ii++) { neighs[imj_atom][ii] = neighs[owned_atom][ii]; }
          num_neighs[imj_atom] = num_neighs[owned_atom];
          imj_atom = sametag[imj_atom];
        }
      }

      //Create bond
      for (int m = 2; m < 4; m++) {
        int owned_atom = atom->map(cell_tags[m]);
        tagint neigh_tag;
        if (m == 2) {
          neigh_tag = cell_tags[3];
        } else if (m == 3) {
          neigh_tag = cell_tags[2];
        }
        //Update info of owned atom
        neighs[owned_atom][num_neighs[owned_atom]] = neigh_tag;
        num_neighs[owned_atom] += 1;
        if (num_neighs[owned_atom] == neighs_MAX) {
          printf("\n Limit reached for adding new neighs to atom %d\n", tag[owned_atom]);
          error->one(FLERR, "EXITED");
        }
        arrange_cyclic(neighs[owned_atom], num_neighs[owned_atom], owned_atom);
        //Now update neighs list of ghost images
        int imj_atom = sametag[owned_atom];
        while (imj_atom != -1) {
          for (int ii = 0; ii < neighs_MAX; ii++) { neighs[imj_atom][ii] = neighs[owned_atom][ii]; }
          num_neighs[imj_atom] = num_neighs[owned_atom];
          imj_atom = sametag[imj_atom];
        }
      }

      //Now get the cell_geometry data for all 4 (their ghost images as well) atoms after updated topology

      for (int m = 0; m < 4; m++) {
        int owned_atom = atom->map(cell_tags[m]);
        get_cell_data1(vertices[owned_atom], cell_shape[owned_atom], neighs[owned_atom],
                       num_neighs[owned_atom], owned_atom);
        int imj_atom = sametag[owned_atom];
        while (imj_atom != -1) {
          get_cell_data1(vertices[imj_atom], cell_shape[imj_atom], neighs[imj_atom],
                         num_neighs[imj_atom], imj_atom);
          imj_atom = sametag[imj_atom];
        }
      }

      //Now get the cell mechanics data for all 4 atoms after updated topology

      for (int m = 0; m < 4; m++) {
        int owned_atom = atom->map(cell_tags[m]);
        get_cell_data2(owned_atom, vertices[owned_atom], num_neighs[owned_atom], cell_shape,
                       neighs[owned_atom]);
        int imj_atom = sametag[owned_atom];
        while (imj_atom != -1) {
          for (int ii = 5; ii < 10; ii++) { cell_shape[imj_atom][ii] = cell_shape[owned_atom][ii]; }
          imj_atom = sametag[imj_atom];
        }
      }

      // Check to see if any of the new polygons formed is self-intersecting or not

      int flag_self_intersect = 1;

      for (int i = 0; i < 4; i++) {
        int cell = atom->map(cell_tags[i]);
        if (isSelfIntersecting(vertices[cell], 2 * num_neighs[cell])) {
          flag_self_intersect = 0;
          break;
        }
      }

      // Check if the T1 swap results in valid triangulation or not

      int flag_valid_tri = 1;

      for (int i = 0; i < 4; i++) {
        int cell = atom->map(cell_tags[i]);
        int num_neighs_cell = num_neighs[cell];
        int neigh_prev, neigh_next;
        for (int j = 0; j < num_neighs_cell; j++) {
          int curr_neigh = atom->map(neighs[cell][j]);    //local id
          if (j == 0) {
            neigh_prev = neighs[cell][num_neighs_cell - 1];
            neigh_next = neighs[cell][j + 1];
          } else if (j == num_neighs_cell - 1) {
            neigh_prev = neighs[cell][j - 1];
            neigh_next = neighs[cell][0];
          } else {
            neigh_prev = neighs[cell][j - 1];
            neigh_next = neighs[cell][j + 1];
          }
          //check the previous and next atoms for the current neigh j
          int num_neighs_curr = num_neighs[curr_neigh];
          int neigh_prev_curr, neigh_next_curr;
          int jj = 0;
          while (tag[cell] != neighs[curr_neigh][jj]) {
            jj++;
            // error check
            if (jj == num_neighs_curr) {
              error->one(FLERR, "\n tag[cell] not found in the neighs list of curr neigh \n");
            }
            //error check
          }
          if (jj == 0) {
            neigh_prev_curr = neighs[curr_neigh][num_neighs_curr - 1];
            neigh_next_curr = neighs[curr_neigh][jj + 1];
          } else if (jj == num_neighs_curr - 1) {
            neigh_prev_curr = neighs[curr_neigh][jj - 1];
            neigh_next_curr = neighs[curr_neigh][0];
          } else {
            neigh_prev_curr = neighs[curr_neigh][jj - 1];
            neigh_next_curr = neighs[curr_neigh][jj + 1];
          }
          if (neigh_prev != neigh_next_curr || neigh_next != neigh_prev_curr) {
            flag_valid_tri = 0;
            break;
          }
        }
        if (flag_valid_tri == 0) { break; }
      }

      // If all criteria are satisfied then check the energy criteria for swap

      double E_after = cell_shape[cell_i][5] + cell_shape[cell_j][5] + cell_shape[cell_k][5] +
          cell_shape[cell_l][5];

      // find the Probability of swapping

      double P_swap;

      if ((E_before >= E_TS) && (E_before >= E_after)) {
        P_swap = 1.0;
      } else if ((E_TS >= E_before) && (E_TS >= E_after)) {
        P_swap = exp(-1.0 * (E_TS - E_before) / KT);
      } else if ((E_after >= E_TS) && (E_after >= E_before)) {
        P_swap = exp(-1.0 * (E_after - E_before) / KT);
      }

      int flag_energy_crit = 1;                        //for energy criteria
      double rand_num = (double) rand() / RAND_MAX;    //generate a random no between 0 and 1

      if (flag_valid_tri == 1 && flag_self_intersect == 1) {
        //Metropolis scheme
        if (P_swap >= rand_num) {
          //Update bond_life array with cell k and l (and their lifetime)
          bond_life[bond_id][0] = cell_tags[2];
          bond_life[bond_id][1] = cell_tags[3];
          bond_life[bond_id][2] = 1;
          //update info for all the atoms
          // for (int ii = 0; ii < nall; ii++) {
          //   get_cell_data1(vertices[ii], cell_shape[ii], neighs[ii], num_neighs[ii], ii);
          // }      
          // Update cell stress
          for (int ii = 0; ii < nlocal; ii++) {
            get_cell_data2(ii, vertices[ii], num_neighs[ii], cell_shape, neighs[ii]);
          }
          num_cellrearr += 1;
          //OUTPUTTING
          if (n == 0) {
            etotal[n] = etotal_0 + E_after - E_before;
          } else {
            etotal[n] = etotal[n - 1] + E_after - E_before;
          }
          //OUTPUTTING
          continue;
        } else {
          flag_energy_crit = 0;    //swapping attempt failed as not energetically favoured
        }
      }

      //Reverse swapping if not fulfilled all the criteria

      if ((flag_energy_crit == 0) || (flag_valid_tri == 0) || (flag_self_intersect == 0)) {
        //OUTPUTTING
        if (n == 0) {
          etotal[n] = etotal_0;
        } else {
          etotal[n] = etotal[n - 1];
        }
        //OUTPUTTING

        cell_tags[0] = tag[cell_k];
        cell_tags[1] = tag[cell_l];
        cell_tags[2] = tag[cell_i];
        cell_tags[3] = tag[cell_j];

        //First delete bond between atom k and l (only if they were previously unbonded)
        //Delete bond
        for (int m = 0; m < 2; m++) {
          int owned_atom = atom->map(cell_tags[m]);
          tagint neigh_tag;
          if (m == 0) {
            neigh_tag = cell_tags[1];
          } else if (m == 1) {
            neigh_tag = cell_tags[0];
          }
          //First update info on owned atom
          remove_neigh(owned_atom, neigh_tag, neighs[owned_atom], num_neighs[owned_atom],
                       neighs_MAX);
          num_neighs[owned_atom] -= 1;
          if (num_neighs[owned_atom] == 3) { error->one(FLERR, "number of neighs fell below 4"); }
          arrange_cyclic(neighs[owned_atom], num_neighs[owned_atom], owned_atom);
          //Now update info on its ghost images
          int imj_atom = sametag[owned_atom];
          while (imj_atom != -1) {
            for (int ii = 0; ii < neighs_MAX; ii++) {
              neighs[imj_atom][ii] = neighs[owned_atom][ii];
            }
            num_neighs[imj_atom] = num_neighs[owned_atom];
            imj_atom = sametag[imj_atom];
          }
        }

        //Then Create bond between i and j
        for (int m = 2; m < 4; m++) {
          int owned_atom = atom->map(cell_tags[m]);
          tagint neigh_tag;
          if (m == 2) {
            neigh_tag = cell_tags[3];
          } else if (m == 3) {
            neigh_tag = cell_tags[2];
          }
          //Update info of owned atom
          neighs[owned_atom][num_neighs[owned_atom]] = neigh_tag;
          num_neighs[owned_atom] += 1;
          if (num_neighs[owned_atom] == neighs_MAX) {
            printf("\n Limit reached for adding new neighs to atom %d\n", tag[owned_atom]);
            error->one(FLERR, "EXITED");
          }
          arrange_cyclic(neighs[owned_atom], num_neighs[owned_atom], owned_atom);
          //Now update info on ghost images
          int imj_atom = sametag[owned_atom];
          while (imj_atom != -1) {
            for (int ii = 0; ii < neighs_MAX; ii++) {
              neighs[imj_atom][ii] = neighs[owned_atom][ii];
            }
            num_neighs[imj_atom] = num_neighs[owned_atom];
            imj_atom = sametag[imj_atom];
          }
        }

        // Again go back to original information
        for (int m = 0; m < 4; m++) {
          int owned_atom = atom->map(cell_tags[m]);
          get_cell_data1(vertices[owned_atom], cell_shape[owned_atom], neighs[owned_atom],
                         num_neighs[owned_atom], owned_atom);
          int imj_atom = sametag[owned_atom];
          while (imj_atom != -1) {
            get_cell_data1(vertices[imj_atom], cell_shape[imj_atom], neighs[imj_atom],
                           num_neighs[imj_atom], imj_atom);
            imj_atom = sametag[imj_atom];
          }
        }

        for (int m = 0; m < 4; m++) {
          int owned_atom = atom->map(cell_tags[m]);
          get_cell_data2(owned_atom, vertices[owned_atom], num_neighs[owned_atom], cell_shape,
                         neighs[owned_atom]);
          int imj_atom = sametag[owned_atom];
          while (imj_atom != -1) {
            for (int ii = 5; ii < 10; ii++) {
              cell_shape[imj_atom][ii] = cell_shape[owned_atom][ii];
            }
            imj_atom = sametag[imj_atom];
          }
        }
      }
    }
  }

  /*We have made sure that ghost atoms info is correct so don't need communication now 
  Might help save some time and make code a bit faster*/

  // // DEBUGGER
  // fp[1] << endl << endl;
  // // DEBUGGER

  /*Write the output txt "files" for ts > 0*/
  if ((update->ntimestep > 0) && (update->ntimestep % nevery_output == 0)) {
    fp1 = fopen(file1.c_str(), "a");
    fp2 = fopen(file2.c_str(), "a");
    fp3 = fopen(file3.c_str(), "a");
    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < num_neighs[i]; j++) {
        fprintf(fp1, "%d  ", neighs[i][j]);
        fprintf(fp2, "%f  %f  ", vertices[i][2 * j], vertices[i][2 * j + 1]);
      }
      fprintf(fp1, "\n");
      fprintf(fp2, "\n");
    }
    fprintf(fp1, "\n\n");
    fclose(fp1);
    fprintf(fp2, "\n\n");
    fclose(fp2);
    fprintf(fp3, "%f   %d \n", update->ntimestep * update->dt, num_cellrearr);
    // for (int n = 0; n < num_MC; n++){
    //   fprintf(fp3, " %f \n", etotal[n]);
    // }
    fclose(fp3);
  }

  /************************************************************************************************************* */
  /*               Initial relaxation step : Apply fix berendsen using current pressure p(xn, Tn)
                   and then remap coordinates xn to xn* based on changed box dimensions
                      again find cell data based on remapped positions 
                      then find forces fn based on (xn*, Tn) and new cell data (area and peri)     */
  /************************************************************************************************************* */  
  //Apply fix berendesn for timesteps > begintimestep_berendsen

  // take press_flag as input which will be 1 if berendsen is to be applied

  if (berflagx == 1 || berflagy == 1) {
    
    //Let us implemenent a simpler version
      double sigma_xx = 0.0;
      double sigma_yy = 0.0;
      double area_sum = 0.0;
      for (int i = 0; i < nlocal; i++) {
        area_sum += cell_shape[i][0];
        sigma_xx += cell_shape[i][0] * cell_shape[i][6];
        sigma_yy += cell_shape[i][0] * cell_shape[i][7];
      }
      press_current[0] = sigma_xx / area_sum;
      press_current[1] = sigma_yy / area_sum;

      //Update x dimensions
      if (berflagx == 1) {
        dilation_MDN[0] = 1 + update->dt / berdamp * (press_stop[0] - press_current[0]);
      } else {
        dilation_MDN[0] = 1.0;
      }

      //Update y dimensions
      if (berflagy == 1) {
        dilation_MDN[1] = 1 + update->dt / berdamp * (press_stop[1] - press_current[1]);
      } else {
        dilation_MDN[1] = 1.0;
      }

    resize_box();

    // find updated neighs list based on new box dimensions

    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < neighs_MAX; j++) { neighs_rearr[i][j] = neighs[i][j]; }
      //arrange the neighs cyclically based on x*
      arrange_cyclic(neighs_rearr[i], num_neighs[i], i);
    }

    // Compare the original and rearranged neighs list to see if it is a cyclic permutation or not

    for (int i = 0; i < nlocal; i++) {
      if (is_cyclic_perm(neighs[i], neighs_rearr[i], num_neighs[i])) {
        continue;
      } else {
        print_neighs_list(neighs_rearr[i], num_neighs[i], i);
        print_neighs_list(neighs[i], num_neighs[i], i);
        error->one(FLERR,
                   "\n !!!! The new neighs list (after fix berendsen) is not a "
                   "cyclic permutation of the previous one !!!! \n");
      }
    }

    // Update neighs list (now cyclically permuted based on x_updated)

    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < neighs_MAX; j++) { neighs[i][j] = neighs_rearr[i][j]; }
    }

    // Communicate neighs (to have ordered neighs list of ghost atoms as well)
    commflag = 2;
    comm->forward_comm(this, neighs_MAX);

    // Now that cells have moved, tiling has changed. Update cell-shape data based on x_updated
    // Num of neighs has not changed

    for (int i = 0; i < nall; i++) {
      get_cell_data1(vertices[i], cell_shape[i], neighs[i], num_neighs[i], i);
    }

    // Update cell stress
    for (int i = 0; i < nlocal; i++) {
      get_cell_data2(i, vertices[i], num_neighs[i], cell_shape, neighs[i]);
    }

    //CHECK for self intersecting polygons

    for (int i = 0; i < nlocal; i++) {
      if (isSelfIntersecting(vertices[i], 2 * num_neighs[i])) {
        error->one(FLERR, " !!!!!! Self intersecting polygon %d detected after fix berendsen !!!!! \n", tag[i]);
      }
    }
  }
  /************************************************************************************************************* */
  /*               STEP 3 : find forces f_n on atoms now using x_n* and T_n               */
  /************************************************************************************************************* */

  /*Access neighbors list for pairwise soft repulsive forces*/

  // int inum = list->inum;
  // int *ilist = list->ilist;
  // int *numneigh = list->numneigh;
  // int **firstneigh = list->firstneigh;

  // //DEBUGGER
  // for (int i = 0; i < nlocal; i++) {
  //   int jnum = numneigh[i];
  //   printf("neigh list for cell %d is: ", tag[i]);
  //   for (int j = 0; j < jnum; j++) { printf("%d,  ", tag[firstneigh[i][j]]); }
  //   printf("\n");
  // }
  // //DEBUGGER

  double Jac = 1.0 / 3.0;

  double nu[nlocal] = {nu0};    //initializes everthing to eta_0
  double fx[nlocal] = {0.0};
  double fy[nlocal] = {0.0};

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      //coords of cell i
      double x0 = x[i][0];
      double y0 = x[i][1];

      double F_t1[2] = {0.0};     //term 1 of vertex/voronoi model
      double F_t2[2] = {0.0};     //term 2 of vertex/voronoi model
      double F_rep[2] = {0.0};    //soft repulsion forces

      //Force contri from neighs
      for (int j = 0; j < num_neighs[i]; j++) {

        int current_neigh = domain->closest_image(i, atom->map(neighs[i][j]));

        //coords of cell j
        double x1 = x[current_neigh][0];
        double y1 = x[current_neigh][1];

        //soft repulsion force
        double rij = sqrt(pow(x1 - x0, 2.0) + pow(y1 - y0, 2.0));
        if (rij < rcut) {
          F_rep[0] += ksoft * (rcut - rij) * (x0 - x1);
          F_rep[1] += ksoft * (rcut - rij) * (y0 - y1);
        }

        int num_neighsj = num_neighs[current_neigh];

        int cyclic_neighsj[num_neighsj + 4];

        cyclic_neighsj[0] = neighs[current_neigh][num_neighsj - 2];
        cyclic_neighsj[1] = neighs[current_neigh][num_neighsj - 1];

        for (int n = 2; n < num_neighsj + 2; n++) {
          cyclic_neighsj[n] = neighs[current_neigh][n - 2];
        }

        cyclic_neighsj[num_neighsj + 2] = neighs[current_neigh][0];
        cyclic_neighsj[num_neighsj + 3] = neighs[current_neigh][1];

        // First term values needed
        double ar = cell_shape[current_neigh][0];
        double elasticity_area = 0.5 * (ar - 1);

        // Second term values needed
        double pe = cell_shape[current_neigh][1];
        double elasticity_peri = kp * (pe - p0);

        /*Now there are 2 vertices shared by cell i and j*/

        double nu[2][2] = {0.0};    //each row is a vertex and columns are nu_x and nu_y
        double nu_prev[2][2] = {0.0};
        double nu_next[2][2] = {0.0};

        int k = 2;
        while (cyclic_neighsj[k] != tag[i]) {
          k++;
          if (k == num_neighsj + 2) {
            error->one(FLERR,
                       " \n Inside force code: cell j did not find cell i as its neighbor \n");
          }
        }

        int cell_l = domain->closest_image(current_neigh, atom->map(cyclic_neighsj[k - 1]));
        int cell_k = domain->closest_image(current_neigh, atom->map(cyclic_neighsj[k + 1]));

        //error check
        int cell_l_from_i = domain->closest_image(i, atom->map(cyclic_neighsj[k - 1]));
        int cell_k_from_i = domain->closest_image(i, atom->map(cyclic_neighsj[k + 1]));

        if (cell_l != cell_l_from_i || cell_k != cell_k_from_i) {
          error->one(FLERR,
                     "Inside force code: nearest images of common neighbors for cell i and j do "
                     "not match");
        }

        /* for nu1 ---> (j,i,l): */
        nu[0][0] = (x[i][0] + x[current_neigh][0] + x[cell_l][0]) / 3.0;
        nu[0][1] = (x[i][1] + x[current_neigh][1] + x[cell_l][1]) / 3.0;

        //nu1_prev ---> (j,l,n)
        int cell_n = domain->closest_image(current_neigh, atom->map(cyclic_neighsj[k - 2]));
        nu_prev[0][0] = (x[current_neigh][0] + x[cell_l][0] + x[cell_n][0]) / 3.0;
        nu_prev[0][1] = (x[current_neigh][1] + x[cell_l][1] + x[cell_n][1]) / 3.0;

        //nu1_next ---> (j,i,k)
        nu_next[0][0] = (x[current_neigh][0] + x[i][0] + x[cell_k][0]) / 3.0;
        nu_next[0][1] = (x[current_neigh][1] + x[i][1] + x[cell_k][1]) / 3.0;

        /* for nu2 ---> (k,i,l) */
        nu[1][0] = (x[i][0] + x[current_neigh][0] + x[cell_k][0]) / 3.0;
        nu[1][1] = (x[i][1] + x[current_neigh][1] + x[cell_k][1]) / 3.0;

        //nu1_prev ---> (j,i,l)

        nu_prev[1][0] = (x[current_neigh][0] + x[i][0] + x[cell_l][0]) / 3.0;
        nu_prev[1][1] = (x[current_neigh][1] + x[i][1] + x[cell_l][1]) / 3.0;

        //nu1_next ---> (j,k,m)
        int cell_m = domain->closest_image(current_neigh, atom->map(cyclic_neighsj[k + 2]));
        nu_next[1][0] = (x[current_neigh][0] + x[cell_k][0] + x[cell_m][0]) / 3.0;
        nu_next[1][1] = (x[current_neigh][1] + x[cell_k][1] + x[cell_m][1]) / 3.0;

        double vertex_force_sum_t1[2] = {0.0};
        double vertex_force_sum_t2[2] = {0.0};

        //Now loop through both the vertices

        for (int n = 0; n < 2; n++) {

          //first term stuff
          double r_next_prev[3] = {nu_next[n][0] - nu_prev[n][0], nu_next[n][1] - nu_prev[n][1],
                                   0.0};
          double cp[3] = {0.0};
          double N[3] = {0, 0, 1};    // normal vector to the plane of cell layer (2D)
          getCP(cp, r_next_prev, N);

          //second term stuff
          double r_curr_prev[2] = {nu[n][0] - nu_prev[n][0], nu[n][1] - nu_prev[n][1]};
          double r_next_curr[2] = {nu_next[n][0] - nu[n][0], nu_next[n][1] - nu[n][1]};
          normalize(r_curr_prev);
          normalize(r_next_curr);
          double r_hatdiff_t2[2] = {r_curr_prev[0] - r_next_curr[0],
                                    r_curr_prev[1] - r_next_curr[1]};

          // Term 1 forces
          vertex_force_sum_t1[0] += cp[0] * Jac;
          vertex_force_sum_t1[1] += cp[1] * Jac;

          // Term 2 forces
          vertex_force_sum_t2[0] += r_hatdiff_t2[0] * Jac;
          vertex_force_sum_t2[1] += r_hatdiff_t2[1] * Jac;
        }

        F_t1[0] += elasticity_area * vertex_force_sum_t1[0];
        F_t1[1] += elasticity_area * vertex_force_sum_t1[1];
        F_t2[0] += elasticity_peri * vertex_force_sum_t2[0];
        F_t2[1] += elasticity_peri * vertex_force_sum_t2[1];
      }

      /*~~~~~~~~~~~~~~~~~ Force contribution from self ~~~~~~~~~~~~~~~~~*/

      double vertex_force_sum_t1[2] = {0.0};
      double vertex_force_sum_t2[2] = {0.0};

      // First term values needed
      double ai = cell_shape[i][0];
      double elasticity_area = 0.5 * (ai - 1);    //ka/2

      // Second term values needed
      double pi = cell_shape[i][1];
      double elasticity_peri = kp * (pi - p0);

      //Also find the minimum distance of cell from its bonds
      double dist_bond;
      double min_dist_bond = INFINITY;    //initialize to infinity

      int cell2, cell3;    //store the ids of neighs for which bond distance is minimum

      //Loop through the vertices/neighbors
      for (int n = 0; n < num_neighs[i]; n++) {
        double vert[2] = {0.0};
        double vert_next[2] = {0.0};
        double vert_prev[2] = {0.0};

        //current vertex
        vert[0] = vertices[i][2 * n];
        vert[1] = vertices[i][2 * n + 1];

        int neigh_j1 = domain->closest_image(i, atom->map(neighs[i][n]));    //local id of the neigh
        int neigh_j2;    //local id of next neighbor

        if (n == 0) {
          vert_next[0] = vertices[i][2 * (n + 1)];
          vert_next[1] = vertices[i][2 * (n + 1) + 1];
          vert_prev[0] = vertices[i][2 * (num_neighs[i] - 1)];
          vert_prev[1] = vertices[i][2 * (num_neighs[i] - 1) + 1];
          neigh_j2 = domain->closest_image(i, atom->map(neighs[i][1]));
        } else if (n == num_neighs[i] - 1) {
          vert_next[0] = vertices[i][0];
          vert_next[1] = vertices[i][1];
          vert_prev[0] = vertices[i][2 * (n - 1)];
          vert_prev[1] = vertices[i][2 * (n - 1) + 1];
          neigh_j2 = domain->closest_image(i, atom->map(neighs[i][0]));
        } else {
          vert_next[0] = vertices[i][2 * (n + 1)];
          vert_next[1] = vertices[i][2 * (n + 1) + 1];
          vert_prev[0] = vertices[i][2 * (n - 1)];
          vert_prev[1] = vertices[i][2 * (n - 1) + 1];
          neigh_j2 = domain->closest_image(i, atom->map(neighs[i][n + 1]));
        }

        //first term stuff
        double r_next_prev[3] = {vert_next[0] - vert_prev[0], vert_next[1] - vert_prev[1], 0.0};
        double cp[3] = {0.0};
        double N[3] = {0, 0, 1};    // normal vector to the plane of cell layer (2D)
        getCP(cp, r_next_prev, N);

        //second term stuff
        double r_curr_prev[2] = {vert[0] - vert_prev[0], vert[1] - vert_prev[1]};
        double r_next_curr[2] = {vert_next[0] - vert[0], vert_next[1] - vert[1]};
        normalize(r_curr_prev);
        normalize(r_next_curr);
        double rhatdiff_t2[2] = {r_curr_prev[0] - r_next_curr[0], r_curr_prev[1] - r_next_curr[1]};

        // Term 1 forces
        vertex_force_sum_t1[0] += cp[0] * Jac;
        vertex_force_sum_t1[1] += cp[1] * Jac;

        // Term 2 forces
        vertex_force_sum_t2[0] += rhatdiff_t2[0] * Jac;
        vertex_force_sum_t2[1] += rhatdiff_t2[1] * Jac;

        /*Viscosity stuff*/
        double xj1[2] = {x[neigh_j1][0], x[neigh_j1][1]};
        double xj2[2] = {x[neigh_j2][0], x[neigh_j2][1]};

        double A_line = xj1[1] - xj2[1];                      //A = y1-y2
        double B_line = xj2[0] - xj1[0];                      //B = x2-x1
        double C_line = xj1[0] * xj2[1] - xj2[0] * xj1[1];    //C = x1*y2 - x2*y1

        dist_bond = fabs(A_line * x[i][0] + B_line * x[i][1] + C_line) /
            sqrt(A_line * A_line + B_line * B_line);

        if (dist_bond < min_dist_bond) {
          min_dist_bond = dist_bond;    //update the minimum value
          cell2 = atom->map(
              tag[neigh_j1]);    //store local cell ids of the neighs forming bonds as well
          cell3 = atom->map(tag[neigh_j2]);
        }
      }

      //vertex model forces

      F_t1[0] += elasticity_area * vertex_force_sum_t1[0];
      F_t1[1] += elasticity_area * vertex_force_sum_t1[1];
      F_t2[0] += elasticity_peri * vertex_force_sum_t2[0];
      F_t2[1] += elasticity_peri * vertex_force_sum_t2[1];

      /*self force contri end*/

      /*Active force contribution*/

      double F_active[2] = {0.0};
      //random direction for active cell motion
      double n_ac_x = 2 * ((double) rand() / (RAND_MAX)) - 1;
      double n_ac_y = 2 * ((double) rand() / (RAND_MAX)) - 1;

      F_active[0] = fa * n_ac_x / sqrt(n_ac_x * n_ac_x + n_ac_y * n_ac_y);
      F_active[1] = fa * n_ac_y / sqrt(n_ac_x * n_ac_x + n_ac_y * n_ac_y);

      /*Update viscosities*/

      double nu_i;
      if (min_dist_bond >= c1)
        nu_i = nu0;
      else {
        nu_i = nu0 *
            (1 +
             pow(min_dist_bond - c1, 2.0) / pow(min_dist_bond, c2));    //eta in normal direction
      }

      nu[i] = max(nu_i, nu[i]);
      nu[cell2] = max(nu[cell2], nu_i);
      nu[cell3] = max(nu[cell3], nu_i);

      // Add all the force contributions
      fx[i] = -F_t1[0] - F_t2[0] + F_rep[0] + F_active[0];
      fy[i] = -F_t1[1] - F_t2[1] + F_rep[1] + F_active[1];
    }
  }

  // Update the force values (which become velocities now)

  for (int i = 0; i < nlocal; i++) {
    f[i][0] = fx[i] / nu[i];
    f[i][1] = fy[i] / nu[i];
    f[i][2] = 0.0;
  }

  /*********************************************************************** */
  /*CHECK: See if cells' motion results in self-intersecting polygons*/
  /************************************************************************ */

  double x_moved[nall][2];

  for (int i = 0; i < nlocal; i++) {
    x_moved[i][0] = x[i][0] + f[i][0] * dt;
    x_moved[i][1] = x[i][1] + f[i][1] * dt;
    int imj_atom = sametag[i];
    while (imj_atom != -1) {
      x_moved[imj_atom][0] = x[imj_atom][0] + f[i][0] * dt;
      x_moved[imj_atom][1] = x[imj_atom][1] + f[i][1] * dt;
      imj_atom = sametag[imj_atom];
    }
  }

  double vertices_moved[nlocal][2 * neighs_MAX];
  tagint neighs_moved[nlocal][neighs_MAX];

  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < neighs_MAX; j++) { neighs_moved[i][j] = neighs[i][j]; }
    //arrange the neighs cyclically based on x_moved
    arrange_cyclic_moved(neighs_moved[i], num_neighs[i], i, x_moved);
  }

  for (int i = 0; i < nlocal; i++) {
    get_cell_config_moved(vertices_moved[i], neighs_moved[i], num_neighs[i], i, x_moved);
  }

  //CHECK for self intersecting polygons

  for (int i = 0; i < nlocal; i++) {
    if (isSelfIntersecting(vertices_moved[i], 2 * num_neighs[i])) {
      //First freeze atom itself
      f[i][0] = 0.0;
      f[i][1] = 0.0;
      //Now freeze the neighbors
      for (int j = 0; j < num_neighs[i]; j++) {
        int neigh_j = atom->map(neighs[i][j]);
        f[neigh_j][0] = 0.0;
        f[neigh_j][1] = 0.0;
      }
    }
  }

  /******************************************************************************* */
  /*Outputting data for dump file*/
  /******************************************************************************* */

  // Reinitialize entries to 0
  if (peratom_flag) {
    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < size_peratom_cols; j++) { array_atom[i][j] = 0.0; }
    }
  }

  // Write cell shape data corresponding to current time step (if zeroth ts then for ts = 0)
  if (peratom_flag) {
    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < per_atom_out; j++) { array_atom[i][j] = cell_shape[i][j]; }
    }
  }

  // Write coordinate data
  if (peratom_flag) {
    for (int i = 0; i < nlocal; i++) {
      for (int j = per_atom_out; j < per_atom_out + 2; j++) { array_atom[i][j] = x[i][j - per_atom_out]; }
    }
  }

  /*~~~~~~~~~~~~~~~~~END of Post Force~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
}

/* <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*<<<<<<<<<<<<<<<<<<<<<< HELPER FUNCTIONS (BEGIN) >>>>>>>>>>>>>>>>>>>>>>>>>*/
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/*~~~~~~~~~~~~~~ FUNCTION 1: New arrange_cyclic for different position of atoms x_pos
                  (we need to use our own x_pos) ~~~~~~~~~~~~~~~~*/

void FixMorphoDynamic::arrange_cyclic(tagint *celli_neighs, int num_faces, int icell)
{
  double **x_pos = atom->x;
  // Temporary array to store angles
  double angles[num_faces];

  // Compute the angles for each neighbor
  for (int n = 0; n < num_faces; n++) {
    int cell_j = atom->map(celli_neighs[n]);          //local id of neighbor
    cell_j = domain->closest_image(icell, cell_j);    //nearest image of neighbor
    angles[n] = atan2(x_pos[cell_j][1] - x_pos[icell][1], x_pos[cell_j][0] - x_pos[icell][0]);
    if (isnan(angles[n])) { error->one(FLERR, "angle in sort neighbors function returned nan"); }
  }

  bool swapped;

  // Sort the neighbors based on their angles using a simple bubble sort

  for (int i = 0; i < num_faces - 1; i++) {
    swapped = false;
    for (int j = 0; j < num_faces - i - 1; j++) {
      if (angles[j] > angles[j + 1]) {
        std::swap(angles[j], angles[j + 1]);                // Swap angles
        std::swap(celli_neighs[j], celli_neighs[j + 1]);    // Swap corresponding neighs
        swapped = true;
      }
    }
    if (!swapped) { break; }
  }
}

/*~~~~~~FUNCTION 3 This function finds cell data entirely related to cell shape itself~~~~~~~*/
/*~~~~~~~~~~~~~~~~These include area, perimeter and shape tensor~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void FixMorphoDynamic::get_cell_data1(double *celli_vertices, double *celli_geom, tagint *celli_neighs,
                                   int num_faces, int icell)
{
  tagint *tag = atom->tag;
  double **x_pos = atom->x;

  //vertices array has to be clean everytime you rewrite it
  for (int j = 0; j < neighs_MAX * 2; j++) { celli_vertices[j] = 0.0; }

  // Coordinates of the vertices
  double vert[num_faces][2] = {0.0};

  // Find coordinates of each vertex
  for (int n = 0; n < num_faces; n++) {

    // Indices of current triangulation
    int mu1 = n;
    int mu2 = n + 1;

    // Wrap back to first vertex for final term
    if (mu2 == num_faces) { mu2 = 0; }

    // local id of mu1
    int j_mu1 = domain->closest_image(icell, atom->map(celli_neighs[mu1]));

    // local id of mu2
    int j_mu2 = domain->closest_image(icell, atom->map(celli_neighs[mu2]));

    // Coordinates of current triangulation
    double xn[3] = {x_pos[icell][0], x_pos[j_mu1][0], x_pos[j_mu2][0]};
    double yn[3] = {x_pos[icell][1], x_pos[j_mu1][1], x_pos[j_mu2][1]};

    // Find centroid
    vert[n][0] = (xn[0] + xn[1] + xn[2]) / 3.0;
    vert[n][1] = (yn[0] + yn[1] + yn[2]) / 3.0;
  }

  /*Force calculation depends on order of neighs and vertices must match*/
  /*Hence don't sort the vertices again!!!!!*/
  /*Also, if neighs are ordered cyclically, then vertices must also be ordered cyclically*/

  //sortverticesCCW(vert, num_faces, x_pos[icell][0], x_pos[icell][1]);

  double area = 0.0;
  double peri = 0.0;
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  double eig_vx, eig_vy;

  for (int n = 0; n < num_faces; n++) {

    celli_vertices[2 * n] = vert[n][0];
    celli_vertices[2 * n + 1] = vert[n][1];

    // Indices of current and next vertex
    int mu1 = n;
    int mu2 = n + 1;

    // Wrap back to first vertex for final term
    if (mu2 == num_faces) { mu2 = 0; }

    // Sum the area contribution
    area += 0.5 * (vert[mu1][0] * vert[mu2][1] - vert[mu1][1] * vert[mu2][0]);
    peri += sqrt(pow(vert[mu2][0] - vert[mu1][0], 2.0) + pow(vert[mu2][1] - vert[mu1][1], 2.0));

    // Find shape tensor
    A += pow(vert[mu2][0] - vert[mu1][0], 2.0);
    B += (vert[mu2][0] - vert[mu1][0]) * (vert[mu2][1] - vert[mu1][1]);
    C += pow(vert[mu2][1] - vert[mu1][1], 2.0);
  }

  double lam1 = 0.5 * (A + C + sqrt(pow(A - C, 2) + 4 * B * B));
  double lam2 = 0.5 * (A + C - sqrt(pow(A - C, 2) + 4 * B * B));
  double lam_max = max(lam1, lam2);
  double lam_min = min(lam1, lam2);

  if (B != 0) {
    eig_vx = B / sqrt(B * B + (lam_max - A) * (lam_max - A));
    eig_vy = (lam_max - A) / sqrt(B * B + (lam_max - A) * (lam_max - A));
  } else if (B == 0 && lam_max == A) {
    eig_vx = 1.0;
    eig_vy = 0.0;
  } else {
    eig_vx = 0.0;
    eig_vy = 1.0;
  }

  celli_geom[0] = area;
  celli_geom[1] = peri;
  celli_geom[2] = sqrt(lam_min / lam_max);
  celli_geom[3] = eig_vx;
  celli_geom[4] = eig_vy;
}

/*~~~~~~~~~~~~FUNCTION 4 This function finds cell data dependent on neighbors info~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~These include stress, energy and coordination~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void FixMorphoDynamic::get_cell_data2(int celli, double *celli_vertices, int num_faces,
                                   double **cells_geom, tagint *neighs_i)
{

  double vert[num_faces][2];

  double ai = cells_geom[celli][0];
  double pi = cells_geom[celli][1];
  double press_i = ai - 1;
  double eig_vx_i = cells_geom[celli][3];
  double eig_vy_i = cells_geom[celli][4];

  double T_xx = 0.0;
  double T_yy = 0.0;
  double T_xy = 0.0;

  double crdn = 0.0;    //coordination of cell i with its neighbors

  for (int n = 0; n < num_faces; n++) {

    // Contribution from neighbors

    int neighj = atom->map(neighs_i[n]);
    double eig_vx_j = cells_geom[neighj][3];
    double eig_vy_j = cells_geom[neighj][4];
    double cos_theta = eig_vx_i * eig_vx_j + eig_vy_i * eig_vy_j;
    crdn += cos_theta * cos_theta;
    double pj = cells_geom[neighj][1];

    int mui = n;
    int muj = n - 1;

    if (n == 0) { muj = num_faces - 1; }

    double Tij = kp * (pi + pj - 2 * p0);
    double lij_x = celli_vertices[2 * mui] - celli_vertices[2 * muj];
    double lij_y = celli_vertices[2 * mui + 1] - celli_vertices[2 * muj + 1];
    double lij = sqrt(lij_x * lij_x + lij_y * lij_y);
    T_xx += (Tij / lij) * lij_x * lij_x;
    T_yy += (Tij / lij) * lij_y * lij_y;
    T_xy += (Tij / lij) * lij_x * lij_y;
  }

  crdn = crdn / num_faces;

  cells_geom[celli][5] = 0.5 * pow((ai - 1), 2.0) + 0.5 * kp * pow((pi - p0), 2.0) -
      Js * pow(eig_vx_i, 2.0) - Jn * crdn;
  cells_geom[celli][6] = press_i + (0.5 / ai) * T_xx;
  cells_geom[celli][7] = press_i + (0.5 / ai) * T_yy;
  cells_geom[celli][8] = (0.5 / ai) * T_xy;
  cells_geom[celli][9] = crdn;
}

/*~~~~~~~~~~~FUNCTION: GET cell_data in the transition state~~~~~~~~*/

void FixMorphoDynamic::get_trans_state_geom(double *TS_geom_data, double vert[][2], int num_verts)
{
  double area = 0.0;
  double peri = 0.0;
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  double eig_vx, eig_vy;    //principle eigen vector

  //Find cell-shape data

  for (int n = 0; n < num_verts; n++) {

    // Indices of current and next vertex
    int mu1 = n;
    int mu2 = n + 1;

    // Wrap back to first vertex for final term
    if (mu2 == num_verts) { mu2 = 0; }

    // Sum the area contribution
    area += 0.5 * (vert[mu1][0] * vert[mu2][1] - vert[mu1][1] * vert[mu2][0]);
    peri += sqrt(pow(vert[mu2][0] - vert[mu1][0], 2.0) + pow(vert[mu2][1] - vert[mu1][1], 2.0));

    // Find shape tensor
    A += pow(vert[mu2][0] - vert[mu1][0], 2.0);
    B += (vert[mu2][0] - vert[mu1][0]) * (vert[mu2][1] - vert[mu1][1]);
    C += pow(vert[mu2][1] - vert[mu1][1], 2.0);
  }

  double lam1 = 0.5 * (A + C + sqrt(pow(A - C, 2) + 4 * B * B));
  double lam2 = 0.5 * (A + C - sqrt(pow(A - C, 2) + 4 * B * B));
  double lam_max = max(lam1, lam2);

  if (B != 0) {
    eig_vx = B / sqrt(B * B + (lam_max - A) * (lam_max - A));
    eig_vy = (lam_max - A) / sqrt(B * B + (lam_max - A) * (lam_max - A));
  } else if (B == 0 && lam_max == A) {
    eig_vx = 1.0;
    eig_vy = 0.0;
  } else {
    eig_vx = 0.0;
    eig_vy = 1.0;
  }

  TS_geom_data[0] = area;
  TS_geom_data[1] = peri;
  TS_geom_data[2] = eig_vx;
  TS_geom_data[3] = eig_vy;
}

/*~~~~~~~~~~~FUNCTION: get energy of cells in the transition state~~~~~~~~*/

double FixMorphoDynamic::get_trans_state_energy(int id, tagint *neighs_cell, int num_neighs,
                                             double TS_data[][4], double **cells_data,
                                             tagint celltags[4])
{
  double area_cell = TS_data[id][0];
  double peri_cell = TS_data[id][1];
  double eig_vx_cell = TS_data[id][2];
  double eig_vy_cell = TS_data[id][3];
  double crdn = 0.0;
  int neigh1, neigh2, neigh3;
  int id1, id2, id3;
  double eig_vx_neigh, eig_vy_neigh;

  if (id == 0) {
    neigh1 = atom->map(celltags[1]);
    neigh2 = atom->map(celltags[2]);
    neigh3 = atom->map(celltags[3]);
    id1 = 1;
    id2 = 2;
    id3 = 3;
  } else if (id == 1) {
    neigh1 = atom->map(celltags[0]);
    neigh2 = atom->map(celltags[2]);
    neigh3 = atom->map(celltags[3]);
    id1 = 0;
    id2 = 2;
    id3 = 3;
  } else if (id == 2) {
    neigh1 = atom->map(celltags[0]);
    neigh2 = atom->map(celltags[1]);
    neigh3 = atom->map(celltags[3]);
    id1 = 0;
    id2 = 1;
    id3 = 3;
  } else if (id == 3) {
    neigh1 = atom->map(celltags[0]);
    neigh2 = atom->map(celltags[1]);
    neigh3 = atom->map(celltags[2]);
    id1 = 0;
    id2 = 1;
    id3 = 2;
  }

  //Find crdn

  for (int n = 0; n < num_neighs; n++) {
    // Contribution from neighbors
    int neigh = atom->map(neighs_cell[n]);
    if (neigh == neigh1) {
      eig_vx_neigh = TS_data[id1][2];
      eig_vy_neigh = TS_data[id1][3];
    } else if (neigh == neigh2) {
      eig_vx_neigh = TS_data[id2][2];
      eig_vy_neigh = TS_data[id2][3];
    } else if (neigh == neigh3) {
      eig_vx_neigh = TS_data[id3][2];
      eig_vy_neigh = TS_data[id3][3];
    } else {
      eig_vx_neigh = cells_data[neigh][3];
      eig_vy_neigh = cells_data[neigh][4];
    }
    double cos_theta = eig_vx_cell * eig_vx_neigh + eig_vy_cell * eig_vy_neigh;
    crdn += cos_theta * cos_theta;
  }

  crdn = crdn / num_neighs;

  double energy = 0.5 * pow((area_cell - 1), 2.0) + 0.5 * kp * pow((peri_cell - p0), 2.0) -
      Js * pow(eig_vx_cell, 2.0) - Jn * crdn;
  return energy;
}

/*~~~~~~~~~~~~~~ FUNCTIONS: To test self intersection issue ~~~~~~~~~~~~~~~~*/

int FixMorphoDynamic::nearest_image(int i, int j, double x_pos[][2])
{

  if (j < 0) return j;

  int *sametag = atom->sametag;
  int nearest = j;
  double delx = x_pos[i][0] - x_pos[j][0];
  double dely = x_pos[i][1] - x_pos[j][1];
  double rsqmin = delx * delx + dely * dely;
  double rsq;

  while (sametag[j] >= 0) {
    j = sametag[j];
    delx = x_pos[i][0] - x_pos[j][0];
    dely = x_pos[i][1] - x_pos[j][1];
    rsq = delx * delx + dely * dely;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      nearest = j;
    }
  }

  return nearest;
}

void FixMorphoDynamic::arrange_cyclic_moved(tagint *celli_neighs, int num_faces, int icell,
                                         double x_pos[][2])
{
  // Temporary array to store angles
  double angles[num_faces];

  // Compute the angles for each neighbor
  for (int n = 0; n < num_faces; n++) {
    int cell_j = atom->map(celli_neighs[n]);         //local id of neighbor
    cell_j = nearest_image(icell, cell_j, x_pos);    //nearest image of neighbor
    angles[n] = atan2(x_pos[cell_j][1] - x_pos[icell][1], x_pos[cell_j][0] - x_pos[icell][0]);
    if (isnan(angles[n])) { error->one(FLERR, "angle in sort neighbors function returned nan"); }
  }

  bool swapped;

  // Sort the neighbors based on their angles using a simple bubble sort

  for (int i = 0; i < num_faces - 1; i++) {
    swapped = false;
    for (int j = 0; j < num_faces - i - 1; j++) {
      if (angles[j] > angles[j + 1]) {
        std::swap(angles[j], angles[j + 1]);                // Swap angles
        std::swap(celli_neighs[j], celli_neighs[j + 1]);    // Swap corresponding neighs
        swapped = true;
      }
    }
    if (!swapped) { break; }
  }
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void FixMorphoDynamic::get_cell_config_moved(double *celli_vertices, tagint *celli_neighs,
                                          int num_faces, int icell, double x_pos[][2])
{

  //vertices array has to be clean everytime you rewrite it
  for (int j = 0; j < neighs_MAX * 2; j++) { celli_vertices[j] = 0.0; }

  // Coordinates of the vertices
  double vert[num_faces][2] = {0.0};

  // Find coordinates of each vertex
  for (int n = 0; n < num_faces; n++) {

    // Indices of current triangulation
    int mu1 = n;
    int mu2 = n + 1;

    // Wrap back to first vertex for final term
    if (mu2 == num_faces) { mu2 = 0; }

    // local id of mu1
    int j_mu1 = nearest_image(icell, atom->map(celli_neighs[mu1]), x_pos);

    // local id of mu2
    int j_mu2 = nearest_image(icell, atom->map(celli_neighs[mu2]), x_pos);

    // Coordinates of current triangulation
    double xn[3] = {x_pos[icell][0], x_pos[j_mu1][0], x_pos[j_mu2][0]};
    double yn[3] = {x_pos[icell][1], x_pos[j_mu1][1], x_pos[j_mu2][1]};

    // Find centroid
    vert[n][0] = (xn[0] + xn[1] + xn[2]) / 3.0;
    vert[n][1] = (yn[0] + yn[1] + yn[2]) / 3.0;

    celli_vertices[2 * n] = vert[n][0];
    celli_vertices[2 * n + 1] = vert[n][1];
  }
}

/*~~~~~~~~~~~~~~ FUNCTION: Arrange cyclic of vertices ~~~~~~~~~~~~~~~~*/

void FixMorphoDynamic::sortverticesCCW(double vertices[][2], int num_faces, double xi, double yi)
{
  // Temporary array to store angles
  double angles[num_faces];

  // Compute the angles for each vertex
  for (int i = 0; i < num_faces; i++) {
    double x = vertices[i][0];
    double y = vertices[i][1];
    angles[i] = atan2(y - yi, x - xi);

    if (isnan(angles[i])) { error->one(FLERR, "angle in sort vertices function returned nan"); }
  }

  bool swapped;

  // Sort the vertices based on their angles using a simple bubble sort
  for (int i = 0; i < num_faces - 1; i++) {
    swapped = false;
    for (int j = 0; j < num_faces - i - 1; j++) {
      if (angles[j] > angles[j + 1]) {
        std::swap(angles[j], angles[j + 1]);    // Swap angles
        // Swap corresponding vertices
        std::swap(vertices[j][0], vertices[j + 1][0]);
        std::swap(vertices[j][1], vertices[j + 1][1]);
        swapped = true;
      }
    }
    if (!swapped) { break; }
  }
}

/*~~~~~~~~~~~~~~ FUNCTION 4: See if one list is cyclic permutation of other or not ~~~~~~~~~~~~~~~~*/

bool FixMorphoDynamic::is_cyclic_perm(tagint *original, tagint *rearranged, int num_of_neighs)
{

  tagint rearranged_concat[num_of_neighs * 2];

  //Fill in the rearranged_concat list

  for (int i = 0; i < num_of_neighs; i++) {
    rearranged_concat[i] = rearranged[i];
    rearranged_concat[i + num_of_neighs] = rearranged[i];
  }

  //find the location where first element of OG list is found
  int idx;

  for (int i = 0; i < num_of_neighs; i++) {
    if (rearranged_concat[i] == original[0]) {
      idx = i;
      break;
    }
  }

  // Now see if the original list is a sublist of concatenated list or not
  for (int i = 0; i < num_of_neighs; i++) {
    if (original[i] != rearranged_concat[idx + i]) { return false; }
  }

  return true;
}

/*~~~~~~~~~~~~~~~~~~~Function 3: Get common neighs for chosen bonded pair of atoms~~~~~~~~~~~~~~~~*/

int FixMorphoDynamic::get_comm_neigh(tagint *common_neighs, tagint *cella_neighs, tagint *cellb_neighs,
                                  int num_cella_neighs, int num_cellb_neighs)
{

  int num_comm = 0;

  for (int n = 0; n < num_cella_neighs; n++) {
    for (int m = 0; m < num_cellb_neighs; m++) {
      if (cella_neighs[n] == cellb_neighs[m]) {
        num_comm += 1;
        common_neighs[num_comm - 1] = cella_neighs[n];
      }
    }
  }

  return num_comm;
}

/*~~~~~~~~~~~~~~~~~~~ Function 4: see if cells k and l are already bonded or not ~~~~~~~~~~~~~~~~*/

int FixMorphoDynamic::unbonded(tagint *cell_neighs, int num_cell_neighs, tagint tag_neigh)
{
  for (int n = 0; n < num_cell_neighs; n++) {
    if (cell_neighs[n] == tag_neigh) { return 0; }
  }
  return 1;
}

/*~~~~~~~~~~~~~~~~~~~ Function 5: Add and remove neighbors ~~~~~~~~~~~~~~~~*/

void FixMorphoDynamic::remove_neigh(int celli, tagint tag_neigh, tagint *celli_neighs,
                                 int num_celli_neighs, int max_neighs)
{
  int pos_removed = -1;

  for (int n = 0; n < num_celli_neighs; n++) {
    if (celli_neighs[n] == tag_neigh) {
      pos_removed = n;
      break;
    }
  }

  //Now shift all the elements after pos_removed back
  for (int n = pos_removed; n < max_neighs; n++) {
    if (n < num_celli_neighs - 1) {
      celli_neighs[n] = celli_neighs[n + 1];
    } else {
      celli_neighs[n] = 0;
    }
  }
}

/*~~~~~~~~~~~~~~~~~~~ Function 6: Print neighs list as required ~~~~~~~~~~~~~~~~*/

void FixMorphoDynamic::print_neighs_list(tagint *cell_neighs, int num_cell_neighs, int cell)
{
  printf("Neighs List for local %d (global %d)--->", cell, atom->tag[cell]);
  for (int n = 0; n < num_cell_neighs; n++) { printf("%d,  ", cell_neighs[n]); }
  printf("\n");
}

/*~~~~~~~~~~~~~~~~~~~ Function 7: check if the quad is convex or concave ~~~~~~~~~~~~~~~~*/

bool FixMorphoDynamic::isConcave(double *p1, double *p2, double *p3, double *p4)
{
  // Compute cross products for each consecutive triplet of points
  double cross1 = crossProduct(p1, p2, p3);
  double cross2 = crossProduct(p2, p3, p4);
  double cross3 = crossProduct(p3, p4, p1);
  double cross4 = crossProduct(p4, p1, p2);

  // Check if all cross products have the same sign
  if ((cross1 > 0 && cross2 > 0 && cross3 > 0 && cross4 > 0) ||
      (cross1 < 0 && cross2 < 0 && cross3 < 0 && cross4 < 0)) {
    return false;    // Convex
  }
  return true;    // Concave
}

/*~~~~~~~~~~~~~~~~~~~ Function 8: cross product ~~~~~~~~~~~~~~~~*/

double FixMorphoDynamic::crossProduct(double *p1, double *p2, double *p3)
{
  double cross = (p2[0] - p1[0]) * (p3[1] - p2[1]) - (p3[0] - p2[0]) * (p2[1] - p1[1]);
  return cross;
}

/*~~~~~~~~~~~~~~~~~~~~~~FUnction 9: get cross product~~~~~~~~~~~~~~~~*/

void FixMorphoDynamic::getCP(double *cp, double *v1, double *v2)
{
  cp[0] = v1[1] * v2[2] - v1[2] * v2[1];
  cp[1] = v1[2] * v2[0] - v1[0] * v2[2];
  cp[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

/*~~~~~~~~~~~~~~~~~~~~~Function 10: normalize a vector~~~~~~~~~~~*/

// helper function: normalizes a vector

void FixMorphoDynamic::normalize(double *v)
{
  double norm = sqrt(pow(v[0], 2) + pow(v[1], 2));
  v[0] = v[0] / norm;
  v[1] = v[1] / norm;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS TO SEE IF POLYGON is SELF-INTERSECTING~~~~~~~~~~~~~~~~~*/

bool FixMorphoDynamic::onSegment(double x1, double y1, double x2, double y2, double x3, double y3)
{
  return x2 <= std::max(x1, x3) && x2 >= std::min(x1, x3) && y2 <= std::max(y1, y3) &&
      y2 >= std::min(y1, y3);
}

int FixMorphoDynamic::orientation(double x1, double y1, double x2, double y2, double x3, double y3)
{
  double val = (y2 - y1) * (x3 - x2) - (x2 - x1) * (y3 - y2);
  if (val == 0) return 0;
  return (val > 0) ? 1 : 2;
}

bool FixMorphoDynamic::doIntersect(double x1, double y1, double x2, double y2, double x3, double y3,
                                double x4, double y4)
{
  int o1 = orientation(x1, y1, x2, y2, x3, y3);
  int o2 = orientation(x1, y1, x2, y2, x4, y4);
  int o3 = orientation(x3, y3, x4, y4, x1, y1);
  int o4 = orientation(x3, y3, x4, y4, x2, y2);

  // General case
  if (o1 != o2 && o3 != o4) return true;

  // Special cases (collinear points lying on the segment)
  if (o1 == 0 && onSegment(x1, y1, x3, y3, x2, y2)) return true;
  if (o2 == 0 && onSegment(x1, y1, x4, y4, x2, y2)) return true;
  if (o3 == 0 && onSegment(x3, y3, x1, y1, x4, y4)) return true;
  if (o4 == 0 && onSegment(x3, y3, x2, y2, x4, y4)) return true;

  return false;
}

bool FixMorphoDynamic::isSelfIntersecting(double vertices[], int n)
{
  if (n < 8)
    return false;    // At least 4 points are required to form a polygon that can intersect itself

  for (int i = 0; i < n - 2; i += 2) {
    for (int j = i + 4; j < n; j += 2) {
      // Skip adjacent segments
      if ((i == 0 && j == n - 2) || (j == i + 2)) continue;

      if (doIntersect(vertices[i], vertices[i + 1], vertices[i + 2], vertices[i + 3], vertices[j],
                      vertices[j + 1], vertices[(j + 2) % n], vertices[(j + 3) % n])) {
        return true;
      }
    }
  }
  return false;
}

/*~~~~~~Function to deform the box and remap atom positions according to Fix berendsen~~~~~~~*/

void FixMorphoDynamic::resize_box()
{
  int i;
  double oldlo, oldhi, ctr;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // convert pertinent atoms and rigid bodies to lamda coords

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) domain->x2lamda(x[i], x[i]);

  for (auto &ifix : rfix) ifix->deform(0);

  // reset global and local box to new size/shape

  for (i = 0; i < 2; i++) {
    oldlo = domain->boxlo[i];
    oldhi = domain->boxhi[i];
    ctr = 0.5 * (oldlo + oldhi);
    domain->boxlo[i] = (oldlo - ctr) * dilation_MDN[i] + ctr;
    domain->boxhi[i] = (oldhi - ctr) * dilation_MDN[i] + ctr;
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) domain->lamda2x(x[i], x[i]);

  for (auto &ifix : rfix) ifix->deform(1);
}

/* <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*<<<<<<<<<<<<<<<<<<<<<< HELPER FUNCTIONS (END) >>>>>>>>>>>>>>>>>>>>>>>>>*/
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/* ---------------------------------------------------------------------- */

void FixMorphoDynamic::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMorphoDynamic::min_post_force(int vflag)
{
  post_force(vflag);
}

/*------------------------------------------------------------------------*/

int FixMorphoDynamic::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)

{
  int i, j, k, m;

  m = 0;

  if (commflag == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      for (k = 0; k < 4; k++) { buf[m++] = cell_shape[j][k]; }
    }
  } else if (commflag == 2) {
    int tmp1, tmp2;
    index = atom->find_custom("neighs", tmp1, tmp2);
    tagint **neighs = atom->iarray[index];
    for (i = 0; i < n; i++) {
      j = list[i];
      for (k = 0; k < neighs_MAX; k++) { buf[m++] = neighs[j][k]; }
    }
  }
  return m;
}

void FixMorphoDynamic::unpack_forward_comm(int n, int first, double *buf)
{
  int i, j, m, last;

  m = 0;
  last = first + n;

  if (commflag == 1) {
    for (i = first; i < last; i++) {
      for (j = 0; j < 4; j++) { cell_shape[i][j] = buf[m++]; }
    }
  } else if (commflag == 2) {
    int tmp1, tmp2;
    index = atom->find_custom("neighs", tmp1, tmp2);
    tagint **neighs = atom->iarray[index];
    for (i = first; i < last; i++) {
      for (j = 0; j < neighs_MAX; j++) { neighs[i][j] = buf[m++]; }
    }
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixMorphoDynamic::pack_exchange(int i, double *buf)
{
  int n = 0;
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++) buf[n++] = array_atom[i][m];
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixMorphoDynamic::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++) array_atom[nlocal][m] = buf[n++];
  }
  return n;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixMorphoDynamic::grow_arrays(int nmax)
{
  if (peratom_flag) {
    memory->grow(array_atom, nmax, size_peratom_cols, "fix_morphodynamic:array_atom");
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixMorphoDynamic::set_arrays(int i)
{
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++) array_atom[i][m] = 0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixMorphoDynamic::memory_usage()
{
  int maxatom = atom->nmax;
  double bytes = (double) maxatom * 4 * sizeof(double);
  return bytes;
}

/*--------------------------------------------------------------------
                        END OF MAIN CODE

---------------------------------------------------------------------*/

