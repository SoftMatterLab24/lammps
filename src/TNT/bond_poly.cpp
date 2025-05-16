/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "bond_poly.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_bond_history.h"
#include "fix_store_local.h"
#include "fix_update_special_bonds.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "table_file_reader.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondPOLY::BondPOLY(LAMMPS *_lmp) : 
    Bond(_lmp), k(nullptr), fc(nullptr), id_fix_property_bond(nullptr), id_fix_dummy(nullptr), id_fix_dummy2(nullptr), id_fix_update(nullptr),
    id_fix_bond_history(nullptr), id_fix_store_local(nullptr), id_fix_prop_atom(nullptr),
    fix_store_local(nullptr), fix_bond_history(nullptr), fix_update_special_bonds(nullptr),
    pack_choice(nullptr), output_data(nullptr)
{

  break_flag = 1;
  overlay_flag = 0;
  ignore_special_flag = 1;
  prop_atom_flag = 0;
  nvalues = 0;

  writedata = 0;
  ntables = 0;
  tables = nullptr;

  nhistory = 3;
  update_flag = 0;
  id_fix_bond_history = utils::strdup("HISTORY_BPM_SPRING");

  single_extra = 1;
  svector = new double[1];

  nmax = 0;

  comm_forward = 0;
  comm_reverse = 0;

  // create dummy fix as placeholder for FixUpdateSpecialBonds & BondHistory
  // this is so final order of Modify:fix will conform to input script
  // BondHistory technically only needs this if updateflag = 1

  id_fix_dummy = utils::strdup(fmt::format("BPM_DUMMY_{}", instance_total));
  modify->add_fix(fmt::format("{} all DUMMY ", id_fix_dummy));

  id_fix_dummy2 = utils::strdup(fmt::format("BPM_DUMMY2_{}", instance_total));
  modify->add_fix(fmt::format("{} all DUMMY ", id_fix_dummy2));
}

BondPOLY::~BondPOLY()
{

  delete[] pack_choice;

  if (id_fix_dummy) modify->delete_fix(id_fix_dummy);
  if (id_fix_dummy2) modify->delete_fix(id_fix_dummy2);
  if (id_fix_update) modify->delete_fix(id_fix_update);
  if (fix_bond_history) modify->delete_fix(id_fix_bond_history);
  if (id_fix_store_local) modify->delete_fix(id_fix_store_local);
  if (id_fix_prop_atom) modify->delete_fix(id_fix_prop_atom);

  delete[] id_fix_dummy;
  delete[] id_fix_dummy2;
  delete[] id_fix_update;
  delete[] id_fix_bond_history;
  delete[] id_fix_store_local;
  delete[] id_fix_prop_atom;

  memory->destroy(output_data);

  delete[] svector;
  if (id_fix_property_bond && modify->nfix) {
    modify->delete_fix(id_fix_property_bond);
    delete[] id_fix_property_bond;
  }

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);
  
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(tabindex);
    memory->destroy(r0);
    memory->destroy(k);
    memory->destroy(fc);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void BondPOLY::coeff(int narg, char **arg) // *UPDATED
{
  if (narg != 5) error->all(FLERR, "Illegal bond_coeff command: must have 3 arguments");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double f_crit = utils::numeric(FLERR, arg[2], false, lmp);

  tables = (Table *) memory->srealloc(tables, (ntables + 1) * sizeof(Table), "bond:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (comm->me == 0) read_table(tb, arg[3], arg[4]);
  bcast_table(tb);

  // error check on table parameters

  //if (tb->ninput <= 1) error->all(FLERR, "Invalid bond table length: {}", tb->ninput);

  // store ptr to table in tabindex

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    fc[i] = f_crit;
    tabindex[i] = ntables;
    r0[i] = tb->r0;
    setflag[i] = 1;
    count++;
  }
  ntables++;

  if (count == 0) error->all(FLERR, "Illegal bond_coeff command");
}

/* ---------------------------------------------------------------------- */

void BondPOLY::init_style()
{
  if (id_fix_store_local) {
    auto ifix = modify->get_fix_by_id(id_fix_store_local);
    if (!ifix) error->all(FLERR, "Cannot find fix STORE/LOCAL id {}", id_fix_store_local);
    if (strcmp(ifix->style, "STORE/LOCAL") != 0)
      error->all(FLERR, "Incorrect fix style matched, not STORE/LOCAL: {}", ifix->style);
    fix_store_local = dynamic_cast<FixStoreLocal *>(ifix);
    fix_store_local->nvalues = nvalues;
  }

  if (!ignore_special_flag) {
    if (overlay_flag) {
      if (force->special_lj[1] != 1.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0 ||
          force->special_coul[1] != 1.0 || force->special_coul[2] != 1.0 || force->special_coul[3] != 1.0)
        error->all(FLERR,
                   "With overlay/pair yes, POLY bond styles require a value of 1.0 for all special_bonds weights");
      if (id_fix_update) {
        modify->delete_fix(id_fix_update);
        delete[] id_fix_update;
        id_fix_update = nullptr;
      }
    } else {
      // Require atoms know about all of their bonds and if they break
      if (force->newton_bond && break_flag)
        error->all(FLERR, "With overlay/pair no, or break yes, POLY bond styles require Newton bond off");

      // special lj must be 0 1 1 to censor pair forces between bonded particles
      if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
        error->all(FLERR,
                   "With overlay/pair no, POLY bond styles require special LJ weights = 0,1,1");
      // if bonds can break, special coulomb must be 1 1 1 to ensure all pairs are included in the
      //    neighbor list and 1-3 and 1-4 special bond lists are skipped
      if (break_flag && (force->special_coul[1] != 1.0 || force->special_coul[2] != 1.0 ||
          force->special_coul[3] != 1.0))
        error->all(FLERR,
                   "With overlay/pair no, and break yes, BPM bond styles requires special Coulomb weights = 1,1,1");

      if (id_fix_dummy && break_flag) {
        id_fix_update = utils::strdup("BPM_UPDATE_SPECIAL_BONDS");
        fix_update_special_bonds = dynamic_cast<FixUpdateSpecialBonds *>(modify->replace_fix(
            id_fix_dummy, fmt::format("{} all UPDATE_SPECIAL_BONDS", id_fix_update), 1));
        delete[] id_fix_dummy;
        id_fix_dummy = nullptr;
      }
    }

    // special 1-3 and 1-4 weights must be 1 to prevent building 1-3 and 1-4 special bond lists
    if (force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0 || force->special_coul[2] != 1.0 ||
        force->special_coul[3] != 1.0)
      error->all(FLERR, "Bond style bpm requires 1-3 and 1-4 special weights of 1.0");
  }

  if (force->angle || force->dihedral || force->improper)
    error->all(FLERR, "Bond style bpm cannot be used with 3,4-body interactions");
  if (atom->molecular == 2)
    error->all(FLERR, "Bond style bpm cannot be used with atom style template");

  // find all instances of bond history to delete/shift data
  // (bond hybrid may create multiple)
  histories = modify->get_fix_by_style("BOND_HISTORY");
  n_histories = histories.size();

  // If a bond type isn't set, must be using bond style hybrid
  hybrid_flag = 0;
  for (int i = 1; i <= atom->nbondtypes; i++)
    if (!setflag[i]) hybrid_flag = 1;
  fix_bond_history->setflag = setflag;

  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Bond bpm/spring requires ghost atoms store velocity");

  id_fix_property_bond = utils::strdup("BOND_BPM_SPRING_FIX_PROPERTY_ATOM");
  modify->add_fix(fmt::format("{} all property/atom d_vol d_vol0 ghost yes writedata no",id_fix_property_bond));
}

/* ----------------------------------------------------------------------
   global settings
   All args before store/local command are saved for potential args
     for specific bond BPM substyles
   All args after optional stode/local command are variables stored
     in the compute store/local
------------------------------------------------------------------------- */

void BondPOLY::settings(int narg, char **arg)
{
  leftover_iarg.clear();

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "store/local") == 0) {
      nvalues = 0;
      id_fix_store_local = utils::strdup(arg[iarg + 1]);
      store_local_freq = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      pack_choice = new FnPtrPack[narg - iarg - 1];
      iarg += 3;
      while (iarg < narg) {
        if (strcmp(arg[iarg], "id1") == 0) {
          pack_choice[nvalues++] = &BondPOLY::pack_id1;
        } else if (strcmp(arg[iarg], "id2") == 0) {
          pack_choice[nvalues++] = &BondPOLY::pack_id2;
        } else if (strcmp(arg[iarg], "time") == 0) {
          pack_choice[nvalues++] = &BondPOLY::pack_time;
        } else if (strcmp(arg[iarg], "x") == 0) {
          pack_choice[nvalues++] = &BondPOLY::pack_x;
        } else if (strcmp(arg[iarg], "y") == 0) {
          pack_choice[nvalues++] = &BondPOLY::pack_y;
        } else if (strcmp(arg[iarg], "z") == 0) {
          pack_choice[nvalues++] = &BondPOLY::pack_z;
        } else if (strcmp(arg[iarg], "x/ref") == 0) {
          pack_choice[nvalues++] = &BondPOLY::pack_x_ref;
          prop_atom_flag = 1;
        } else if (strcmp(arg[iarg], "y/ref") == 0) {
          pack_choice[nvalues++] = &BondPOLY::pack_y_ref;
          prop_atom_flag = 1;
        } else if (strcmp(arg[iarg], "z/ref") == 0) {
          pack_choice[nvalues++] = &BondPOLY::pack_z_ref;
          prop_atom_flag = 1;
        } else {
          break;
        }
        iarg++;
      }
    } else if (strcmp(arg[iarg], "break") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond poly command, missing option for break");
      break_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      iarg++;
    }
  }

  if (id_fix_store_local) {

    if (nvalues == 0)
      error->all(FLERR, "Storing local data must include at least one value to output");
    memory->create(output_data, nvalues, "bond/poly:output_data");

    auto ifix = modify->get_fix_by_id(id_fix_store_local);
    if (!ifix)
      ifix = modify->add_fix(
          fmt::format("{} all STORE/LOCAL {} {}", id_fix_store_local, store_local_freq, nvalues));
    fix_store_local = dynamic_cast<FixStoreLocal *>(ifix);

    // Use property/atom to save reference positions as it can transfer to ghost atoms
    // This won't work for instances where bonds are added (e.g. fix pour) but in those cases
    // a reference state isn't well defined
    if (prop_atom_flag == 1) {

      id_fix_prop_atom = utils::strdup("BPM_property_atom");
      char *x_ref_id = utils::strdup("BPM_X_REF");
      char *y_ref_id = utils::strdup("BPM_Y_REF");
      char *z_ref_id = utils::strdup("BPM_Z_REF");

      ifix = modify->get_fix_by_id(id_fix_prop_atom);
      if (!ifix)
        ifix = modify->add_fix(fmt::format("{} all property/atom d_{} d_{} d_{} ghost yes",
                                           id_fix_prop_atom, x_ref_id, y_ref_id, z_ref_id));

      int type_flag;
      int col_flag;
      index_x_ref = atom->find_custom(x_ref_id, type_flag, col_flag);
      index_y_ref = atom->find_custom(y_ref_id, type_flag, col_flag);
      index_z_ref = atom->find_custom(z_ref_id, type_flag, col_flag);

      delete[] x_ref_id;
      delete[] y_ref_id;
      delete[] z_ref_id;

      if (ifix->restart_reset) {
        ifix->restart_reset = 0;
      } else {
        double *x_ref = atom->dvector[index_x_ref];
        double *y_ref = atom->dvector[index_y_ref];
        double *z_ref = atom->dvector[index_z_ref];

        double **x = atom->x;
        for (int i = 0; i < atom->nlocal; i++) {
          x_ref[i] = x[i][0];
          y_ref[i] = x[i][1];
          z_ref[i] = x[i][2];
        }
      }
    }
  }

  // Set up necessary history fix
  if (!fix_bond_history) {
    fix_bond_history = dynamic_cast<FixBondHistory *>(modify->replace_fix(
        id_fix_dummy2, fmt::format("{} all BOND_HISTORY {} {}", id_fix_bond_history, update_flag, nhistory), 1));
    delete[] id_fix_dummy2;
    id_fix_dummy2 = nullptr;
  }

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(tabindex);
    memory->destroy(r0);
    memory->destroy(k);
    memory->destroy(fc);
  }
  allocated = 0;

  ntables = 0;
  tables = nullptr;

}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double BondPOLY::store_bond(int n, int i, int j)
{
  double delx, dely, delz, r;
  double N, b;
  int type;
  double **x = atom->x;
  int **bondlist = neighbor->bondlist;
  double **bondstore = fix_bond_history->bondstore;
  tagint *tag = atom->tag;

  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  type = bondlist[n][2];

  // Lookup to find bond information
  bond_lookup(type, n, N, b);

  r = sqrt(delx * delx + dely * dely + delz * delz);
  bondstore[n][0] = r;
  bondstore[n][1] = N;
  bondstore[n][2] = b;

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) { 
        //fix_bond_history->update_atom_value(i, m, 0, r); 
        //fix_bond_history->update_atom_value(i, m, 1, N);
        //fix_bond_history->update_atom_value(i, m, 2, b);
      }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) { 
        //fix_bond_history->update_atom_value(j, m, 0, r); 
        //fix_bond_history->update_atom_value(j, m, 1, N);
        //fix_bond_history->update_atom_value(j, m, 2, b);
        }
    }
  }

  return r, N, b;
}

/* ----------------------------------------------------------------------
  Store data for all bonds called once
------------------------------------------------------------------------- */

void BondPOLY::store_data()
{
  int i, j, m, type;
  double delx, dely, delz, r, N, b;
  double **x = atom->x;
  int **bond_type = atom->bond_type;

  double **bondstore = fix_bond_history->bondstore;

  for (i = 0; i < atom->nlocal; i++) {
    for (m = 0; m < atom->num_bond[i]; m++) {
      type = bond_type[i][m];

      //Skip if bond was turned off
      //if (type <= 0) continue;

      // map to find index n
      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];

      // Get closest image in case bonded with ghost
      domain->minimum_image(delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);

      bondstore[m][0] = r;

      // Lookup to find bond information
      bond_lookup(type, m, N, b);
      bondstore[m][1] = N;
      bondstore[m][2] = b;

      //printf("stored | atomI %i, atomJ %i, ID: %i,r: %f, N: %f, b: %f\n",atom->tag[i],atom->tag[j],m,r,N,b);

      fix_bond_history->update_atom_value(i, m, 0, r);
      fix_bond_history->update_atom_value(i, m, 1, N);
      fix_bond_history->update_atom_value(i, m, 2, b);
    }
  }

  fix_bond_history->post_neighbor();
}

/* ---------------------------------------------------------------------- */

void BondPOLY::compute(int eflag, int vflag) // *UPDATED
{

  // new from BPM
  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    //printf("Storing bond data\n");
    //store_data();  
  }

  int i1, i2, itmp, n, type, ID;
  double delx, dely, delz, ebond, fbond;
  double rsq, r, r0;
  double N, b, Nb, lam, numer, denom, term1, term2;

  ebond = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double **bondstore = fix_bond_history->bondstore;

  //printf("Atom IDs from bond 1: %i, %i,atom1: %i,atom2: %i\n", bondlist[0][0], bondlist[0][1],tag[0],tag[3]);
  //printf("Atom IDs from bond 2: %i, %i,atom1: %i,atom2: %i\n", bondlist[1][0], bondlist[1][1],tag[1],tag[2]);

  for (n = 0; n < nbondlist; n++) {
    
    //printf("Number of bonds: %i\n",nbondlist);
    // skip bond if already broken
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    r0 = bondstore[n][0];
    
    //printf("atom1: %i, atom2: %i\n",tag[i1],tag[i2]);

    // probably faster to get values from the lookup table
    //printf("Compute");
    bond_lookup(type, n, N, b);

    // Ensure pair is always ordered to ensure numerical operations
    // are identical to minimize the possibility that a bond straddling
    // an mpi grid (newton off) doesn't break on one proc but not the other
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    // If bond hasn't been set - should be initialized to zero
    if (r0 < EPSILON || std::isnan(r0)) {
      //r0 = store_bond(n, i1, i2);
    } 

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);

    //printf("C | ID: %i,r: %f, N: %f, b: %f\n",n,r,N,b);

    // Determine stretch ratio //
    Nb = N * b; 
    lam = sqrt(rsq)/Nb;

    // if lam -> 1, then chain is approaching contour length
        // issue a warning
    // if lam > 2 something serious is wrong, abort

    if (lam > 0.99) {
        error->warning(FLERR, "POLY bond too long: {} {:.8}", update->ntimestep, lam);
        if (lam > 2.0) error->one(FLERR, "Bad POLY bond");
        lam = 0.99;
    }

    // Calculate force magnitude

    numer = lam*(3.0 - pow(lam,2.0));
    denom = 1.0 - pow(lam,2.0);
    fbond = -k[type]*numer/denom/b;

    // Calculate energy 

    term1 = pow(lam,2.0)/2.0;
    term2 = log(1.0 - pow(lam,2.0));
    ebond = N*(term1 - term2);

    //printf("force: %f, critical force: %f\n", fbond, fc[type]);
    // See if bond breaks
    if (fabs(fbond) > fc[type] && break_flag) {
      bondlist[n][2] = -1;
      process_broken(i1, i2);
    }

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
        f[i1][0] += delx * fbond;
        f[i1][1] += dely * fbond;
        f[i1][2] += delz * fbond;
      }
  
      if (newton_bond || i2 < nlocal) {
        f[i2][0] -= delx * fbond;
        f[i2][1] -= dely * fbond;
        f[i2][2] -= delz * fbond;
      }
  
      if (evflag) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
  }

  if (hybrid_flag) fix_bond_history->uncompress_history();
}

/* ---------------------------------------------------------------------- */

void BondPOLY::allocate() // *UPDATED!
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(tabindex, np1, "bond:tabindex");
  memory->create(setflag, np1, "bond:setflag");
  memory->create(r0, np1, "bond:r0");
  memory->create(k, np1, "bond:k");
  memory->create(fc, np1, "bond:fc");

  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondPOLY::equilibrium_distance(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondPOLY::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondPOLY::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondPOLY::write_restart_settings(FILE *fp)
{
  fwrite(&tabstyle, sizeof(int), 1, fp);
  fwrite(&tablength, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondPOLY::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &tabstyle, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tablength, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&tabstyle, 1, MPI_INT, 0, world);
  MPI_Bcast(&tablength, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondPOLY::single(int type, double rsq, int i, int j, double &fforce) // *UPDATED
{
  if (type <= 0) {
    fforce = 0.0;
    return 0;
  } 

  double N, b, Nb;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  for (int n = 0; n < nbondlist; n++) {
        int atom1 = bondlist[n][0];
        int atom2 = bondlist[n][1];
        //printf("atom1: %i, atom2: %i, i: %i, j: %i\n",atom1,atom2,i,j);
        if ((atom1 == i && atom2 == j) || (atom1 == j && atom2 == i)) {
            //printf("Single");
            bond_lookup(type, n, N, b);
            //printf("retrieved | i: %i, j: %i, ID: %i, N: %f, b: %f\n",atom->tag[i],atom->tag[j],n,N,b);
            break;
        } else {
        }
  }
  
  Nb = N * b; 
  double r = sqrt(rsq);
  double lam = r/Nb;

  // if lam -> 1, then chain is approaching contour length
    // issue a warning
  // if lam > 2 something serious is wrong, abort

  if (lam > 0.99) {
    printf("Bond data | N: %f\n",N);
    error->warning(FLERR, "POLY bond too long: {} {:.8}", update->ntimestep, lam);
    if (lam > 2.0) error->one(FLERR, "Bad POLY bond");
    lam = 0.99;
  }

  double numer = lam*(3.0 - pow(lam,2.0));
  double denom = 1.0 - pow(lam,2.0);
  fforce = -(k[type]*numer/denom/b) / r ;

  if (r = 0.0) {
    fforce = 0.0;
    return 0.0;
  }

  //printf("N: %f, b: %f, r: %f, force: %f \n", N,b,r,fforce);

  double term1 = pow(lam,2.0)/2.0;
  double term2 = log(1.0 - pow(lam,2.0));
  double eng = N*(term1 - term2);

  return eng;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void BondPOLY::null_table(Table *tb) // *UPDATED
{
  tb->nfile = tb->bfile = nullptr;
  tb->n = tb->b =  nullptr;

}

/* ---------------------------------------------------------------------- */

void BondPOLY::free_table(Table *tb) // *UPDATED
{
  memory->destroy(tb->nfile);
  memory->destroy(tb->bfile);

  memory->destroy(tb->n);
  memory->destroy(tb->b);

}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void BondPOLY::read_table(Table *tb, char *file, char *keyword) // *UPDATED
{
  TableFileReader reader(lmp, file, "bond");

  char *line = reader.find_section_start(keyword);

  if (!line) error->one(FLERR, "Did not find keyword {} in table file", keyword);

  // read args on 2nd line of section
  // allocate table arrays for file values

  line = reader.next_line();
  param_extract(tb, line);
  memory->create(tb->nfile, tb->ninput, "bond:nfile");
  memory->create(tb->bfile, tb->ninput, "bond:bfile");

  // read n,b table values from file

  int r0idx = -1;

  reader.skip_line();
  for (int i = 0; i < tb->ninput; i++) {
    line = reader.next_line();
    if (!line)
      error->one(FLERR, "Data missing when parsing bond table '{}' line {} of {}.", keyword, i + 1,
                 tb->ninput);
    try {
      ValueTokenizer values(line);
      values.next_int();
      tb->nfile[i] = values.next_double(); 
      tb->bfile[i] = values.next_double();
      printf("reading | i: %i, N: %f, b: %f\n", i,tb->nfile[i], tb->bfile[i]);
    } catch (TokenizerException &e) {
      error->one(FLERR, "Error parsing bond table '{}' line {} of {}. {}\nLine was: {}", keyword,
                 i + 1, tb->ninput, e.what(), line);
    }

  }

}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ r0
   N is required, other params are optional
------------------------------------------------------------------------- */

void BondPOLY::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->r0 = 0.0;

  try {
    ValueTokenizer values(line);

    while (values.has_next()) {
      std::string word = values.next_string();

      if (word == "N") {
        tb->ninput = values.next_int();
        printf("Num Bonds to Read: %i\n", tb->ninput);
      } else {
        error->one(FLERR, "Unknown keyword {} in bond table parameters", word);
      }
    }
  } catch (TokenizerException &e) {
    error->one(FLERR, e.what());
  }

  if (tb->ninput == 0) error->one(FLERR, "Bond table parameters did not set N");
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile,ffile,fpflag,fplo,fphi,r0
------------------------------------------------------------------------- */

void BondPOLY::bcast_table(Table *tb) // *UPDATED
{
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);
  MPI_Bcast(&tb->r0, 1, MPI_DOUBLE, 0, world);

  int me;
  MPI_Comm_rank(world, &me);
  if (me > 0) {
    memory->create(tb->nfile, tb->ninput, "bond:nfile");
    memory->create(tb->bfile, tb->ninput, "bond:bfile");
  }

  MPI_Bcast(tb->nfile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->bfile, tb->ninput, MPI_DOUBLE, 0, world);

}

/* void BondPOLY::bond_lookup(int type, int ID, double &N, double &b) // *UPDATED
{
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  int atom1,atom2,id1,id2;

  const Table *tb = &tables[tabindex[type]];

  //printf("ID: %i\n", ID);
  // check if ID is in the table otherwise find minimum id
  if (ID < tb->ninput) {
    N = tb->nfile[ID];
    b = tb->bfile[ID];
  } else {
    
    atom1 = bondlist[ID][0];
    atom2 = bondlist[ID][1];

    for (size_t i = 0; i < tb->ninput; i++) {
      
      id1 = bondlist[i][0];
      id2 = bondlist[i][1];

      //printf("ID: %i, atom1: %i, atom2: %i, id1: %i, id2: %i\n", ID, tag[atom1], tag[atom2], tag[id1], tag[id2]);

      if ((tag[atom1] == tag[id1] && tag[atom2] == tag[id2]) || (tag[atom1] == tag[id2] && tag[atom2] == tag[id1])) {
        N = tb->nfile[i];
        b = tb->bfile[i];
        break;
      } else {
        //error->one(FLERR, "Bond not found in table");
      }
    }
  }

  //printf("Type: %i\n", type);
  //printf("N: %f\n", N);

  // add some error checking here
  //if (ID < 0 || ID >= tb->ninput) {
  //  error->one(FLERR, "Illegal bond in bond style poly: {} {}", type, ID);
  //}
  
} */

/* ---------------------------------------------------------------------- */

void BondPOLY::bond_lookup(int type, int ID, double &N, double &b) // *UPDATED
{
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  int atom1,atom2,id1,id2;

  const Table *tb = &tables[tabindex[type]];

    atom1 = bondlist[ID][0];
    atom2 = bondlist[ID][1];

    for (size_t i = 0; i < tb->ninput; i++) {
      
      id1 = bondlist[i][0];
      id2 = bondlist[i][1];

      //printf("ID: %i, atom1: %i, atom2: %i, id1: %i, id2: %i\n", ID, tag[atom1], tag[atom2], tag[id1], tag[id2]);

      if ((tag[atom1] == tag[id1] && tag[atom2] == tag[id2]) || (tag[atom1] == tag[id2] && tag[atom2] == tag[id1])) {
        N = tb->nfile[i];
        b = tb->bfile[i];
        break;
      } 
    }
  
}

/* ---------------------------------------------------------------------- */

void BondPOLY::process_broken(int i, int j)
{
  if (!break_flag)
    error->one(FLERR, "BPM bond broke with break no option");

  int nlocal = atom->nlocal;
  if (fix_store_local) {
    // If newton off, bond can break on two procs so only record if proc owns lower tag
    //    (BPM bond styles should sort so i -> atom with lower tag)
    if (force->newton_bond || (i < nlocal)) {
      for (int n = 0; n < nvalues; n++) (this->*pack_choice[n])(n, i, j);
      fix_store_local->add_data(output_data, i, j);
    }
  }

  if (fix_update_special_bonds) {
    // If this processor owns two copies of the bond (i.e. if the domain is periodic and 1 proc thick),
    //   skip instance where larger tag (j) owned
    int check = 1;
    if (i >= nlocal) {
      int imap = atom->map(atom->tag[i]);
      if (imap < nlocal) check = 0;
    }
    if (check) fix_update_special_bonds->add_broken_bond(i, j);
  }

  // Manually search and remove from atom arrays
  // need to remove in case special bonds arrays rebuilt

  int m, n;
  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;

  if (i < nlocal) {
    for (m = 0; m < num_bond[i]; m++) {
      if (bond_atom[i][m] == tag[j] && setflag[bond_type[i][m]]) {
        n = num_bond[i];
        bond_type[i][m] = bond_type[i][n - 1];
        bond_atom[i][m] = bond_atom[i][n - 1];
        for (auto &ihistory: histories) {
          auto fix_bond_history2 = dynamic_cast<FixBondHistory *>  (ihistory);
          fix_bond_history2->shift_history(i, m, n - 1);
          fix_bond_history2->delete_history(i, n - 1);
        }
        num_bond[i]--;
        break;
      }
    }
  }

  if (j < nlocal) {
    for (m = 0; m < num_bond[j]; m++) {
      if (bond_atom[j][m] == tag[i] && setflag[bond_type[j][m]]) {
        n = num_bond[j];
        bond_type[j][m] = bond_type[j][n - 1];
        bond_atom[j][m] = bond_atom[j][n - 1];
        for (auto &ihistory: histories) {
          auto fix_bond_history2 = dynamic_cast<FixBondHistory *>  (ihistory);
          fix_bond_history2->shift_history(j, m, n - 1);
          fix_bond_history2->delete_history(j, n - 1);
        }
        num_bond[j]--;
        break;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int BondPOLY::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++); //buf[m++] = vol_current[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void BondPOLY::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    //vol_current[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int BondPOLY::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    //buf[m++] = vol_current[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void BondPOLY::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++); //vol_current[i] = buf[m++];
}


/* ----------------------------------------------------------------------
   one method for every keyword bond bpm can output
   the atom property is packed into array or vector
------------------------------------------------------------------------- */

void BondPOLY::pack_id1(int n, int i, int /*j*/) //transfered
{
  tagint *tag = atom->tag;
  output_data[n] = tag[i];
}

/* ---------------------------------------------------------------------- */

void BondPOLY::pack_id2(int n, int /*i*/, int j) //transfered
{
  tagint *tag = atom->tag;
  output_data[n] = tag[j];
}

/* ---------------------------------------------------------------------- */

void BondPOLY::pack_time(int n, int /*i*/, int /*j*/) //transfered
{
  bigint time = update->ntimestep;
  output_data[n] = time;
}

/* ---------------------------------------------------------------------- */

void BondPOLY::pack_x(int n, int i, int j)
{
  double **x = atom->x;
  output_data[n] = (x[i][0] + x[j][0]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondPOLY::pack_y(int n, int i, int j)
{
  double **x = atom->x;
  output_data[n] = (x[i][1] + x[j][1]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondPOLY::pack_z(int n, int i, int j)
{
  double **x = atom->x;
  output_data[n] = (x[i][2] + x[j][2]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondPOLY::pack_x_ref(int n, int i, int j)
{
  double *x = atom->dvector[index_x_ref];
  output_data[n] = (x[i] + x[j]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondPOLY::pack_y_ref(int n, int i, int j)
{
  double *y = atom->dvector[index_y_ref];
  output_data[n] = (y[i] + y[j]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondPOLY::pack_z_ref(int n, int i, int j)
{
  double *z = atom->dvector[index_z_ref];
  output_data[n] = (z[i] + z[j]) * 0.5;
}