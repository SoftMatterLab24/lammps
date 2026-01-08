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

#include "bond_bpm_gkv.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"
#include "table_file_reader.h"

#include <iostream>
#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondBPMGKV::BondBPMGKV(LAMMPS *_lmp) :
    BondBPM(_lmp), Ks(nullptr), Kj(nullptr), rcrit(nullptr), gamma(nullptr), zeta(nullptr),
    aT(nullptr), aT_temp(nullptr), id_fix_property_bond(nullptr)
{
  partial_flag = 1;
  smooth_flag = 1;
  normalize_flag = 0;
  temperature_flag = 0;
  writedata = 0;

  ntables = 0;
  tables = nullptr;

  nhistory = 3;
  update_flag = 1;
  id_fix_bond_history = utils::strdup("HISTORY_BPM_PRONY");

  single_extra = 5;
  svector = new double[5];

  nmax = 0;

  comm_forward = 1;
  comm_reverse = 1;

  dt_temp = 0;
}

/* ---------------------------------------------------------------------- */

BondBPMGKV::~BondBPMGKV()
{
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
    memory->destroy(Ks);
    memory->destroy(Kj);
    memory->destroy(rcrit);
    memory->destroy(zeta);
    memory->destroy(gamma);
    memory->destroy(aT);
    memory->destroy(aT_temp);
  }

}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double BondBPMGKV::store_bond(int n, int i, int j) // !! This is not updated for GKV !!
{
  int type;
  double delx, dely, delz, r;
  double k_temp, eta_temp, exp_j;
  double **x = atom->x;
  double dt = update->dt;
  double **bondstore = fix_bond_history->bondstore;
  tagint *tag = atom->tag;

  int **bond_type = atom->bond_type;

  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  r = sqrt(delx * delx + dely * dely + delz * delz);

  bondstore[n][0] = r;
  bondstore[n][1] = r;
  bondstore[n][2] = 0;

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) { 
        fix_bond_history->update_atom_value(i, m, 0, r); // rs
        fix_bond_history->update_atom_value(i, m, 1, r); // rn
        fix_bond_history->update_atom_value(i, m, 2, 0); // ep
        
      }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) { 
        fix_bond_history->update_atom_value(j, m, 0, r); //r0
        fix_bond_history->update_atom_value(j, m, 1, r); //rn
        fix_bond_history->update_atom_value(j, m, 2, 0); //ep

      }
    }
  }

  if (r < EPSILON) {
    error->one(FLERR, "Bond reference length too small");
  }

  return r;
}

/* ----------------------------------------------------------------------
  Store data for all bonds called once
------------------------------------------------------------------------- */

void BondBPMGKV::store_data()
{
  int i, j, n, m, type, N;
  double delx, dely, delz, r;
  double b, fn;
  double term1, eta;
  double **x = atom->x;
  double dt = update->dt;
  int **bond_type = atom->bond_type;

  double **bondstore = fix_bond_history->bondstore;

  for (i = 0; i < atom->nlocal; i++) {
    for (m = 0; m < atom->num_bond[i]; m++) {
      type = bond_type[i][m];

      //Skip if bond was turned off
      if (type < 0) continue;

      // map to find index n
      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");

      bond_lookup(type, atom->tag[i], atom->tag[j], N, b); // lookup N and b from tables

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];

      // Get closest image in case bonded with ghost
      domain->minimum_image(FLERR, delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);

      // compute initial bond force
      fn = 0; // for now assume zero initial force

      fix_bond_history->update_atom_value(i, m, 0, r); // rs
      fix_bond_history->update_atom_value(i, m, 1, r); // rn
      fix_bond_history->update_atom_value(i, m, 2, N); // N
      fix_bond_history->update_atom_value(i, m, 3, b); // b
      fix_bond_history->update_atom_value(i, m, 4, fn);// fn

      bondstore[m][0] = r;
      bondstore[m][1] = r;
      bondstore[m][2] = N;
      bondstore[m][3] = b;
      bondstore[m][4] = fn;

      const Table *tb = &tables[tabindex[type]];
      if (r < EPSILON) {
        error->one(FLERR, "Bond reference length too small");
      }

      // Loop through all Kelvin elements and initialize variable 
      for (int n = 0; n < N; n++ ) {

        dt_temp = dt;

        // Set internal viscous history variables to zero
        fix_bond_history->update_atom_value(i, m, n+5, 0);    // qi
        bondstore[m][n+5] = 0;

        // Set intial lengths of Kelvin-Voigt elements to zero
        fix_bond_history->update_atom_value(i, m, n+5+N, 0); // ri
        bondstore[m][n+5+N] = 0;

        // Compute viscosity and set
        term1 = M_PI*(n+1) / (2*N);
        eta = zeta[type] / (4.0*pow(sin(term1),2.0));
        fix_bond_history->update_atom_value(i, m, n+5+2*N, eta); // eta
        bondstore[m][n+5+2*N] = eta; // eta

      }
    }
  }
  fix_bond_history->post_neighbor();
}

/* ---------------------------------------------------------------------- */

void BondBPMGKV::compute(int eflag, int vflag)
{
  int i, bond_change_flag;

  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();
  }

  if (hybrid_flag) fix_bond_history->compress_history();

  int i1, i2, itmp, n, m, type, N;
  double delx, dely, delz, delvx, delvy, delvz;
  double e, ep, rsq, r, rs, rn , rjp , rinv, smooth, fs, fbond, dot;
  double b, eta, eta_temp, fn, rjn, qn, rjn1, qn1;
  double term1, term2, term3, numer, denom, lam;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double dt = update->dt;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double dim = domain->dimension;
  double invdim = 1.0 / dim;

  double **bondstore = fix_bond_history->bondstore;

  for (n = 0; n < nbondlist; n++) {

    // skip bond if already broken
    if (bondlist[n][2] <= 0) {
      continue;
    };

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    rs = bondstore[n][0]; 

    const Table *tb = &tables[tabindex[type]];

    // Update table (exponential constants)
    if (!(dt == dt_temp)) {
      update_table(type); // if the timestep has changed
    }
    
    if (!(aT[type] == aT_temp[type])) {
      update_table(type); // if the shift factor has changed
    }

    // Ensure pair is always ordered to ensure numerical operations
    // are identical to minimize the possibility that a bond straddling
    // an mpi grid (newton off) doesn't break on one proc but not the other 
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    // If bond hasn't been set - should be initialized to zero - (e.g. pour, fix bond/dynamic)
    //if (r0 < EPSILON || std::isnan(r0)) {
    //  r0 =store_bond(n, i1, i2);
    //}

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);    
    //e = (r0 !=0.0) ? (r - r0) / r0 : 0.0;

    //bond break criterion !! update
    if ((fabs(fbond/Ks[type]) > rcrit[type]) && break_flag) {  
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    // Check stability criterion
    int stable = 1;
    for (m = 0; m < tb->ninput; m++ ) {
      N = bondstore[n][2];
      b = bondstore[n][3];
      rjp  = bondstore[n][m+5+N]; 

      if (rjp > 0.90*b) {
        stable = 0;
        break;
      }
    }

    if (stable) {
      direct_solve(r, type, n, dt, fs);
    } else {
      iter_solve(r, type, n, dt, fs);
    }

    fbond = -fs;

    delvx = v[i1][0] - v[i2][0];
    delvy = v[i1][1] - v[i2][1];
    delvz = v[i1][2] - v[i2][2];
    dot = delx * delvx + dely * delvy + delz * delvz;
    fbond -= gamma[type] * dot * rinv;
    fbond *= rinv;

    if (smooth_flag) {
      // Disable for GKV
      smooth = 0;//(r - r0) / (r0 * rcrit[type]);
      smooth *= smooth;
      smooth *= smooth;
      smooth *= smooth;
      smooth = 1 - smooth;
      fbond *= smooth;
    }

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

    if (evflag) ev_tally(i1, i2, nlocal, newton_bond, 0.0, fbond, delx, dely, delz);
  }

  if (hybrid_flag) fix_bond_history->uncompress_history();
}

/* ---------------------------------------------------------------------- */

void BondBPMGKV::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(Ks, np1, "bond:Ks");
  memory->create(Kj, np1, "bond:Kj");
  memory->create(rcrit, np1, "bond:rcrit");
  memory->create(gamma, np1, "bond:gamma");
  memory->create(zeta, np1, "bond:zeta");
  memory->create(aT,np1,"bond:aT"); 
  memory->create(aT_temp,np1,"bond:aT_temp");
  memory->create(tabindex, np1, "bond:tabindex");
  memory->create(setflag, np1, "bond:setflag");

  for (int i = 1; i < np1; i++) setflag[i] = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMGKV::coeff(int narg, char **arg)
{
  if (narg < 7) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);
  
  double ks_zero = utils::numeric(FLERR, arg[1], false, lmp);
  double kj_zero = utils::numeric(FLERR, arg[2], false, lmp);
  double rcrit_one = utils::numeric(FLERR, arg[3], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[4], false, lmp);
  double zeta_one = utils::numeric(FLERR, arg[5], false, lmp);

  tables = (Table *) memory->srealloc(tables, (ntables + 1) * sizeof(Table), "bond:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (comm->me == 0) read_table(tb, arg[6], arg[7]);
  bcast_table(tb);

  // default values
  double aT_one = 1;

  // Parse optional remaining arguments
  int iarg = 8;
  while (iarg < narg) {
    if (temperature_flag) {
      if (iarg+1 > narg)  error->all(FLERR,"Incorrect args for bond coefficients");
      aT_one = utils::numeric(FLERR, arg[iarg], false, lmp);
      iarg += 1;
    } else error->all(FLERR,"Illegal fix bond/dynamic command");
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    Ks[i] = ks_zero;
    Kj[i] = kj_zero;
    rcrit[i] = rcrit_one;
    gamma[i] = gamma_one;
    zeta[i] = zeta_one;
    aT[i] = aT_one;
    aT_temp[i] = aT_one;
    setflag[i] = 1;
    tabindex[i] = ntables;
    
    count++;
  }
   ntables++;

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
  
}

/* ----------------------------------------------------------------------
   check for correct settings and create fix
------------------------------------------------------------------------- */

void BondBPMGKV::init_style()
{
  BondBPM::init_style();

  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Bond bpm/gkv requires ghost atoms store velocity");

}

/* ---------------------------------------------------------------------- */

void BondBPMGKV::settings(int narg, char **arg)
{
  nhistory = 3 * utils::numeric(FLERR, arg[0], false, lmp) + 5;
  
  BondBPM::settings(narg, arg);
  int iarg; 
  for (std::size_t i = 1; i < leftover_iarg.size(); i++) {
    iarg = leftover_iarg[i];
    if (strcmp(arg[iarg], "smooth") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for smooth");
      smooth_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else if (strcmp(arg[iarg], "normalize") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for normalize");
      normalize_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else if (strcmp(arg[iarg], "temp/shift") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for temp/shift");
      temperature_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else {
      error->all(FLERR, "Illegal bond bpm command, invalid argument {}", arg[iarg]);
    }
  }

  comm_forward = 1;
  comm_reverse = 1;

  if (smooth_flag && !break_flag)
    error->all(FLERR, "Illegal bond bpm command, must turn off smoothing with break no option");

}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMGKV::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&Ks[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Kj[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&rcrit[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gamma[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&zeta[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&aT[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&aT_temp[1], sizeof(double), atom->nbondtypes, fp);

  fwrite(&tabstyle, sizeof(int), 1, fp);
  fwrite(&tablength, sizeof(int), 1, fp);
  fwrite(&ntables, sizeof(int), 1, fp);
  fwrite(&nhistory, sizeof(int), 1, fp);

  for (int t = 0; t < ntables; t++) {
    Table &tb = tables[t];

    // write header: ninput and r0
    fwrite(&tb.ninput, sizeof(int), 1, fp);
    fwrite(&tb.r0, sizeof(double), 1, fp);

    fwrite(tb.iatomfile, sizeof(int), tb.ninput, fp);
    fwrite(tb.jatomfile, sizeof(int), tb.ninput, fp);
    fwrite(tb.Nfile, sizeof(int), tb.ninput, fp);
    fwrite(tb.bfile, sizeof(double), tb.ninput, fp);

  }
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMGKV::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  read_restart_settings(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &Ks[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Kj[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &rcrit[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gamma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &zeta[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &aT[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &aT_temp[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gamma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &aT[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &aT_temp[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);

    utils::sfread(FLERR, &tabstyle, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tablength, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &ntables, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &nhistory, sizeof(int), 1, fp, nullptr, error);
    
  }

  MPI_Bcast(&Ks[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Kj[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&rcrit[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&zeta[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&aT[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&aT_temp[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  MPI_Bcast(&tabstyle, 1, MPI_INT, 0, world);
  MPI_Bcast(&tablength, 1, MPI_INT, 0, world);
  MPI_Bcast(&ntables, 1, MPI_INT, 0, world);
  MPI_Bcast(&nhistory, 1, MPI_INT, 0, world);

  // allocate tables array on all procs
  tables = (Table *) memory->srealloc(tables, ntables * sizeof(Table), "bond:tables");

  // Read tables written by write_restart
  if (comm->me == 0) {
    utils::sfread(FLERR, &ntables, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&ntables, 1, MPI_INT, 0, world);

  // allocate tables array on all procs
  if (ntables > 0) {
    tables = (Table *) memory->srealloc(tables, ntables * sizeof(Table), "bond:tables");
  }

  for (int t = 0; t < ntables; t++) {
    Table *tb = &tables[t];
    null_table(tb);

    int ninput_local = 0;
    double r0_local = 0.0;

    // read header: ninput and r0
    if (comm->me == 0) {
      utils::sfread(FLERR, &ninput_local, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &r0_local, sizeof(double), 1, fp, nullptr, error);

      tb->iatomfile = nullptr; tb->jatomfile = nullptr; tb->Nfile = nullptr; tb->bfile = nullptr;
      memory->create(tb->iatomfile, tb->ninput, "bond:iatomfile");
      memory->create(tb->jatomfile, tb->ninput, "bond:jatomfile");
      memory->create(tb->Nfile, tb->ninput, "bond:Nfile");
      memory->create(tb->bfile, tb->ninput, "bond:bfile");

      utils::sfread(FLERR, tb->iatomfile, sizeof(int), tb->ninput, fp, nullptr, error);
      utils::sfread(FLERR, tb->jatomfile, sizeof(int), tb->ninput, fp, nullptr, error);
      utils::sfread(FLERR, tb->Nfile, sizeof(int), tb->ninput, fp, nullptr, error);
      utils::sfread(FLERR, tb->bfile, sizeof(double), tb->ninput, fp, nullptr, error);
    }
    
    bcast_table(tb);
  }

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPMGKV::write_restart_settings(FILE *fp)
{
  fwrite(&smooth_flag, sizeof(int), 1, fp);
  fwrite(&normalize_flag, sizeof(int), 1, fp);
  fwrite(&temperature_flag,sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMGKV::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &smooth_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &normalize_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &temperature_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&smooth_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&normalize_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&temperature_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMGKV::single(int type, double rsq, int i, int j, double &fforce)
{
  if (type <= 0) return 0.0;

  const Table *tb = &tables[tabindex[type]];
  double dt = update->dt;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  double **bondstore = fix_bond_history->bondstore;
  
  double r = sqrt(rsq);
  double rinv = 1.0 / r;

  int N;
  double rs, rn, rjp;
  double b, eta, eta_temp, fn, rjn, qn, rjn1, qn1;
  double term1, term2, term3, numer, denom, lam;

  // rn, ep, hn can be updated, so search bondlist vs. fix_bond_history->get_atom_value()
  tagint tagi = tag[i];
  tagint tagj = tag[j];
  tagint tag1, tag2;

  int n;
  for (n = 0; n < nbondlist; n++) {
    tag1 = tag[bondlist[n][0]];
    tag2 = tag[bondlist[n][1]];
    if ((tag1 == tagi && tag2 == tagj) || (tag1 == tagj && tag2 == tagi))
      break;
  }
   
  // retrieve bond history variables
  fn = bondstore[n][4];

  fforce = -fn;
  
  //double e = (r0 !=0.0) ? (r - r0) / r0 : 0.0;

  double **x = atom->x;
  double **v = atom->v;
  double delx = x[i][0] - x[j][0];
  double dely = x[i][1] - x[j][1];
  double delz = x[i][2] - x[j][2];
  double delvx = v[i][0] - v[j][0];
  double delvy = v[i][1] - v[j][1];
  double delvz = v[i][2] - v[j][2];
  double dot = delx * delvx + dely * delvy + delz * delvz;
  fforce -= gamma[type] * dot * rinv;
  fforce *= rinv;

  if (smooth_flag) {
    double smooth = 0;//(r0 != 0.0) ? (r - r0) / (r0 * ecrit[type]) : 0.0;
    smooth *= smooth;
    smooth *= smooth;
    smooth *= smooth;
    smooth = 1 - smooth;
    fforce *= smooth;
  }

  // set single_extra quantities

  svector[0] = 0;
  svector[1] = 0;
  svector[2] = 0;
  svector[3] = 0;
  svector[4] = 0;

  return 0.0;
}

/* ----------------------------------------------------------------------
    read from table file
 ------------------------------------------------------------------------- */

void BondBPMGKV::null_table(Table *tb)
{
  tb->Nfile = tb->bfile = nullptr;
  tb->N = tb->b =  nullptr;

}

/* ---------------------------------------------------------------------- */

void BondBPMGKV::free_table(Table *tb)
{
  memory->destroy(tb->iatomfile);
  memory->destroy(tb->jatomfile);
  memory->destroy(tb->Nfile);
  memory->destroy(tb->bfile);

  memory->destroy(tb->N);
  memory->destroy(tb->b);

}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void BondBPMGKV::read_table(Table *tb, char *file, char *keyword)
{
  double dt = update->dt;
  
  TableFileReader reader(lmp, file, "bond");

  char *line = reader.find_section_start(keyword);

  if (!line) error->one(FLERR, "Did not find keyword {} in table file", keyword);

  // read args on 2nd line of section
  // allocate table arrays for file values

  line = reader.next_line();
  param_extract(tb, line);
  memory->create(tb->iatomfile, tb->ninput, "bond:iatomfile");
  memory->create(tb->jatomfile, tb->ninput, "bond:jatomfile");
  memory->create(tb->Nfile, tb->ninput, "bond:Nfile");
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
      tb->iatomfile[i] = values.next_int(); 
      tb->jatomfile[i] = values.next_int();
      tb->Nfile[i] = values.next_int(); 
      tb->bfile[i] = values.next_double();
      //printf("reading | i: %i, N: %f, b: %f\n", i,tb->Nfile[i], tb->bfile[i]);
    } catch (TokenizerException &e) {
      error->one(FLERR, "Error parsing bond table '{}' line {} of {}. {}\nLine was: {}", keyword,
                 i + 1, tb->ninput, e.what(), line);
    }

  }

  printf("Read %i parameters from bond table\n",tb->ninput);

}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ r0
   N is required, other params are optional
------------------------------------------------------------------------- */

void BondBPMGKV::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->r0 = 0.0;

  try {
    ValueTokenizer values(line);

    while (values.has_next()) {
      std::string word = values.next_string();

      if (word == "N") {
        tb->ninput = values.next_int();
      } else {
        error->one(FLERR, "Unknown keyword {} in bond table parameters", word);
      }
    }
  } catch (TokenizerException &e) {
    error->one(FLERR, e.what());
  }

  if (tb->ninput == 0) error->one(FLERR, "Bond table parameters did not set N");
  //if (tb->ninput > nhistory - 3) error->one(FLERR, "New element exceeded elements per bond in table file");
  
}

/* ---------------------------------------------------------------------- */

 void BondBPMGKV::update_table(int type)
{   
  double dt = update->dt;

  //const Table *tb = &tables[tabindex[type]];
  // This is no longer needed 
  //  for (int m = 0; m < tb->ninput; m++ ) {

  //    k_temp = tb->kfile[m];
  //    eta_temp = aT[type] * tb->etafile[m]; 

  //    exp_j = exp(-dt * k_temp / eta_temp);
  //    tb->expfile[m] = exp_j;
  //  }

  dt_temp = dt;
  aT_temp[type] = aT[type]; 
}

/* ---------------------------------------------------------------------- */

void BondBPMGKV::bond_lookup(int type, int i, int j, int &N, double &b)
{
  int iatom, jatom;
  const Table *tb = &tables[tabindex[type]];

  for (int n = 0; n < tb->ninput; n++) {
    iatom = tb->iatomfile[n];
    jatom = tb->jatomfile[n];
    if ((iatom == i && jatom == j) || (iatom == j && jatom == i)) {
      N = tb->Nfile[n];
      b = tb->bfile[n];
      return; 
    }
  }
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile,ffile,fpflag,fplo,fphi,r0
------------------------------------------------------------------------- */

void BondBPMGKV::bcast_table(Table *tb) // *UPDATED
{
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);
  MPI_Bcast(&tb->r0, 1, MPI_DOUBLE, 0, world);

  int me;
  MPI_Comm_rank(world, &me);

  if (me > 0) {
    memory->create(tb->iatomfile, tb->ninput, "bond:iatomfile");
    memory->create(tb->jatomfile, tb->ninput, "bond:jatomfile");
    memory->create(tb->Nfile, tb->ninput, "bond:Nfile");
    memory->create(tb->bfile, tb->ninput, "bond:bfile");
  }

  MPI_Bcast(tb->iatomfile, tb->ninput, MPI_INT, 0, world);
  MPI_Bcast(tb->jatomfile, tb->ninput, MPI_INT, 0, world);
  MPI_Bcast(tb->Nfile, tb->ninput, MPI_INT, 0, world);
  MPI_Bcast(tb->bfile, tb->ninput, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   solvers for bond force given the total bond stretch
   if below stability limit, use fast direct solve
   else use slower but more accurate iterative solve
   
------------------------------------------------------------------------- */

void BondBPMGKV::direct_solve(double r, int type, int n, double dt , double &f)
{
  int    m, N;
  double rs, rn, rj_sum, b, eta, eta_temp, fn, rjn, qn, rjn1, qn1;
  double term1, term2, term3, numer, denom, lam;
  double fpred, fcor;

  double **bondstore = fix_bond_history->bondstore;

  // retrieve bond history variables
  rs = bondstore[n][0];
  rn = bondstore[n][1];
  N  = bondstore[n][2];
  b  = bondstore[n][3];
  fn = bondstore[n][4];

  double kj[N], exp_j[N], alph[N], qn_pred[N], rj_pred[N];
   
  // update bond length in bondstore
  bondstore[n][1] = r;
    
  term1 = 0.0; term2 = 0.0;
  for (int m = 0; m < N; m++ ) {

    // Get element specific params
    qn  = bondstore[n][m+5];        // old qi
    rjn = bondstore[n][m+5+N];      // old ri
    eta = bondstore[n][m+5+2*N];    // eta

    // update stiffness and exponential terms
    lam = rjn/b;

    numer = (pow(lam,2.0)- 3.0);
    denom = (pow(lam,2.0)- 1.0);
    kj[m] = Kj[type]*numer/denom/pow(b,2.0); // new stiffness
    eta_temp = aT[type] * eta;               // viscosity

    exp_j[m] = exp(-dt * kj[m] / eta_temp);  // exponential term

    if (dt/eta_temp < 1e-10){
        alph[m] = 1; // for small dt/eta take limit directly: alpha -> 1
    } else {
        alph[m] = (eta_temp / kj[m]) * (1 - exp_j[m]) / dt;
    }

    term1 = term1 + (qn*exp_j[m] - alph[m]*fn) / kj[m];
    term2 = term2 + (1 - alph[m]) / kj[m];
  }
 
  // Compute trial bond force
  fpred = (r + term1) / (1 / Ks[type] + term2);

  rj_sum = 0.0;
  for (m = 0; m < N; m++ ) {
      // update history variable
      qn1 = exp_j[m] * qn + (alph[m]) * (fpred - fn);
       
      rjn1 = (fpred - qn1) / kj[m];
        
      qn_pred[m] = qn1;
      rj_pred[m] = rjn1;
      rj_sum = rj_sum + rjn1; // total length of KV elements
  }

  // Corrector step
  rs = r - rj_sum;
  fcor = Ks[type] * rs;
 
  for (m = 0; m < N; m++ ) {
       
      rjn1 = (fcor - qn1) / kj[m];
      rjn1 = std::min(rjn1, rj_pred[m]);
      bondstore[n][m+5]   = qn1;  // qi
      bondstore[n][m+5+N] = rjn1; // ri

      lam = rjn1/b;
      numer = (pow(lam,2.0)- 3.0);
      denom = (pow(lam,2.0)- 1.0);
      kj[m] = Kj[type]*numer/denom/pow(b,2.0); // new stiffness
      
      qn1 = fcor - rjn1 * kj[m];

      bondstore[n][m+5]   = qn1;  // update qi
      bondstore[n][m+5+N] = rjn1; // update ri
  }

  f = fcor; 

  bondstore[n][0] = rs;    // update rs
  bondstore[n][4] = f;     // update fn

  return;
}

void BondBPMGKV::iter_solve(double r, int type, int n, double dt , double &f)
{
  int    m, N;
  int    iter_max, iter_local_max;
  double rs, rn, rjp, rj_sum, rj_low, rj_high, rj_mid;
  double b, lam, eta, eta_temp, qn1, qnp;
  double fn, f_low, f_high, f_mid;
  double R, G;
  double numer, denom;
  
  double **bondstore = fix_bond_history->bondstore;

  // retrieve bond history variables
  rs = bondstore[n][0];
  rn = bondstore[n][1];
  N  = bondstore[n][2];
  b  = bondstore[n][3];
  fn = bondstore[n][4];

  double kj[N], exp_j[N], alph[N], qn[N], rj_pred[N];
   
  // update bond length in bondstore
  bondstore[n][1] = r;

  // bracket stress
  f_low = 0.5 * fn;
  f_high = 1.1 * Ks[type] * rs;
  
  // Global bisection for bond force
  iter_max = 100;
  for (int iter = 0; iter < iter_max; iter++) {
    f_mid = 0.5 * (f_low + f_high);

    rj_sum = 0.0;
    // Solve each KV element length at this trial force
    for (m = 0; m < N; m++ ) {

      // Get element specific params
      //qnp  = bondstore[n][m+5];        // old qi
      rjp  = bondstore[n][m+5+N];      // old ri
      eta  = bondstore[n][m+5+2*N];    // eta

      // Local bisection for element length
      iter_local_max = 100;

      // bracket stretch
      rj_low = 0.95 * rjp;
      rj_high = b; // max stretch of element

      for (int iter_local = 0; iter_local < iter_local_max; iter_local++) {
        rj_mid = 0.5 * (rj_low + rj_high);

        // get stiffness
        lam = rj_mid/b;
        numer = (pow(lam,2.0)- 3.0);
        denom = (pow(lam,2.0)- 1.0);
        kj[m]    = Kj[type]*numer/denom/pow(b,2.0); // new stiffness

        eta_temp = aT[type] * eta;               // viscosity

        // residual
        G = rj_mid + dt/eta_temp * kj[m] * rj_mid - (rjp + dt/eta_temp*f_mid);

        if (G > 0) {
          rj_high = rj_mid;
        } else {
          rj_low = rj_mid;
        }

        // Convergence check
        if (fabs(G) < 1e-6) break;
      }

      rj_pred[m] = rj_mid; // predicted element length at mid force
      rj_sum = rj_sum + rj_mid;
    }

    // Compute residual
    R = f_mid/Ks[type] + rj_sum - r;

    if (R > 0) {
      f_high = f_mid;
    } else {
      f_low = f_mid;
    }

    // Convergence check
    if (fabs(R) < 1e-6) break;

  }
  
  // Update bond force and elastic spring length
  rs = r - rj_sum;
  f  = f_mid;

  bondstore[n][0] = rs;    // update rs
  bondstore[n][4] = f;     // update fn

  // Update bond history variables at converged force
  for (m = 0; m < N; m++ ) {

      // Get element specific params
      bondstore[n][m+5+N] = rj_pred[m];      // old ri
      
      // get stiffness
      lam = rj_pred[m]/b;
      numer = (pow(lam,2.0)- 3.0);
      denom = (pow(lam,2.0)- 1.0);
      kj[m] = Kj[type]*numer/denom/pow(b,2.0); // new stiffness

      // update history variable
      qn1 = f_mid - rj_pred[m] * kj[m];
      bondstore[n][m+5] = qn1;
  }

  return;
}


/* ---------------------------------------------------------------------- */

void *BondBPMGKV::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "aT") == 0) return (void *) aT;
  return nullptr;
}
