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

#include "bond_bpm_prony.h"

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

BondBPMProny::BondBPMProny(LAMMPS *_lmp) :
    BondBPM(_lmp), k0(nullptr), ecrit(nullptr), gamma(nullptr), lamc(nullptr), eplastic(nullptr), 
    aT(nullptr), aT_temp(nullptr), id_fix_property_bond(nullptr)
{
  partial_flag = 1;
  smooth_flag = 1;
  normalize_flag = 0;
  nonlinear_flag = 0;
  plastic_flag = 0;
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

BondBPMProny::~BondBPMProny()
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
    memory->destroy(k0);
    memory->destroy(ecrit);
    memory->destroy(gamma);
    memory->destroy(lamc);
    memory->destroy(eplastic);
    memory->destroy(aT);
    memory->destroy(aT_temp);
  }

}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double BondBPMProny::store_bond(int n, int i, int j)
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
        fix_bond_history->update_atom_value(i, m, 0, r); // r0
        fix_bond_history->update_atom_value(i, m, 1, r); // rn
        fix_bond_history->update_atom_value(i, m, 2, 0); // ep
        
        type = bond_type[i][m];
        const Table *tb = &tables[tabindex[type]];
        for (int l = 0; l < tb->ninput; l++ ) {

        // Compute exponential terms
        k_temp = tb->kfile[l];
        eta_temp = aT[type] * tb->etafile[l];

        exp_j = exp(-dt * k_temp / eta_temp);
        tb->expfile[l] = exp_j;
        dt_temp = dt;

        // Internal stress variable
        fix_bond_history->update_atom_value(i, m, l+3, 0);

        }  
      }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) { 
        fix_bond_history->update_atom_value(j, m, 0, r); //r0
        fix_bond_history->update_atom_value(j, m, 1, r); //rn
        fix_bond_history->update_atom_value(j, m, 2, 0); //ep

        type = bond_type[j][m];
        const Table *tb = &tables[tabindex[type]];
        for (int l = 0; l < tb->ninput; l++ ) {

        // Compute exponential terms
        k_temp = tb->kfile[l];
        eta_temp = aT[type] * tb->etafile[l];

        exp_j = exp(-dt * k_temp / eta_temp);
        tb->expfile[l] = exp_j;
        dt_temp = dt;

        // Internal stress variable
        fix_bond_history->update_atom_value(j, m, l+3, 0);     

        }
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

void BondBPMProny::store_data()
{
  int i, j, n, m, type;
  double delx, dely, delz, r;
  double k_temp, eta_temp, exp_j;
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

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];

      // Get closest image in case bonded with ghost
      domain->minimum_image(FLERR,delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);

      fix_bond_history->update_atom_value(i, m, 0, r);
      fix_bond_history->update_atom_value(i, m, 1, r);
      fix_bond_history->update_atom_value(i, m, 2, 0);

      bondstore[m][0] = r;
      bondstore[m][1] = r;
      bondstore[m][2] = 0;

      const Table *tb = &tables[tabindex[type]];
      if (r < EPSILON) {
        error->one(FLERR, "Bond reference length too small");
      }

      // Loop through all Maxwell elements and initialize variable 
      for (int n = 0; n < tb->ninput; n++ ) {

        // Compute exponential terms
        k_temp = tb->kfile[n];
        eta_temp = aT[type] * tb->etafile[n];
        
        exp_j = exp(-dt * k_temp / eta_temp);
        tb->expfile[n] = exp_j;
        dt_temp = dt;

        // Internal stress variable
        fix_bond_history->update_atom_value(i, m, n+3, 0);
        bondstore[m][n+3] = 0;

      }

    }
  }
  fix_bond_history->post_neighbor();
}

/* ---------------------------------------------------------------------- */

void BondBPMProny::compute(int eflag, int vflag)
{
  int i, bond_change_flag;

  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();
  }

  if (hybrid_flag) fix_bond_history->compress_history();

  int i1, i2, itmp, n, m, type;
  double delx, dely, delz, delvx, delvy, delvz;
  double e, ep, rsq, r, r0, rn , r0p, rc , rinv,  smooth, fbond, dot;
  double k_temp, eta_temp, exp_j, Hn, term1, term2, term3;
  double N, b, Nb, lam, numer, denom;

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
    r0 = bondstore[n][0]; 

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
    if (r0 < EPSILON || std::isnan(r0)) {
      r0 =store_bond(n, i1, i2);
    }

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);    
    e = (r0 !=0.0) ? (r - r0) / r0 : 0.0;

    rn = bondstore[n][1]; // This needs to be after bonds have been initialized
    ep = bondstore[n][2];
   
    // update bond length in bondstore
    bondstore[n][1] = r;
    
    //bond break criterion
    if ((fabs(e) > ecrit[type]) && break_flag) {  
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    //plastic calculations
    if (plastic_flag) {
      if (e > (ep + eplastic[type])) {
        ep = e - eplastic[type];
        bondstore[n][2] = ep;
      }

      if (e < (ep - eplastic[type])) {
        ep = e + eplastic[type];
        bondstore[n][2] = ep;
      }

      r0p = (1.0 + ep) * r0;
      bondstore[n][2] = ep;
    } else {
      r0p = r0;
    }

    // rate-independent part of bond force
    rinv = 1.0 / r;
    if (normalize_flag) {
      fbond = -k0[type] * (e - ep);
    } else if (nonlinear_flag) {
      //if (r > r0p) {
      //  rc = r0p * lamc[type]; // if bond is in tension
      //} else {
      //  rc = 0; // if bond is in compression
      // }
      //double lam = (r - r0p) / (rc - r0p);
      //fbond = -k0[type] * (r - r0p) / ( 2 * (1 - (lam * lam)));

      // New nonlinear - langevin
      N = tb->Nfile[n];
      b = tb->bfile[n];

      Nb = N * b;
      lam = r/Nb;
      numer = lam*(3.0 - pow(lam,2.0));
      denom = 1.0 - pow(lam,2.0);

      fbond = -k0[type]*numer/denom/b;

      printf("Nonlinear force %f %f %f\n",fbond,N,b);

    } else {
      fbond = k0[type] * (r0p - r);
    }

    // rate-dependent part of bond force
    // Loop through Maxwell elements
    for (m = 0; m < tb->ninput; m++ ) {

      //Get element specific params
      k_temp = tb->kfile[m];
      eta_temp = aT[type] * tb->etafile[m];
      exp_j = tb->expfile[m];

      // Get bond history variable
      Hn = bondstore[n][m+3];

      if (normalize_flag) {
        term1 = exp_j * Hn;
        term2 =  k_temp * ((rn - r) / r0) * (1 - exp_j) / (dt * k_temp / eta_temp);
      } else {
        term1 = exp_j * Hn;
        term2 =  k_temp * (rn - r) * (1 - exp_j) / (dt * k_temp / eta_temp);
      }
      fbond += (term1 + term2);
      
      // Update bond history variable
      Hn = term1 + term2;
      bondstore[n][m+3] = Hn;
    }

    delvx = v[i1][0] - v[i2][0];
    delvy = v[i1][1] - v[i2][1];
    delvz = v[i1][2] - v[i2][2];
    dot = delx * delvx + dely * delvy + delz * delvz;
    fbond -= gamma[type] * dot * rinv;
    fbond *= rinv;

    if (smooth_flag) {
      smooth = (r - r0) / (r0 * ecrit[type]);
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

void BondBPMProny::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(k0, np1, "bond:k0");
  memory->create(ecrit, np1, "bond:ecrit");
  memory->create(gamma, np1, "bond:gamma");
  memory->create(lamc,np1,"bond:lamc");
  memory->create(eplastic,np1,"bond:eplastic");
  memory->create(aT,np1,"bond:aT"); 
  memory->create(aT_temp,np1,"bond:aT_temp");
  memory->create(tabindex, np1, "bond:tabindex");
  memory->create(setflag, np1, "bond:setflag");

  for (int i = 1; i < np1; i++) setflag[i] = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMProny::coeff(int narg, char **arg)
{
  if (!(narg >= 8)) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);
  
  double k_zero = utils::numeric(FLERR, arg[1], false, lmp);
  double ecrit_one = utils::numeric(FLERR, arg[2], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[3], false, lmp);

  // grab file information
  char *file_nl = nullptr;
  char *keyword_nl = nullptr;
  if (nonlinear_flag) {
    try {
      file_nl = arg[9];
      keyword_nl = arg[10];
    } catch (...) {
      error->all(FLERR, "Incorrect bond coefficient for nonlinear option");
    }
  }

  tables = (Table *) memory->srealloc(tables, (ntables + 1) * sizeof(Table), "bond:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (comm->me == 0) read_table(tb, arg[4], arg[5],file_nl, keyword_nl);
  bcast_table(tb);

  double lamc_one = utils::numeric(FLERR, arg[6], false, lmp);
  double eplastic_one = utils::numeric(FLERR, arg[7], false, lmp);
  double aT_one = 1; 

  //if ((nonlinear_flag) && (lamc_one <= 1)) {
  //  error->all(FLERR, "Incorrect bond coefficient maximum extension must be greater than one");
  //}

  // Parse optional remaining arguments
  int iarg = 8;
  if (temperature_flag) {
    if (iarg+1 > narg)  error->all(FLERR,"Incorrect args for bond coefficients");
    aT_one = utils::numeric(FLERR, arg[8], false, lmp);
    iarg += 1;
  } 


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k0[i] = k_zero;
    ecrit[i] = ecrit_one;
    gamma[i] = gamma_one;
    lamc[i] = lamc_one;
    eplastic[i] = eplastic_one;
    aT[i] = aT_one;
    aT_temp[i] = aT_one;
    setflag[i] = 1;
    tabindex[i] = ntables;
    
    count++;

    if (1.0 + ecrit[i] > max_stretch) max_stretch = 1.0 + ecrit[i];
  }
   ntables++;

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
  
}

/* ----------------------------------------------------------------------
   check for correct settings and create fix
------------------------------------------------------------------------- */

void BondBPMProny::init_style()
{
  BondBPM::init_style();

  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Bond bpm/prony requires ghost atoms store velocity");

}

/* ---------------------------------------------------------------------- */

void BondBPMProny::settings(int narg, char **arg)
{
  nhistory = utils::numeric(FLERR, arg[0], false, lmp) + 3;
  
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
    }  else if (strcmp(arg[iarg], "plastic") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for plastic");
      plastic_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else if (strcmp(arg[iarg], "nonlinear") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for nonlinear");
      nonlinear_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
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

  if (nonlinear_flag && !break_flag)
    error->all(FLERR, "Illegal bond bpm command, must turn on breaking with nonlinear yes option");

  if (nonlinear_flag && plastic_flag)
    error->all(FLERR, "Illegal bond bpm command, must turn off plasticity with nonlinear yes option");

}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMProny::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&k0[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&ecrit[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gamma[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&lamc[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&eplastic[1], sizeof(double), atom->nbondtypes, fp);

  fwrite(&tabstyle, sizeof(int), 1, fp);
  fwrite(&tablength, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMProny::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  read_restart_settings(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &ecrit[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gamma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &lamc[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &eplastic[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &aT[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &aT_temp[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);

    utils::sfread(FLERR, &tabstyle, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tablength, sizeof(int), 1, fp, nullptr, error);
    
  }

  MPI_Bcast(&k0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ecrit[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&lamc[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&eplastic[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&aT[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&aT_temp[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  MPI_Bcast(&tabstyle, 1, MPI_INT, 0, world);
  MPI_Bcast(&tablength, 1, MPI_INT, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPMProny::write_restart_settings(FILE *fp)
{
  fwrite(&smooth_flag, sizeof(int), 1, fp);
  fwrite(&normalize_flag, sizeof(int), 1, fp);
  fwrite(&nonlinear_flag, sizeof(int), 1, fp);
  fwrite(&plastic_flag,sizeof(int), 1, fp);
  fwrite(&temperature_flag,sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMProny::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &smooth_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &normalize_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &nonlinear_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &plastic_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &temperature_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&smooth_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&normalize_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&nonlinear_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&plastic_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&temperature_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMProny::single(int type, double rsq, int i, int j, double &fforce)
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

  double r0, rn, r0p, rc, ep;
  double k_temp, eta_temp, exp_j, Hn, term1, term2;
  double fel, fint;
  double N, b, Nb, lam, numer, denom;

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
   
  r0 = bondstore[n][0];
  rn = bondstore[n][1];
  ep = bondstore[n][2];
   
  fforce = 0;
  // Loop through Maxwell elements (rate-dependent)
  for (int m = 0; m < tb->ninput; m++ ) {

    //Get element specific params
    k_temp = tb->kfile[m];
    eta_temp = aT[type] * tb->etafile[m];
    exp_j = tb->expfile[m];

    Hn = bondstore[n][m+3];

    if (normalize_flag) { 
      term1 = exp_j * Hn;
      term2 =  k_temp * ((rn - r) / r0) * (1 - exp_j) / (dt * k_temp / eta_temp);
    } else {
      term1 = exp_j * Hn;
      term2 =  k_temp * (rn - r) * (1 - exp_j) / (dt * k_temp / eta_temp);
    }

    fforce += (term1 + term2);
  }
  
  double e = (r0 !=0.0) ? (r - r0) / r0 : 0.0;

  //plastic calculations
  if (plastic_flag) {
    r0p = (1.0 + ep) * r0;
  } else
    r0p = r0;

  //rate-independent
  if (normalize_flag) {
    fel = -k0[type] * (e - ep);
    fforce += fel;
  } else if (nonlinear_flag) {
    //if (r > r0p) {
    //  rc = r0p * lamc[type]; // if bond is in tension
    //} else {
    //  rc = 0; // if bond is in compression
    //}
    //double lam = (r - r0p) / (rc - r0p);
    //fel = -k0[type] * (r - r0p) / ( 2 * (1 - (lam * lam)));
    //fforce += fel;

    // New nonlinear - langevin
    N = tb->Nfile[n];
    b = tb->bfile[n];

    Nb = N * b;
    lam = r/Nb;
    numer = lam*(3.0 - pow(lam,2.0));
    denom = 1.0 - pow(lam,2.0);

    fel = -k0[type]*numer/denom/b;
    fforce += fel;
  } else {
    fel = k0[type] * (r0p - r);
    fforce += fel;
  }

  fint = fforce - fel;

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
    double smooth = (r0 != 0.0) ? (r - r0) / (r0 * ecrit[type]) : 0.0;
    smooth *= smooth;
    smooth *= smooth;
    smooth *= smooth;
    smooth = 1 - smooth;
    fforce *= smooth;
  }

  // set single_extra quantities

  svector[0] = r0;
  svector[1] = (1.0 + ep) * r0;
  svector[2] = fel;
  svector[3] = fint;
  svector[4] = Hn;

  return 0.0;
}

/* ----------------------------------------------------------------------
    read from table file
 ------------------------------------------------------------------------- */

void BondBPMProny::null_table(Table *tb)
{
  tb->kfile = tb->etafile = tb->expfile = tb->Nfile = tb->bfile = nullptr;
  tb->k = tb->eta = tb->expj = tb->N = tb->b =  nullptr;

}

/* ---------------------------------------------------------------------- */

void BondBPMProny::free_table(Table *tb)
{
  memory->destroy(tb->kfile);
  memory->destroy(tb->etafile);
  memory->destroy(tb->expfile);
  memory->destroy(tb->Nfile);
  memory->destroy(tb->bfile);

  memory->destroy(tb->k);
  memory->destroy(tb->eta);
  memory->destroy(tb->expj);
  memory->destroy(tb->N);
  memory->destroy(tb->b);

}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void BondBPMProny::read_table(Table *tb, char *file, char *keyword, char *file_nl, char *keyword_nl)
{
  double dt = update->dt;
  
  TableFileReader reader(lmp, file, "bond");
  char *line = reader.find_section_start(keyword);

  if (!line) error->one(FLERR, "Did not find keyword {} in table file", keyword);

  // read args on 2nd line of section
  // allocate table arrays for file values

  line = reader.next_line();
  param_extract(tb, line);
  memory->create(tb->kfile, tb->ninput, "bond:kfile");
  memory->create(tb->etafile, tb->ninput, "bond:etafile");
  memory->create(tb->expfile, tb->ninput, "bond:expfile");
  
  // read table values from file

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
      tb->kfile[i] = values.next_double(); 
      tb->etafile[i] = values.next_double();
      tb->expfile[i] = 0;
      
      if (tb->kfile[i] <= 0) error->one(FLERR, "Bond parameter must positive non-zero");

    } catch (TokenizerException &e) {
      error->one(FLERR, "Error parsing bond table '{}' line {} of {}. {}\nLine was: {}", keyword,
                 i + 1, tb->ninput, e.what(), line);
    }

  }

  printf("Read %i parameters from bond table\n",tb->ninput);

  // if nonlinear option enabled, read additional table parameters
  if (nonlinear_flag) {

    TableFileReader reader(lmp, file_nl, "bond");
    char *line_nl = reader.find_section_start(keyword_nl);
    
    if (!line_nl) error->one(FLERR, "Did not find keyword {} in nonlinear table file", keyword_nl);
    // read table values from file
    // allocate memory for nonlinear

    line_nl = reader.next_line();
    nonlinear_param_extract(tb, line_nl);
    memory->create(tb->Nfile, tb->nninput, "bond:Nfile");
    memory->create(tb->bfile, tb->nninput, "bond:bfile");

    // need the new file name and keyword

    int r0idx = -1;

    reader.skip_line();
    for (int i = 0; i < tb->nninput; i++) {
      line_nl = reader.next_line();
      if (!line_nl)
      error->one(FLERR, "Data missing when parsing bond table '{}' line {} of {}.", keyword_nl, i + 1,
                 tb->nninput);
      try {
        ValueTokenizer values(line_nl);
        values.next_int();
        tb->Nfile[i] = values.next_double(); 
        tb->bfile[i] = values.next_double();

      } catch (TokenizerException &e) {
        error->one(FLERR, "Error parsing bond table '{}' line {} of {}. {}\nLine was: {}", keyword_nl,
                   i + 1, tb->nninput, e.what(), line_nl);
      }
    }

    printf("Read parameters for %i bonds from nonlinear bond table\n",tb->nninput);
  }

}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ r0
   N is required, other params are optional
------------------------------------------------------------------------- */

void BondBPMProny::param_extract(Table *tb, char *line)
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

  //if (!(tb->ninput == nhistory - 3)) error->one(FLERR, "Mismatched args for bond table parameter N");
  if (tb->ninput > nhistory - 3) error->one(FLERR, "New element exceeded elements per bond in table file");
  
}

/* ----------------------------------------------------------------------
   extract attributes for nonlinear parameter line in table section
   format of line: N value FP fplo fphi EQ r0
   N is required, other params are optional
------------------------------------------------------------------------- */

void BondBPMProny::nonlinear_param_extract(Table *tb, char *line_nl)
{
  tb->nninput = 0;

  try {
    ValueTokenizer values(line_nl);

    while (values.has_next()) {
      std::string word = values.next_string();

      if (word == "N") {
        tb->nninput = values.next_int();
      } else {
        error->one(FLERR, "Unknown keyword {} in bond table parameters", word);
      }
    }
  } catch (TokenizerException &e) {
    error->one(FLERR, e.what());
  }

  if (tb->nninput == 0) error->one(FLERR, "Bond table parameters did not set N");

  //if (!(tb->ninput == nhistory - 3)) error->one(FLERR, "Mismatched args for bond table parameter N");
  //if (tb->ninput > nhistory - 3) error->one(FLERR, "New element exceeded elements per bond in table file");
  
}

/* ---------------------------------------------------------------------- */

void BondBPMProny::update_table(int type)
{   
  double dt = update->dt;
  double k_temp, eta_temp, exp_j;
  const Table *tb = &tables[tabindex[type]];
    for (int m = 0; m < tb->ninput; m++ ) {

      k_temp = tb->kfile[m];
      eta_temp = aT[type] * tb->etafile[m]; 

      exp_j = exp(-dt * k_temp / eta_temp);
      tb->expfile[m] = exp_j;
    }
  dt_temp = dt;
  aT_temp[type] = aT[type]; 
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile,ffile,fpflag,fplo,fphi,r0
------------------------------------------------------------------------- */

void BondBPMProny::bcast_table(Table *tb) // *UPDATED
{
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);
  MPI_Bcast(&tb->r0, 1, MPI_DOUBLE, 0, world);

  int me;
  MPI_Comm_rank(world, &me);
  if (me > 0) {
    memory->create(tb->kfile, tb->ninput, "bond:kfile");
    memory->create(tb->etafile, tb->ninput, "bond:etafile");
    memory->create(tb->expfile, tb->ninput, "bond:expfile");
    memory->create(tb->Nfile, tb->nninput, "bond:Nfile");
    memory->create(tb->bfile, tb->nninput, "bond:bfile");
  }

  MPI_Bcast(tb->kfile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->etafile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->expfile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->Nfile, tb->nninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->bfile, tb->nninput, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

void *BondBPMProny::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "aT") == 0) return (void *) aT;
  return nullptr;
}
