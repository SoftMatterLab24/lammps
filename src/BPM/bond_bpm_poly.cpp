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

#include "bond_bpm_poly.h"

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

BondBPMPoly::BondBPMPoly(LAMMPS *_lmp) :
    BondBPM(_lmp), k0(nullptr), fcrit(nullptr), gamma(nullptr), 
    id_fix_property_bond(nullptr)
{
  partial_flag = 1;
  smooth_flag = 1;
  normalize_flag = 0;
  writedata = 0;

  ntables = 0;
  tables = nullptr;

  nhistory = 3;
  update_flag = 1;
  id_fix_bond_history = utils::strdup("HISTORY_BPM_POLY");

  //single_extra = 5;
  //svector = new double[5];

  nmax = 0;

  comm_forward = 1;
  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

BondBPMPoly::~BondBPMPoly()
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
    memory->destroy(fcrit);
    memory->destroy(gamma);
  }

}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double BondBPMPoly::store_bond(int n, int i, int j)
{
  int type;
  double delx, dely, delz, r;
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

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) { 
        fix_bond_history->update_atom_value(i, m, 0, r); // r0
      }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) { 
        fix_bond_history->update_atom_value(j, m, 0, r); //r0
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

void BondBPMPoly::store_data()
{
  int i, j, n, m, type;
  double delx, dely, delz, r;
  double N, b;
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
      domain->minimum_image(FLERR,delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);

      fix_bond_history->update_atom_value(i, m, 0, r);

      bondstore[m][0] = r;
      bondstore[m][1] = N;
      bondstore[m][2] = b;

      const Table *tb = &tables[tabindex[type]];
      if (r < EPSILON) {
        error->one(FLERR, "Bond reference length too small");
      }

    }
  }
  fix_bond_history->post_neighbor();
}

/* ---------------------------------------------------------------------- */

void BondBPMPoly::compute(int eflag, int vflag)
{
  int i, bond_change_flag;

  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();
  }

  if (hybrid_flag) fix_bond_history->compress_history();

  int i1, i2, itmp, n, m, type;
  double delx, dely, delz, delvx, delvy, delvz;
  double e, ep, rsq, r, r0, rn , rc , rinv, ebond, fbond, dot;
  double N, b, Nb, lam, numer, denom, term1, term2;

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

    //printf("Processing bond %d of %d\n", n+1, nbondlist);
    // skip bond if already broken
    if (bondlist[n][2] <= 0) {
      continue;
    };

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    r0 = bondstore[n][0]; 

    const Table *tb = &tables[tabindex[type]];

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
      error->one(FLERR, "This bond style does not support dynamic bond creation");
      r0 =store_bond(n, i1, i2);
    }

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);    
    e = (r0 !=0.0) ? (r - r0) / r0 : 0.0;

    // rate-independent part of bond force
    rinv = 1.0 / r;
    
    N = bondstore[n][1];
    b = bondstore[n][2];

    Nb = N * b;
    lam = r/Nb;

    // Calculate bond force

    numer = lam*(3.0 - pow(lam,2.0));
    denom = 1.0 - pow(lam,2.0);
    fbond = -k0[type]*numer/denom/b;

    // Calculate energy 

    term1 = pow(lam,2.0)/2.0;
    term2 = log(1.0 - pow(lam,2.0));
    ebond = N*(term1 - term2);

    //bond break criterion //disable for nonlinear
    if ((fabs(fbond) > fcrit[type]) && break_flag) {  
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    if (lam > 0.99) error->warning(FLERR, "POLY bond too long: {} {:.8}", update->ntimestep, lam);
    if (lam > 2.0) error->one(FLERR, "Bad PRONY bond"); 

    delvx = v[i1][0] - v[i2][0];
    delvy = v[i1][1] - v[i2][1];
    delvz = v[i1][2] - v[i2][2];
    dot = delx * delvx + dely * delvy + delz * delvz;
    fbond -= gamma[type] * dot * rinv;
    fbond *= rinv;

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

void BondBPMPoly::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(k0, np1, "bond:k0");
  memory->create(fcrit, np1, "bond:fcrit");
  memory->create(gamma, np1, "bond:gamma");
  memory->create(tabindex, np1, "bond:tabindex");
  memory->create(setflag, np1, "bond:setflag");

  for (int i = 1; i < np1; i++) setflag[i] = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMPoly::coeff(int narg, char **arg)
{
  if (!(narg >= 6)) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);
  
  double k_zero = utils::numeric(FLERR, arg[1], false, lmp);
  double fcrit_one = utils::numeric(FLERR, arg[2], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[3], false, lmp);

  tables = (Table *) memory->srealloc(tables, (ntables + 1) * sizeof(Table), "bond:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (comm->me == 0) read_table(tb, arg[4], arg[5]);
  bcast_table(tb);


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k0[i] = k_zero;
    fcrit[i] = fcrit_one;
    gamma[i] = gamma_one;
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

void BondBPMPoly::init_style()
{
  BondBPM::init_style();

  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Bond bpm/prony requires ghost atoms store velocity");

}

/* ---------------------------------------------------------------------- */

void BondBPMPoly::settings(int narg, char **arg)
{
  
  BondBPM::settings(narg, arg);

  int iarg; 
  for (std::size_t i = 0; i < leftover_iarg.size(); i++) {
    iarg = leftover_iarg[i];
    if (strcmp(arg[iarg], "smooth") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for smooth");
      smooth_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else if (strcmp(arg[iarg], "normalize") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for normalize");
      normalize_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
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

void BondBPMPoly::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&k0[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&fcrit[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gamma[1], sizeof(double), atom->nbondtypes, fp);
  
  fwrite(&tabstyle, sizeof(int), 1, fp);
  fwrite(&tablength, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMPoly::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  read_restart_settings(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &fcrit[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gamma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);

    utils::sfread(FLERR, &tabstyle, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tablength, sizeof(int), 1, fp, nullptr, error);
    
  }

  MPI_Bcast(&k0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&fcrit[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  MPI_Bcast(&tabstyle, 1, MPI_INT, 0, world);
  MPI_Bcast(&tablength, 1, MPI_INT, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPMPoly::write_restart_settings(FILE *fp)
{
  fwrite(&smooth_flag, sizeof(int), 1, fp);
  fwrite(&normalize_flag, sizeof(int), 1, fp);

}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMPoly::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &smooth_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &normalize_flag, sizeof(int), 1, fp, nullptr, error);
  }

  MPI_Bcast(&smooth_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&normalize_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMPoly::single(int type, double rsq, int i, int j, double &fforce)
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
   
  fforce = 0;
  
  double e = (r0 !=0.0) ? (r - r0) / r0 : 0.0;

  // New nonlinear - langevin
  N = bondstore[n][1]; // This needs to be from the bondstore !!!
  b = bondstore[n][2];

  Nb = N * b;
  lam = r/Nb;
  numer = lam*(3.0 - pow(lam,2.0));
  denom = 1.0 - pow(lam,2.0);

  fforce = -k0[type]*numer/denom/b;

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
    double smooth = (r0 != 0.0) ? (r - r0) / (r0 * fcrit[type]) : 0.0;
    smooth *= smooth;
    smooth *= smooth;
    smooth *= smooth;
    smooth = 1 - smooth;
    fforce *= smooth;
  }

  if (r = 0.0) {
    fforce = 0.0;
    return 0.0;
  }

  term1 = pow(lam,2.0)/2.0;
  term2 = log(1.0 - pow(lam,2.0));
  double eng = N*(term1 - term2);

  // set single_extra quantities

  //svector[0] = r0;
  //svector[1] = (1.0 + ep) * r0;
  //svector[2] = fel;
  //svector[3] = fint;
  //svector[4] = Hn;

  return eng;
}

/* ----------------------------------------------------------------------
    read from table file
 ------------------------------------------------------------------------- */

void BondBPMPoly::null_table(Table *tb)
{
  tb->Nfile = tb->bfile = nullptr;
  tb->N = tb->b =  nullptr;

}

/* ---------------------------------------------------------------------- */

void BondBPMPoly::free_table(Table *tb)
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

void BondBPMPoly::read_table(Table *tb, char *file, char *keyword)
{
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
      tb->Nfile[i] = values.next_double(); 
      tb->bfile[i] = values.next_double();
      //printf("reading | i: %i, N: %f, b: %f\n", i,tb->nfile[i], tb->bfile[i]);
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

void BondBPMPoly::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->r0 = 0.0;

  try {
    ValueTokenizer values(line);

    while (values.has_next()) {
      std::string word = values.next_string();

      if (word == "N") {
        tb->ninput = values.next_int();
        //printf("Num Bonds to Read: %i\n", tb->ninput);
      } else {
        error->one(FLERR, "Unknown keyword {} in bond table parameters", word);
      }
    }
  } catch (TokenizerException &e) {
    error->one(FLERR, e.what());
  }

  if (tb->ninput == 0) error->one(FLERR, "Bond table parameters did not set N");
}

/* ---------------------------------------------------------------------- */

void BondBPMPoly::bond_lookup(int type, int i, int j, double &N, double &b)
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

void BondBPMPoly::bcast_table(Table *tb) // *UPDATED
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
  MPI_Bcast(tb->Nfile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->bfile, tb->ninput, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */
