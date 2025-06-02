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
#include "table_file_reader.h"

#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondBPMProny::BondBPMProny(LAMMPS *_lmp) :
    BondBPM(_lmp), k0(nullptr), k1(nullptr), eta1(nullptr), ecrit(nullptr), gamma(nullptr),
    id_fix_property_bond(nullptr), vol_current(nullptr), dvol0(nullptr)
{
  partial_flag = 1;
  smooth_flag = 1;
  normalize_flag = 0;
  writedata = 0;

  ntables = 0;
  tables = nullptr;

  nhistory = 1;
  id_fix_bond_history = utils::strdup("HISTORY_BPM_PRONY");

  single_extra = 1;
  svector = new double[1];

  nmax = 0;

  comm_forward = 0;
  comm_reverse = 0;
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
    memory->destroy(k1);
    memory->destroy(eta1);
    memory->destroy(ecrit);
    memory->destroy(gamma);
    
  }

}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double BondBPMProny::store_bond(int n, int i, int j)
{
  double delx, dely, delz, r;
  double **x = atom->x;
  double **bondstore = fix_bond_history->bondstore;
  tagint *tag = atom->tag;

  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  r = sqrt(delx * delx + dely * dely + delz * delz);
  bondstore[n][0] = r;

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) { fix_bond_history->update_atom_value(i, m, 0, r); }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) { fix_bond_history->update_atom_value(j, m, 0, r); }
    }
  }

  return r;
}

/* ----------------------------------------------------------------------
  Store data for all bonds called once
------------------------------------------------------------------------- */

void BondBPMProny::store_data()
{
  int i, j, m, type;
  double delx, dely, delz, r;
  double **x = atom->x;
  int **bond_type = atom->bond_type;

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
      domain->minimum_image(delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);

      fix_bond_history->update_atom_value(i, m, 0, r);
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

  int i1, i2, itmp, n, type;
  double delx, dely, delz, delvx, delvy, delvz;
  double e, rsq, r, r0, rinv, smooth, fbond, dot;
  double tau1, term1, term2, h1, h_j;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double dt = update->dt
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double dim = domain->dimension;
  double invdim = 1.0 / dim;

  double **bondstore = fix_bond_history->bondstore;


  // First Maxell element
  tau1 = eta1[type] / k1[type];
  h1 = 0;

  term1 = exp(-1*dt / tau1) * h1;
  term2 = (1 - exp(-1*dt / tau1)) / (dt / tau1);

  // Loop through remaining Maxwell Elements
  for (m = 0; m < tb->ninput; m++) {

    h_j = 0;

    param_lookup(type, m, tau_j, eta_j);

    term1 += exp(-1*dt / tau_j) * h;
    term2 += gamma_j * (1 - exp(-1*dt / tau_j)) / (dt / tau_j);

  }

  for (n = 0; n < nbondlist; n++) {

    // skip bond if already broken
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    r0 = bondstore[n][0];

    // Ensure pair is always ordered to ensure numerical operations
    // are identical to minimize the possibility that a bond straddling
    // an mpi grid (newton off) doesn't break on one proc but not the other
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    // If bond hasn't been set - should be initialized to zero
    if (r0 < EPSILON || std::isnan(r0)) r0 = store_bond(n, i1, i2);

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);
    e = (r - r0) / r0;

    if ((fabs(e) > ecrit[type]) && break_flag) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);

      continue;
    }

    rinv = 1.0 / r;
    if (normalize_flag)
      fbond = -k0[type] * e;
    else
      fbond = k0[type] * (r0 - r);

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
  memory->create(k1, np1, "bond:k1");
  memory->create(eta1, np1, "bond:eta1");
  memory->create(ecrit, np1, "bond:ecrit");
  memory->create(gamma, np1, "bond:gamma");

  memory->create(tabindex, np1, "bond:tabindex");
  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMProny::coeff(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);
  
  double k_zero = utils::numeric(FLERR, arg[1], false, lmp);
  double k_one = utils::numeric(FLERR, arg[2], false, lmp);
  double eta_one = utils::numeric(FLERR, arg[3], false, lmp);
  double ecrit_one = utils::numeric(FLERR, arg[4], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[5], false, lmp);

  tables = (Table *) memory->srealloc(tables, (ntables + 1) * sizeof(Table), "bond:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (comm->me == 0) read_table(tb, arg[6], arg[7]);
  bcast_table(tb);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k0[i] = k_zero;
    k1[i] = k_one;
    ecrit[i] = ecrit_one;
    gamma[i] = gamma_one;
    setflag[i] = 1;
    count++;

    if (1.0 + ecrit[i] > max_stretch) max_stretch = 1.0 + ecrit[i];
  }

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

  if (smooth_flag && !break_flag)
    error->all(FLERR, "Illegal bond bpm command, must turn off smoothing with break no option");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMProny::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&k0[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&k1[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&eta1[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&ecrit[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gamma[1], sizeof(double), atom->nbondtypes, fp);

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
    utils::sfread(FLERR, &k0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &eta1[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &ecrit[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gamma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);

    utils::sfread(FLERR, &tabstyle, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tablength, sizeof(int), 1, fp, nullptr, error);
  }

  MPI_Bcast(&k0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k1[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&eta1[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ecrit[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

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
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMProny::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &smooth_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &normalize_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&smooth_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&normalize_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMProny::single(int type, double rsq, int i, int j, double &fforce)
{
  if (type <= 0) return 0.0;

  double r0;
  for (int n = 0; n < atom->num_bond[i]; n++) {
    if (atom->bond_atom[i][n] == atom->tag[j]) r0 = fix_bond_history->get_atom_value(i, n, 0);
  }

  double r = sqrt(rsq);
  double rinv = 1.0 / r;
  double e = (r - r0) / r0;

  if (normalize_flag)
    fforce = -k0[type] * e;
  else
    fforce = k0[type] * (r0 - r);

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
    double smooth = (r - r0) / (r0 * ecrit[type]);
    smooth *= smooth;
    smooth *= smooth;
    smooth *= smooth;
    smooth = 1 - smooth;
    fforce *= smooth;
  }

  // set single_extra quantities

  svector[0] = r0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

//int BondBPMProny::pack_reverse_comm(int n, int first, double *buf)
//{
//  int i, m, last;
//  m = 0;
//  last = first + n;
//  for (i = first; i < last; i++) buf[m++] = vol_current[i];
//  return m;
//}

/* ---------------------------------------------------------------------- */

//void BondBPMProny::unpack_reverse_comm(int n, int *list, double *buf)
//{
//  int i, j, m;
//  m = 0;
//  for (i = 0; i < n; i++) {
//    j = list[i];
//    vol_current[j] += buf[m++];
//  }
//}

/* ---------------------------------------------------------------------- */

//int BondBPMProny::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
//{
//  int i, j, m;
//  m = 0;
//  for (i = 0; i < n; i++) {
//    j = list[i];
//    buf[m++] = vol_current[j];
//  }
//  return m;
//}

/* ---------------------------------------------------------------------- */

//void BondBPMProny::unpack_forward_comm(int n, int first, double *buf)
//{
//  int i, m, last;
//  m = 0;
//  last = first + n;
//  for (i = first; i < last; i++) vol_current[i] = buf[m++];
//}

/* ----------------------------------------------------------------------
    read from table file
 ------------------------------------------------------------------------- */

void BondBPMProny::null_table(Table *tb) // *UPDATED
{
  tb->kfile = tb->etafile = nullptr;
  tb->k = tb->eta =  nullptr;

}

/* ---------------------------------------------------------------------- */

void BondBPMProny::free_table(Table *tb) // *UPDATED
{
  memory->destroy(tb->kfile);
  memory->destroy(tb->etafile);

  memory->destroy(tb->k);
  memory->destroy(tb->eta);

}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void BondBPMProny::read_table(Table *tb, char *file, char *keyword) // *UPDATED
{
  TableFileReader reader(lmp, file, "bond");

  char *line = reader.find_section_start(keyword);

  if (!line) error->one(FLERR, "Did not find keyword {} in table file", keyword);

  // read args on 2nd line of section
  // allocate table arrays for file values

  line = reader.next_line();
  param_extract(tb, line);
  memory->create(tb->kfile, tb->ninput, "bond:kfile");
  memory->create(tb->etafile, tb->ninput, "bond:etafile");

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
      tb->kfile[i] = values.next_double(); 
      tb->etafile[i] = values.next_double();
      
      if (tb->kfile[i] <= 0) error->one(FLERR, "Bond parameter must positive non-zero");

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

//double dt = (update->dt)*nevery;

 void BondBPMProny::param_lookup(int type, int ID, double &tau_j, double &gamma_j)
{
    double k_temp, eta_temp;
    const Table *tb = &tables[tabindex[type]];
    
    // Grab properties for (ID + 1)th Maxwell element (first element read in as coeff)
    k_temp = tb->kfile[ID];
    eta_temp = tb->etafile[ID];

    tau_j = eta_temp / k_temp;
    gamma_j = k_temp / k0[type];
  
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
  }

  MPI_Bcast(tb->kfile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->etafile, tb->ninput, MPI_DOUBLE, 0, world);

}

/* ---------------------------------------------------------------------- */
