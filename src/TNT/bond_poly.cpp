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
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include "table_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondPOLY::~BondPOLY()
{
    for (int m = 0; m < ntables; m++) free_table(&tables[m]);
    memory->sfree(tables);
  
    if (allocated) {
      memory->destroy(setflag);
      memory->destroy(r0);
      memory->destroy(tabindex);
    }
}

/* ---------------------------------------------------------------------- */

void BondPOLY::compute(int eflag, int vflag) // *UPDATED
{
    int i1, i2, n, type, ID;
    double delx, dely, delz, ebond, fbond;
    double rsq, r;
    double N, b, Nb, lam, numer, denom, term1, term2;

    ebond = 0.0;
    ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);
    ID = n;
    
    // need to get the correct bond info from bond table 
    bond_lookup(type, ID, N, b);

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
    fbond = -numer/denom/b;

    //printf("force: %f\n", fbond);

    // Calculate energy 

    term1 = pow(lam,2.0)/2.0;
    term2 = log(1.0 - pow(lam,2.0));
    ebond = N*(term1 - term2);

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
}

/* ---------------------------------------------------------------------- */

void BondPOLY::allocate() // *UPDATED
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(tabindex, np1, "bond:tabindex");
  memory->create(r0, np1, "bond:r0");
  //memory->create(b, np1, "bond:b");
  //memory->create(N, np1, "bond:N");
  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void BondPOLY::settings(int narg, char **arg) // *UPDATED
{
  if (narg != 0) error->all(FLERR, "Illegal bond_style command: must have 0 arguments");

  //tabstyle = NONE;
  //if (strcmp(arg[0], "linear") == 0)
  //  tabstyle = LINEAR;
  //else if (strcmp(arg[0], "spline") == 0)
  //  tabstyle = SPLINE;
  //else
  //  error->all(FLERR, "Unknown table style {} in bond style table", arg[0]);

  //tablength = utils::inumeric(FLERR, arg[1], false, lmp);
  //if (tablength < 2) error->all(FLERR, "Illegal number of bond table entries: {}", arg[1]);

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(tabindex);
  }
  allocated = 0;

  ntables = 0;
  tables = nullptr;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void BondPOLY::coeff(int narg, char **arg) // *UPDATED
{
  if (narg != 3) error->all(FLERR, "Illegal bond_coeff command: must have 3 arguments");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  tables = (Table *) memory->srealloc(tables, (ntables + 1) * sizeof(Table), "bond:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (comm->me == 0) read_table(tb, arg[1], arg[2]);
  bcast_table(tb);

  // error check on table parameters

  //if (tb->ninput <= 1) error->all(FLERR, "Invalid bond table length: {}", tb->ninput);

  //tb->lo = tb->rfile[0];
  //tb->hi = tb->rfile[tb->ninput - 1];
  //if (tb->lo >= tb->hi) error->all(FLERR, "Bond table values are not increasing");

  // spline read-in and compute r,e,f vectors within table

  //spline_table(tb);
  //compute_table(tb);

  // store ptr to table in tabindex

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    tabindex[i] = ntables;
    r0[i] = tb->r0;
    setflag[i] = 1;
    count++;
  }
  ntables++;

  if (count == 0) error->all(FLERR, "Illegal bond_coeff command");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
   should not be used, since don't know minimum of tabulated function
------------------------------------------------------------------------- */

double BondPOLY::equilibrium_distance(int i)
{
  return r0[i];
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

double BondPOLY::single(int type, double rsq, int /*i*/, int /*j*/, double &fforce) // *UPDATED
{
  double r = sqrt(rsq);
  int ID = 0; // ID is not used in this case, but needed for bond_lookup
  double N;
  double b;
  bond_lookup(type, ID, N, b);

  double Nb = N * b;
  double lam = sqrt(rsq)/Nb;

  // if lam -> 1, then chain is approaching contour length
    // issue a warning
  // if lam > 2 something serious is wrong, abort

  if (lam > 0.99) {
    error->warning(FLERR, "POLY bond too long: {} {:.8}", update->ntimestep, lam);
    if (lam > 2.0) error->one(FLERR, "Bad POLY bond");
    lam = 0.99;
  }

  double numer = lam*(3.0 - pow(lam,2.0));
  double denom = 1.0 - pow(lam,2.0);
  fforce = -(numer/denom/b) / r ;

  //printf("N: %f, b: %f, r: %f, force: %f \n", N,b,r,fforce);

  double term1 = pow(lam,2.0)/2.0;
  double term2 = log(1.0 - pow(lam,2.0));
  double eng = N*(term1 - term2);

  return eng;
}

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
  //memory->create(tb->ffile, tb->ninput, "bond:ffile");

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
      //tb->ffile[i] = values.next_double();
    } catch (TokenizerException &e) {
      error->one(FLERR, "Error parsing bond table '{}' line {} of {}. {}\nLine was: {}", keyword,
                 i + 1, tb->ninput, e.what(), line);
    }

    //if (tb->efile[i] < emin) {
    //  emin = tb->efile[i];
    //  r0idx = i;
    //}
  }

  // infer r0 from minimum of potential, if not given explicitly

  //if ((tb->r0 == 0.0) && (r0idx >= 0)) tb->r0 = tb->rfile[r0idx];

  // warn if force != dE/dr at any point that is not an inflection point
  // check via secant approximation to dE/dr
  // skip two end points since do not have surrounding secants
  // inflection point is where curvature changes sign

  //double r, e, f, rprev, rnext, eprev, enext, fleft, fright;

  //int ferror = 0;
  //for (int i = 1; i < tb->ninput - 1; i++) {
  //  n = tb->nfile[i];
  //  b = tb->bfile[i];
  //  rprev = tb->rfile[i - 1];
  //  rnext = tb->rfile[i + 1];
  //  e = tb->efile[i];
  // eprev = tb->efile[i - 1];
  //  enext = tb->efile[i + 1];
  //  f = tb->ffile[i];
  //  fleft = -(e - eprev) / (r - rprev);
  //  fright = -(enext - e) / (rnext - r);
  //  if (f < fleft && f < fright) ferror++;
  //  if (f > fleft && f > fright) ferror++;
  //}

  //if (ferror)
  //  error->warning(FLERR,
  //                 "{} of {} force values in table are inconsistent with -dE/dr.\n"
  //                 "WARNING:  Should only be flagged at inflection points",
  //                 ferror, tb->ninput);
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ r0
   N is required, other params are optional
------------------------------------------------------------------------- */

void BondPOLY::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->fpflag = 0;
  tb->r0 = 0.0;

  try {
    ValueTokenizer values(line);

    while (values.has_next()) {
      std::string word = values.next_string();

      if (word == "N") {
        tb->ninput = values.next_int();
      //} else if (word == "FP") {
      //  tb->fpflag = 1;
      //  tb->fplo = values.next_double();
      //  tb->fphi = values.next_double();
      //} else if (word == "EQ") {
      //  tb->r0 = values.next_double();
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
    //memory->create(tb->ffile, tb->ninput, "bond:ffile");
  }

  MPI_Bcast(tb->nfile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->bfile, tb->ninput, MPI_DOUBLE, 0, world);
  //MPI_Bcast(tb->ffile, tb->ninput, MPI_DOUBLE, 0, world);

  MPI_Bcast(&tb->fpflag, 1, MPI_INT, 0, world);
  //if (tb->fpflag) {
  //  MPI_Bcast(&tb->fplo, 1, MPI_DOUBLE, 0, world);
  //  MPI_Bcast(&tb->fphi, 1, MPI_DOUBLE, 0, world);
  //}
}

void BondPOLY::bond_lookup(int type, int ID, double &N, double &b) // *UPDATED
{
  //if (!std::isfinite(x)) { error->one(FLERR, "Illegal bond in bond style table"); }

  const Table *tb = &tables[tabindex[type]];

  // add some error checking here
  if (ID < 0 || ID >= tb->ninput) {
    //error->one(FLERR, "Illegal bond in bond style poly: {} {}", type, ID);
  }
 
  //printf("ID: %i\n", ID);
  //printf("Type: %i\n", type);
  //printf("N: %f\n", tb->nfile[0]);

  N = tb->nfile[ID];
  b = tb->bfile[ID];
  
}
