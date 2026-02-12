// clang-format off
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

#include "fix_bond_rupture.h"
#include "fix_bond_history.h"

#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "respa.h"
#include "update.h"
#include "table_file_reader.h"

#include <iostream>
#include <cstring>
#include <utility>
#include "math_const.h"
#include "random_mars.h"

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBondRupture::FixBondRupture(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  random(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal fix bond/rupture command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  dynamic_group_allow = 1;
  force_reneighbor = 1;
  next_reneighbor = -1;

  ntables = 0;
  tables = nullptr;

  // Default settings
  seed = 12345;

  // style_flags:
  flag_dist = 0; flag_fraction = 0; flag_slip = 0; flag_slip_catch = 0; flag_rate = 0;
  
  // keyword flags:
  flag_table = 0; flag_distibution = 0; flag_crit = 0;

  btype = utils::inumeric(FLERR,arg[3],false,lmp);
  const char *style_name = arg[4];

  if (btype < 1 || btype > atom->nbondtypes)
  error->all(FLERR,"Invalid bond type in fix bond/rupture command");

  int iarg = 5;
  // Parse style-specific values
    if (strcmp(style_name,"dist") == 0) {
        if (iarg+1 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_dist = 1; 
        double rcrit = utils::numeric(FLERR,arg[iarg],false,lmp);
        if (rcrit < 0.0) error->all(FLERR,"Illegal fix bond/rupture command");
        rcritsq = rcrit*rcrit;
        ndata = 1;
        iarg += 1;
    } else if (strcmp(style_name,"prob/fraction") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_fraction = 1; 
        p_fraction = utils::numeric(FLERR,arg[iarg],false,lmp);
        seed = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (p_fraction < 0.0 || p_fraction > 1.0) error->all(FLERR,"Illegal fix bond/rupture command");
        ndata = 1;
        iarg += 2;
    } else if (strcmp(style_name,"prob/slip") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_slip = 1; 
        k0 = utils::numeric(FLERR,arg[iarg],false,lmp);
        f0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (k0 < 0.0) error->all(FLERR,"Illegal fix bond/rupture command");
        ndata = 2;
        iarg += 2;
    } else if (strcmp(style_name,"prob/slip/catch") == 0) {
        if (iarg+4 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_slip_catch = 1; 
        ks0 = utils::numeric(FLERR,arg[iarg],false,lmp);
        kc0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        fs0 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        fc0 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        if (ks0 < 0.0 || kc0 < 0.0) error->all(FLERR,"Illegal fix bond/rupture command");
        ndata = 4;
        iarg += 4;
    } else if (strcmp(style_name,"prob/rate") == 0) {
        if (iarg+1 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_rate = 1; 
        k0 = utils::numeric(FLERR,arg[iarg],false,lmp);
        ndata = 1;
        iarg += 1;
    } else error->all(FLERR,"Illegal fix bond/rupture style");

  // Parse remaining keyword arguments
  while (iarg < narg) {
    if (strcmp(arg[iarg],"bond/table") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
      flag_table = 1;
      tables = (Table *) memory->srealloc(tables, (ntables + 1) * sizeof(Table), "bond:tables");
      Table *tb = &tables[ntables];
      null_table(tb);
      if (comm->me == 0) read_table(tb, arg[iarg+1], arg[iarg+2]);
      bcast_table(tb);
      iarg += 3;
    } else if (strcmp(arg[iarg],"bond/distribution") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
      flag_distibution = 1;
      // for now this style does not take any additional args
      iarg += 1;
    } else if (strcmp(arg[iarg],"bond/crit") == 0) {
     if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
      flag_crit = 1;
      double rcrit = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (rcrit < 0.0) error->all(FLERR,"Illegal fix bond/rupture command");
      rcritsq = rcrit*rcrit;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix bond/rupture command");
  }

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + me); 

  if (flag_table && flag_distibution)
    error->all(FLERR,"Cannot use argument bond/table with argument bond/distribution");
 
  // Set forward communication size
  comm_forward = 1+atom->maxspecial;

  // create a unique id for this fix's private bond history instance
  update_flag = 1;
  id_fix_bond_history = utils::strdup(fmt::format("HISTORY_BOND_RUPTURE_{}", instance_total));
  
}

/* ---------------------------------------------------------------------- */

FixBondRupture::~FixBondRupture()
{
  if (fix_bond_history && modify->nfix) modify->delete_fix(id_fix_bond_history);
  delete[] id_fix_bond_history;
  if (setflag) memory->destroy(setflag);
  if (tabindex) memory->destroy(tabindex);

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);
  
}

/* ---------------------------------------------------------------------- */

int FixBondRupture::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::init()
{
  // Allocate and set type flags for FixBondHistory on first call
  if (!setflag) {
    int np1 = atom->nbondtypes + 1;
    memory->create(setflag, np1, "fix_bond_rupture:setflag");
    memory->create(tabindex, np1, "fix_bond_rupture:tabindex");
    
    for (int i = 0; i < np1; i++) setflag[i] = 0;
    if (btype >= 1 && btype <= atom->nbondtypes) setflag[btype] = 1;
  }

  // Create private BOND_HISTORY fix if not yet created
  if (!fix_bond_history) {
    Fix *f = modify->add_fix(fmt::format("{} all BOND_HISTORY {} {}", id_fix_bond_history, update_flag, ndata), 1);
    fix_bond_history = dynamic_cast<FixBondHistory *>(f);
  }

  if (fix_bond_history) fix_bond_history->setflag = setflag;
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double FixBondRupture::store_bond(int n, int i, int j)
{
  int type;
  double delx, dely, delz, r;
  double **x = atom->x;
  double dt = update->dt;
  double **bondstore = fix_bond_history->bondstore;
  tagint *tag = atom->tag;

  int **bond_type = atom->bond_type;

  if (flag_table) error->all(FLERR,"Cannot add bonds with fix bond/rupture using bond/table style");
  
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  r = sqrt(delx * delx + dely * dely + delz * delz);
  
  // Disable store_bond
  //bondstore[n][0] = r; // temporary to test

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) { 
        //fix_bond_history->update_atom_value(i, m, 0, r); // r0
      }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) { 
        //fix_bond_history->update_atom_value(j, m, 0, r); //r0
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

void FixBondRupture::store_data()
{
  int i, j, n, m, type;
  int iatom, jatom, ia, ja;
  double delx, dely, delz, r;
  double **x = atom->x;
  int **bond_type = atom->bond_type;
  double **bondstore = fix_bond_history->bondstore;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;

  // if per/bond data store
  if (flag_table || flag_distibution) {

    // Initialize bondstore by looping over the bondlist
    for (i = 0; i < atom->nlocal; i++) {
        for (m = 0; m < atom->num_bond[i]; m++) {
        type = bond_type[i][m];

        //Skip if bond was turned off
        if (type < 0) continue;

        // map to find index n
        j = atom->map(atom->bond_atom[i][m]);
        if (j == -1) error->one(FLERR, "Atom missing in BPM bond");
        
        // if table style, grab data from table
        // need to lookup the bond in table
        ia = atom->tag[i];
        ja = atom->tag[j];
        
        const Table *tb = &tables[tabindex[type]];

        for (n = 0; n < tb->ninput; n++) {
            iatom = tb->iatomfile[n];
            jatom = tb->jatomfile[n];
            if ((iatom == ia && jatom == ja) || (iatom == ja && jatom == ia)) {
                break; 
            }
        }

        for (int l = 0; l < ndata; l++) {
          double value = tb->datafile[n*ndata + l]; // pull value from table file for this bond and this data index
          fix_bond_history->update_atom_value(j, m, l, value);  
        }

        // if distribution style, grab from distribution
        if (flag_distibution) { 
            // 
        }
        
        //delx = x[i][0] - x[j][0];
        //dely = x[i][1] - x[j][1];
        //delz = x[i][2] - x[j][2];

        // Get closest image in case bonded with ghost
        //domain->minimum_image(FLERR,delx, dely, delz);
        //r = sqrt(delx * delx + dely * dely + delz * delz);
        //r = 100; // temporary to test

        //fix_bond_history->update_atom_value(i, m, 0, r);

        //bondstore[m][0] = r;

        }
    }

    fix_bond_history->post_neighbor();
  }
  
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::post_integrate()
{
  // On first call, initialize stored bond data
  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();
  }

  double **x = atom->x;
  Bond *bond = force->bond;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  double **bondstore = fix_bond_history->bondstore;
  double r0;
  double dt = (update->dt);

  comm->forward_comm();

  int break_count = 0;
  for (int n = 0; n < nbondlist; n++) {
    if (bondlist[n][2] <= 0) continue;

    int i1 = bondlist[n][0];
    int i2 = bondlist[n][1];
    int type = bondlist[n][2];

    if (type != btype) continue;

    if (atom->tag[i2] < atom->tag[i1]) {
      int itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    // If per/bond data, get style modifiers from bondstore
    if (flag_table || flag_distibution) {
        // get style modifiers from bondstore

        if (flag_dist) {
            rcritsq = bondstore[n][0]*bondstore[n][0];
        } else if (flag_fraction){
            p_fraction = bondstore[n][0]; 
        } else if (flag_slip){
            k0 = bondstore[n][0];
            f0 = bondstore[n][1];
        } else if (flag_slip_catch){
            ks0 = bondstore[n][0];
            kc0 = bondstore[n][1];
            fs0 = bondstore[n][2];
            fc0 = bondstore[n][3];
        } else if (flag_rate){
            k0 = bondstore[n][0];
        }     
    }

    // Bond length
    double delx = x[i1][0] - x[i2][0];
    double dely = x[i1][1] - x[i2][1];
    double delz = x[i1][2] - x[i2][2];
    domain->minimum_image(FLERR,delx, dely, delz);
    double rsq = delx*delx + dely*dely + delz*delz;

    // Rupture probability
    double p_rupture = 0.0; // default to no rupture

    // Random number for stochastic rupture
    double probability = random->uniform();
    
    // Modifiers to rupture probability based on style
    if (flag_dist){
        if (rsq >= rcritsq) {
            p_rupture = 1.0;
        }
    }
    if (flag_fraction){
        p_rupture = p_fraction;
    }
    if (flag_slip){
        // Find force in bond
        double fbond; // fbond is returned as f/r
        double engpot = bond->single(btype,rsq,i1,i2,fbond);
        double r = sqrt(rsq);
        double bondforce = fabs(fbond)*r;

        double kr = k0*exp(fabs(bondforce)/f0);
        p_rupture = 1.0 - exp(-kr*dt);
    }
    if (flag_slip_catch){
        // Find force in bond
        double fbond; // fbond is returned as f/r
        double engpot = bond->single(btype,rsq,i1,i2,fbond);
        double r = sqrt(rsq);
        double bondforce = fabs(fbond)*r;

        double kr = ks0*exp(fabs(bondforce)/fs0) + kc0*exp(-fabs(bondforce)/fc0);
        p_rupture = 1.0 - exp(-kr*dt);
    }
    if (flag_rate){
        p_rupture = 1.0 - exp(-k0*dt);
    }

    // if dist overide is on, see if dist crit is met
    if (flag_crit){
        if (rsq >= rcritsq) {
            p_rupture = 1.0;
        }
    }
    
    if (p_rupture < probability) continue; // bond does not rupture

    break_count = 1;
    bondlist[n][2] = 0;
    process_broken(i1, i2);
  }

  int break_all = 0;
  MPI_Allreduce(&break_count, &break_all, 1, MPI_INT, MPI_MAX, world);

  if (break_all == 1) next_reneighbor = update->ntimestep;
  if (break_all == 0) return;

  update_special();
  comm->forward_comm(this);
  update_topology();
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::update_special()
{
  int i, j, m, n1;
  tagint tagi, tagj;
  int nlocal = atom->nlocal;

  tagint *slist;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  for (auto const &it : new_broken_pairs) {
    tagi = it.first;
    tagj = it.second;
    i = atom->map(tagi);
    j = atom->map(tagj);

    if (i < 0 || j < 0) {
      error->one(FLERR,"Fix bond/rupture needs ghost atoms from further away");
    }

    if (i < nlocal) {
      slist = special[i];
      n1 = nspecial[i][0];
      for (m = 0; m < n1; m++)
        if (slist[m] == tagj) break;
      for (; m < n1 - 1; m++) slist[m] = slist[m + 1];
      nspecial[i][0]--;
      nspecial[i][1] = nspecial[i][2] = nspecial[i][0];
    }

    if (j < nlocal) {
      slist = special[j];
      n1 = nspecial[j][0];
      for (m = 0; m < n1; m++)
        if (slist[m] == tagi) break;
      for (; m < n1 - 1; m++) slist[m] = slist[m + 1];
      nspecial[j][0]--;
      nspecial[j][1] = nspecial[j][2] = nspecial[j][0];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::update_topology()
{
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;

  for (int ilist = 0; ilist < neighbor->nlist; ilist++) {
    NeighList *list = neighbor->lists[ilist];
    if (list->copy) continue;

    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    for (auto const &it : new_broken_pairs) {
      tagint tag1 = it.first;
      tagint tag2 = it.second;
      int i1 = atom->map(tag1);
      int i2 = atom->map(tag2);

      if (i1 < 0 || i2 < 0) {
        error->one(FLERR,"Fix bond/rupture needs ghost atoms from further away");
      }

      if (i1 < nlocal) {
        int *jlist = firstneigh[i1];
        int jnum = numneigh[i1];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= SPECIALMASK;
          if (tag[j] == tag2) jlist[jj] = j;
        }
      }

      if (i2 < nlocal) {
        int *jlist = firstneigh[i2];
        int jnum = numneigh[i2];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= SPECIALMASK;
          if (tag[j] == tag1) jlist[jj] = j;
        }
      }
    }
  }

  new_broken_pairs.clear();
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::process_broken(int i, int j)
{
  auto tag_pair = std::make_pair(atom->tag[i], atom->tag[j]);
  new_broken_pairs.push_back(tag_pair);

  int m, n;
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;

  // Manually search and remove from atom arrays
  // When a bond is removed, shift bond history data and delete from the tail

  if (i < nlocal) {
    for (m = 0; m < num_bond[i]; m++) {
      if (bond_atom[i][m] == tag[j] && setflag[bond_type[i][m]]) {
        n = num_bond[i];
        bond_type[i][m] = bond_type[i][n - 1];
        bond_atom[i][m] = bond_atom[i][n - 1];
        // Shift and delete from private bond history
        if (fix_bond_history) {
          fix_bond_history->shift_history(i, m, n - 1);
          fix_bond_history->delete_history(i, n - 1);
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
        // Shift and delete from private bond history
        if (fix_bond_history) {
          fix_bond_history->shift_history(j, m, n - 1);
          fix_bond_history->delete_history(j, n - 1);
        }
        num_bond[j]--;
        break;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixBondRupture::pack_forward_comm(int n, int *list, double *buf,
                                    int /*pbc_flag*/, int * /*pbc*/)
{
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    int ns = nspecial[j][0];
    buf[m++] = ubuf(ns).d;
    for (int k = 0; k < ns; k++) {
        buf[m++] = ubuf(special[j][k]).d;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::unpack_forward_comm(int n, int first, double *buf)
{
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    int ns = (int) ubuf(buf[m++]).i;
    nspecial[i][0] = ns;
    for (int j = 0; j < ns; j++) {
        special[i][j] = (tagint) ubuf(buf[m++]).i;
    }
  }
}

/* ----------------------------------------------------------------------
    read from table file
 ------------------------------------------------------------------------- */

void FixBondRupture::null_table(Table *tb)
{
  tb->datafile = nullptr; tb->iatomfile = nullptr; tb->jatomfile = nullptr;
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::free_table(Table *tb)
{
  memory->destroy(tb->iatomfile);
  memory->destroy(tb->jatomfile);
  memory->destroy(tb->datafile);

  memory->destroy(tb->data);
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void FixBondRupture::read_table(Table *tb, char *file, char *keyword)
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
  memory->create(tb->datafile, ndata * tb->ninput, "bond:datafile");

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
      tb->iatomfile[i] = values.next_int();
      tb->jatomfile[i] = values.next_int();
      for (int j = 0; j < ndata; j++) {
        tb->datafile[i*ndata + j] = values.next_double();
      }
    } catch (TokenizerException &e) {
      error->one(FLERR, "Error parsing bond table '{}' line {} of {}. {}\nLine was: {}", keyword,
                 i + 1, tb->ninput, e.what(), line);
    }
  }

  printf("Read parameters for %i bonds from table\n",tb->ninput);

}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ r0
   N is required, other params are optional
------------------------------------------------------------------------- */

void FixBondRupture::param_extract(Table *tb, char *line)
{
  int nbondlist = neighbor->nbondlist;
  tb->ninput = 0;

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
  if (tb->ninput > nbondlist) error->one(FLERR, "Bond table N cannot exceed number of bonds");
  
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile,ffile,fpflag,fplo,fphi,r0
------------------------------------------------------------------------- */

void FixBondRupture::bcast_table(Table *tb) // *UPDATED
{
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);

  int me;
  MPI_Comm_rank(world, &me);
  if (me > 0) {
    memory->create(tb->iatomfile, tb->ninput, "bond:iatomfile");
    memory->create(tb->jatomfile, tb->ninput, "bond:jatomfile");
    memory->create(tb->datafile, ndata * tb->ninput, "bond:datafile");
  }
  MPI_Bcast(tb->iatomfile, tb->ninput, MPI_INT, 0, world);
  MPI_Bcast(tb->jatomfile, tb->ninput, MPI_INT, 0, world);
  MPI_Bcast(tb->datafile, ndata * tb->ninput, MPI_DOUBLE, 0, world);
}
