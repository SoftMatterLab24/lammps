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
#include <unordered_map>
#include "math_const.h"
#include "random_mars.h"

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBondRupture::FixBondRupture(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  bprob(nullptr), random(nullptr), 
  id_fix_bond_history_rupture(nullptr), fix_bond_history(nullptr)
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
  flag_table = 0; flag_distribution = 0; flag_crit = 0, flag_prob = 0;

  btype = utils::inumeric(FLERR,arg[3],false,lmp);
  // nevery -> future version
  const char *style_name = arg[4];
  std::vector<std::string> param_names;

  if (btype < 1 || btype > atom->nbondtypes)
  error->all(FLERR,"Invalid bond type in fix bond/rupture command");

    int iarg = 5;
    // Parse style-specific values
    if (strcmp(style_name,"dist") == 0) {
        if (iarg+1 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_dist = 1; param_names = {"rcrit"};
        rcrit = utils::numeric(FLERR,arg[iarg],false,lmp);
        if (rcrit < 0.0) error->all(FLERR,"Illegal fix bond/rupture command");
        rcritsq = rcrit*rcrit;
        ndata = 1;
        iarg += 1;
    } else if (strcmp(style_name,"prob/fraction") == 0) {
        if (iarg+1 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_fraction = 1; flag_prob = 1; param_names = {"fraction"};
        p_fraction = utils::numeric(FLERR,arg[iarg],false,lmp);
        if (p_fraction < 0.0 || p_fraction > 1.0) error->all(FLERR,"Illegal fix bond/rupture command");
        ndata = 1;
        iarg += 1;
    } else if (strcmp(style_name,"prob/slip") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_slip = 1; flag_prob = 1; param_names = {"ks0","f0"};
        k0 = utils::numeric(FLERR,arg[iarg],false,lmp);
        f0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (k0 < 0.0) error->all(FLERR,"Illegal fix bond/rupture command");
        ndata = 2;
        iarg += 2;
    } else if (strcmp(style_name,"prob/slip/catch") == 0) {
        if (iarg+4 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_slip_catch = 1; flag_prob = 1; param_names = {"ks0","kc0","fs0","fc0"};
        ks0 = utils::numeric(FLERR,arg[iarg],false,lmp);
        kc0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        fs0 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        fc0 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        if (ks0 < 0.0 || kc0 < 0.0) error->all(FLERR,"Illegal fix bond/rupture command");
        ndata = 4;
        iarg += 4;
    } else if (strcmp(style_name,"prob/rate") == 0) {
        if (iarg+1 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
        flag_rate = 1; flag_prob = 1; param_names = {"kr"};
        k0 = utils::numeric(FLERR,arg[iarg],false,lmp);
        ndata = 1;
        iarg += 1;
    } else error->all(FLERR,"Illegal fix bond/rupture style");

    // initialze size of pdf parameter arrays based on style
    lo = new double[ndata];
    hi = new double[ndata];
    mu = new double[ndata];
    sigma = new double[ndata];
    lambda = new double[ndata];
    alpha = new double[ndata];
    beta = new double[ndata];

    use_dist = new int[ndata];
    for (int i = 0; i < ndata; ++i) use_dist[i] = 0;   // 0 = global deterministic, 1 = draw from distribution

    // Create a mapping of parameter name to index for lookup when parsing keyword args
    std::unordered_map<std::string,int> pindex;
    for (int i = 0; i < (int)param_names.size(); ++i) pindex[param_names[i]] = i;

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
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
      flag_distribution = 1;
      
      dist_type = utils::strdup(arg[iarg+1]);
      
      int need = 0; // the number of distribution parameters needed for this distribution style
      if (strcmp(dist_type,"uniform") == 0 ) {
        need = 2;
      } else if (strcmp(dist_type,"gauss") == 0) {
        need = 2;
      } else if (strcmp(dist_type,"exponential") == 0) {
        need = 1;
      } else if (strcmp(dist_type,"weibull") == 0) {
        need = 2;
      } else error->all(FLERR,"Illegal fix bond/rupture command");

      iarg += 2;

      // parse distribution parameters for each style parameter index until run out of args or hit another keyword
      while (iarg < narg && !(strcmp(arg[iarg],"bond/table") == 0 || strcmp(arg[iarg],"critical") == 0 || strcmp(arg[iarg],"seed") == 0)) {
    
        // must be a known parameter name for this rupture style
        auto it = pindex.find(arg[iarg]);
        if (it == pindex.end())
            error->all(FLERR,"Illegal fix bond/rupture bond/distribution parameter name");

        int idx = it->second;
        if (use_dist[idx])
            error->all(FLERR,"Duplicate bond/distribution entry for same parameter");

        if (iarg + 1 + need > narg)
        error->all(FLERR,"Not enough arguments for bond/distribution parameter block");

        use_dist[idx] = 1;
        // parse the distribution parameters for this parameter index
        if (strcmp(dist_type,"uniform") == 0 ) {
            lo[idx] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
            hi[idx] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
            if (hi[idx] <= lo[idx]) error->all(FLERR,"Illegal uniform distribution bounds");
            iarg += 3;
        } else if (strcmp(dist_type,"gauss") == 0) {
            mu[idx]    = utils::numeric(FLERR,arg[iarg+1],false,lmp);
            sigma[idx] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
            if (sigma[idx] <= 0.0) error->all(FLERR,"Illegal gauss sigma");
            iarg += 3;
        } else if (strcmp(dist_type,"exponential") == 0) {
            lambda[idx] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
            if (lambda[idx] <= 0.0) error->all(FLERR,"Illegal exponential lambda");
            iarg += 2;
        } else if (strcmp(dist_type,"weibull") == 0) {
            alpha[idx] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // shape or scale—your convention
            beta[idx]  = utils::numeric(FLERR,arg[iarg+2],false,lmp);
            if (alpha[idx] <= 0.0 || beta[idx] <= 0.0) error->all(FLERR,"Illegal weibull params");
            iarg += 3;
        } else error->all(FLERR,"Illegal fix bond/rupture command");
      }
    } else if (strcmp(arg[iarg],"critical") == 0) {
     if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/rupture command");
      flag_crit = 1;
      double rcrit_g = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (rcrit < 0.0) error->all(FLERR,"Illegal fix bond/rupture command");
      rcritsq_g = rcrit_g*rcrit_g;
      iarg += 2;
    } else if (strcmp(arg[iarg],"seed") == 0) {
      seed = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix bond/rupture command"); 
  }

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + me); 

  if (flag_table && flag_distribution)
    error->all(FLERR,"Cannot use argument bond/table with argument bond/distribution");
 
  // Set forward communication size (defaults)
  comm_forward = 1;
  comm_reverse = 1;

  // create a unique id for this fix's private bond history instance
  update_flag = 1;
  id_fix_bond_history_rupture = utils::strdup(fmt::format("HISTORY_BOND_RUPTURE_{}", instance_total));

  // allocate
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

FixBondRupture::~FixBondRupture()
{
  delete random;

  memory->destroy(bprob);
  memory->destroy(break_partner);
  
  if (fix_bond_history) modify->delete_fix(id_fix_bond_history_rupture);
  delete[] id_fix_bond_history_rupture;

  if (dist_type) delete[] dist_type;
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

  // Create fix if not yet created
  if (!fix_bond_history) {
    Fix *f1 = modify->add_fix(fmt::format("{} all BOND_HISTORY {} {}", id_fix_bond_history_rupture, update_flag, ndata), 1);
    fix_bond_history = dynamic_cast<FixBondHistory *>(f1);
  }

  // Get all histories
  histories = modify->get_fix_by_style("BOND_HISTORY");

  if (fix_bond_history) fix_bond_history->setflag = setflag;
}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour) (DEPRECIATED)
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
  int i1, i2;
  int iatom, jatom, ia, ja;
  double delx, dely, delz, r;
  double **x = atom->x;
  int **bond_type = atom->bond_type;
  double **bondstore = fix_bond_history->bondstore;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  tagint itag, jtag;
  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int nbonds = atom->nbonds;

  long int natoms = atom->natoms;
  long int key;

  // Now if per/bond, store data
  if (flag_table || flag_distribution) {

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
        if (flag_table) {
          ia = atom->tag[i];
          ja = atom->tag[j];
        
          const Table *tb = &tables[tabindex[type]];

          int table_idx = -1;
          for (n = 0; n < tb->ninput; n++) {
            iatom = tb->iatomfile[n];
            jatom = tb->jatomfile[n];
            if ((iatom == ia && jatom == ja) || (iatom == ja && jatom == ia)) {
              table_idx = n;
              break; 
            }
          }
          
          if (table_idx >= 0) {
            // Find the bond in the global bondlist to get correct index
            int bn = -1;
            for (n = 0; n < nbondlist; n++) {
              i1 = bondlist[n][0];
              i2 = bondlist[n][1];
              if ((tag[i1] == ia && tag[i2] == ja) || (tag[i1] == ja && tag[i2] == ia)) {
                bn = n;
                break;
              }
            }
            
            if (bn >= 0) {
              for (int l = 0; l < ndata; l++) {
                double value = tb->datafile[table_idx*ndata + l];
                fix_bond_history->update_atom_value(i, m, l, value);
                bondstore[bn][l] = value;  
              }
            }
          }
        }
        
        // if distribution style, grab from distribution
        if (flag_distribution) { 
          // tags of atoms 
          ia = atom->tag[i];
          ja = atom->tag[j];

          for (int l = 0; l < ndata; l++) {
            double value = 0.0;
            if (use_dist[l]) {
              key = std::min(ia,ja)*natoms + std::max(ia,ja);
              value = sample_cdf(key, l);
            } else {
              // if not using distribution for this parameter, use the default value specified in args
              if (flag_dist){
                value = rcrit;
              } else if (flag_fraction){
                value = p_fraction;
              } else if (flag_slip){
                value = k0;
              } else if (flag_slip_catch){
                value = ks0;
              } else if (flag_rate){
                value = k0;
              }
            }
            fix_bond_history->update_atom_value(i, m, l, value);
          }
        }
      }
    }
    fix_bond_history->post_neighbor();
  }
  if (me == 0) printf("Finished storing bond data in fix_bond_rupture\n");
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::post_integrate()
{
  // On first call, initialize stored per/bond data
  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();
  }

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  
  tagint **bond_atom = atom->bond_atom;
  tagint *tag = atom->tag;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  Bond *bond = force->bond;
  
  int i, j, m, n;
  int i1, i2, itmp, type;
  double **x = atom->x;
  double **bondstore = fix_bond_history->bondstore;

  double dt = (update->dt);

  // acquire updated ghost atom positions
  // necessary b/c are calling this after integrate, but before Verlet comm

  comm->forward_comm();

  // Resize break partner list and initialize it
  // needs to be atom->nmax in length
  if (atom->nmax > nmax) {
    memory->destroy(break_partner);
    nmax = atom->nmax;
    memory->create(break_partner, nmax, "fix_bond_rupture:break_partner");
  }

  for (i = 0; i < nall; i++) {
    break_partner[i] = 0;
  }

  // loop over bondlist to  set probability of rupture for each bond if using probabilistic rupture
  if (flag_prob) {
    maxbond = atom->bond_per_atom;

    // Resize bond probability list and initialize it
    // Must be size [nmax][maxbond]
    if (atom->nmax > nmax) {
      memory->destroy(bprob);
      nmax = atom->nmax;
      memory->create(bprob, nmax, maxbond, "fix_bond_rupture:bprob");
    }
  
    for (i = 0; i < nall; i++) {
      for (int m = 0; m < maxbond; m++)
        bprob[i][m] = 0.0;
    }
  
    comm_forward = 1 + 2*maxbond;
  
    // For each locally owned atom, 
    // generate random value for each atom with a bond
    // forward comm of partner and random value, so ghosts have it
    for (i = 0; i < nlocal; i++){
      for (int m = 0; m < num_bond[i]; m++) {
        int type = bond_type[i][m];
        if (type < 0) continue;          // turned off
        if (type != btype) continue;
  
        i1 = tag[i];
        i2 = bond_atom[i][m];
  
        if (i1 < i2) {
          bprob[i][m] = random->uniform();
        } else {
        }
      }
    }
  
    commflag = 1;
    comm->forward_comm(this);
    commflag = 0;
  }
  
  for (n = 0; n < nbondlist; n++) {
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    if (type != btype) continue;

    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    if (i1 >= atom->nlocal) continue;  // Only process if i1 is local

    // Get bond probability if using
    probability = 0.0;
    if (flag_prob) {
      
      int b = -1;
      int nb = atom->num_bond[i1];
      if (nb > maxbond) nb = maxbond;
  
      for (m = 0; m < nb; m++) {
        if (atom->bond_atom[i1][m] == tag[i2]) {
          b = m; break;
        }
      }
  
      if (b < 0) continue;
      probability = bprob[i1][b];
    }

    // If per/bond data, get style modifiers from bondstore
    if (flag_table || flag_distribution) {
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
    domain->minimum_image(FLERR, delx, dely, delz);
    double rsq = delx*delx + dely*dely + delz*delz;

    // Rupture probability
    double p_rupture = 0.0; // default to no rupture

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
      if (rsq >= rcritsq_g) {
        p_rupture = 1.0;
      }
    }
    
    // Check the probability constraint
    if (p_rupture <= probability) continue; // bond does not rupture
    
    // Mark bond for potential breaking - mark BOTH atoms symmetrically
    bondlist[n][2] = 0;  // Mark bond as dead in neighbor list
    break_partner[i1] = tag[i2];
    break_partner[i2] = tag[i1];

  }
  
  // Synchronize break flags across processors via reverse communication
  // This ensures both atoms agree on breaking when one is ghost on the other processor
  // Pack 1 values per atom (break_partner)
  if (force->newton_bond) comm->reverse_comm(this);

  // foward comm of break_partner so ghosts have it for next step
  commflag = 2;
  comm->forward_comm(this);
  commflag = 0;
  
  // Process bond breaking
  for (i = 0; i < nlocal; i++) {
    if (break_partner[i] == 0) continue;  // No partner, so skip
    j = atom->map(break_partner[i]);
    if (break_partner[j] != tag[i]) continue;  // Both atoms agree to break, so delete the bond
    
    process_broken(i, j);
  }

  update_topology();
  
  // Force reneighboring if any bonds broke
  next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::process_broken(int i, int j) 
{
  
  // Then add the pair to new_broken_pairs
  auto tag_pair = std::make_pair(atom->tag[i], atom->tag[j]);
  new_broken_pairs.push_back(tag_pair);

  int m, n, l, nmax;

  // Manually search and remove from atom arrays
  // need to remove in case special bonds arrays rebuilt
  int nlocal = atom->nlocal;

  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  
  if (i < nlocal) {
    
    for (m = (num_bond[i] - 1); m >= 0; m--) {
      if (bond_atom[i][m] == tag[j]) {

        nmax = num_bond[i] - 1;
        if (m == nmax) {
          if (n_histories > 0)
            for (auto &ihistory : histories) {
              auto fix_bond_history = dynamic_cast<FixBondHistory *>(ihistory);
              fix_bond_history->delete_history(i, m);
            }
        } else {
          bond_type[i][m] = bond_type[i][nmax];
          bond_atom[i][m] = bond_atom[i][nmax];
          if (n_histories > 0) {
            for (auto &ihistory : histories) {
              auto fix_bond_history = dynamic_cast<FixBondHistory *>(ihistory);
              fix_bond_history->shift_history(i, m, nmax);
              fix_bond_history->delete_history(i, nmax);
            }
          }
        }
        bond_type[i][nmax] = 0;
        num_bond[i]--;
        break;
      }
    }
  }

  if (j < nlocal) {
    
    for (n = (num_bond[j] - 1); n >= 0; n--) {
      if (bond_atom[j][n] == tag[i]) {

        nmax = num_bond[j] - 1;
        if (n == nmax) {
          if (n_histories > 0)
            for (auto &ihistory : histories) {
              auto fix_bond_history = dynamic_cast<FixBondHistory *>(ihistory);
              fix_bond_history->delete_history(j, n);
            }
        } else {
          bond_type[j][n] = bond_type[j][nmax];
          bond_atom[j][n] = bond_atom[j][nmax];
          if (n_histories > 0) {
            for (auto &ihistory : histories) {
              auto fix_bond_history = dynamic_cast<FixBondHistory *>(ihistory);
              fix_bond_history->shift_history(j, n, nmax);
              fix_bond_history->delete_history(j, nmax);
            }
          }
        }
        bond_type[j][nmax] = 0;
        num_bond[j]--;
        break;
      }
    }
  }

  // Update special neighbor list
  tagint *slist;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  // remove i from special bond list for atom j and vice versa
  // ignore n2, n3 since 1-3, 1-4 special factors required to be 1.0
  if (i < nlocal) {
    slist = special[i];
    int n1 = nspecial[i][0];
    int m;
    for (m = 0; m < n1; m++)
      if (slist[m] == tag[j]) break;
    for (; m < n1 - 1; m++) slist[m] = slist[m + 1];
    nspecial[i][0]--;
    nspecial[i][1] = nspecial[i][2] = nspecial[i][0];
  }

  if (j < nlocal) {
    slist = special[j];
    int n1 = nspecial[j][0];
    int m;
    for (int m = 0; m < n1; m++)
      if (slist[m] == tag[i]) break;
    for (; m < n1 - 1; m++) slist[m] = slist[m + 1];
    nspecial[j][0]--;
    nspecial[j][1] = nspecial[j][2] = nspecial[j][0];
  }
}

/* ----------------------------------------------------------------------
  Update special lists for recently broken/created bonds
  Assumes appropriate atom/bond arrays were updated, e.g. had called
      neighbor->add_temporary_bond(i1, i2, btype);
------------------------------------------------------------------------- */

void FixBondRupture::update_topology()
{

  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;

  // In theory could communicate a list of broken bonds to neighboring processors here
  // to remove restriction that users use Newton bond off

  for (int ilist = 0; ilist < neighbor->nlist; ilist++) {
    NeighList *list = neighbor->lists[ilist];

    // Skip copied lists, will update original
    if (list->copy) continue;

    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    for (auto const &it : new_broken_pairs) {
      tagint tag1 = it.first;
      tagint tag2 = it.second;
      int i1 = atom->map(tag1);
      int i2 = atom->map(tag2);

      if (i1 < 0 || i2 < 0) {
        error->one(FLERR,"Fix bond/dynamic needs ghost atoms "
                    "from further away 4");
      }

      // Loop through atoms of owned atoms i j
      if (i1 < nlocal) {
        int *jlist = firstneigh[i1];
        int jnum = numneigh[i1];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= SPECIALMASK;    // Clear special bond bits
          if (tag[j] == tag2) jlist[jj] = j;
        }
      }

      if (i2 < nlocal) {
        int *jlist = firstneigh[i2];
        int jnum = numneigh[i2];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= SPECIALMASK;    // Clear special bond bits
          if (tag[j] == tag1) jlist[jj] = j;
        }
      }
    }
  }

  for (int ilist = 0; ilist < neighbor->nlist; ilist++) {
    NeighList *list = neighbor->lists[ilist];

    // Skip copied lists, will update original
    if (list->copy) continue;

    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    for (auto const &it : new_created_pairs) {
      tagint tag1 = it.first;
      tagint tag2 = it.second;
      int i1 = atom->map(tag1);
      int i2 = atom->map(tag2);

      if (i1 < 0 || i2 < 0) {
        error->one(FLERR,"Fix bond/dynamic needs ghost atoms "
                    "from further away 5");
      }

      // Loop through atoms of owned atoms i j
      if (i1 < nlocal) {
        int *jlist = firstneigh[i1];
        int jnum = numneigh[i1];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          if (((j >> SBBITS) & 3) != 0) continue;               // Skip bonded pairs
          if (tag[j] == tag2) jlist[jj] = j ^ (1 << SBBITS);    // Add 1-2 special bond bits
        }
      }

      if (i2 < nlocal) {
        int *jlist = firstneigh[i2];
        int jnum = numneigh[i2];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          if (((j >> SBBITS) & 3) != 0) continue;               // Skip bonded pairs
          if (tag[j] == tag1) jlist[jj] = j ^ (1 << SBBITS);    // Add 1-2 special bond bits
        }
      }
    }
  }

  new_broken_pairs.clear();
  new_created_pairs.clear();

}

/* ---------------------------------------------------------------------- */

int FixBondRupture::pack_forward_comm(int n, int *list, double *buf,
                                      int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;

  if (commflag == 1) {
    int *num_bond = atom->num_bond;
    tagint **bond_atom = atom->bond_atom;

    for (int ii = 0; ii < n; ii++) {
      int i = list[ii];
      int nb = num_bond[i];
      if (nb > maxbond) nb = maxbond;

      buf[m++] = ubuf(nb).d;
      for (int k = 0; k < nb; k++) {
        buf[m++] = ubuf(bond_atom[i][k]).d; // partner tag
        buf[m++] = bprob[i][k];             // probability draw (or 0 if not owner)
      }
    }
  }

  if (commflag == 2) {
    // forward comm of break_partner for ghosts
    for (int ii = 0; ii < n; ii++) {
      int i = list[ii];
      buf[m++] = ubuf(break_partner[i]).d;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int FixBondRupture::pack_reverse_comm_size(int n, int nswap)
{
  // Each atom sends 1 value: break_partner
  return n;
}

/* ---------------------------------------------------------------------- */

int FixBondRupture::pack_reverse_comm(int n, int ndata, double *buf)
{
  // Pack break_partner for reverse communication
  // This synchronizes which bonds should be broken across processors
  int m = 0;
  int last = ndata + n;
  for (int i = ndata; i < last; i++) {
    buf[m++] = ubuf(break_partner[i]).d;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::unpack_reverse_comm(int n, int *list, double *buf)
{
  // Unpack break_partner from reverse communication
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    break_partner[j] = (tagint) ubuf(buf[m++]).i;
  }
}

/* ---------------------------------------------------------------------- */

void FixBondRupture::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;

  if (commflag == 1) {
    int *num_bond = atom->num_bond;
    tagint **bond_atom = atom->bond_atom;

    int last = first + n;
    for (int i = first; i < last; i++) {
      int nb = (int) ubuf(buf[m++]).i;
      if (nb > maxbond) nb = maxbond;

      // reset
      for (int k = 0; k < maxbond; k++) bprob[i][k] = 0.0;

      // for each received (partner,prob), match to our local bond list
      for (int kk = 0; kk < nb; kk++) {
        tagint ptag = (tagint) ubuf(buf[m++]).i;
        double p = buf[m++];

        // find matching bond slot
        int nb_here = num_bond[i];
        if (nb_here > maxbond) nb_here = maxbond;
        for (int k = 0; k < nb_here; k++) {
          if (bond_atom[i][k] == ptag) {
            bprob[i][k] = p;
            break;
          }
        }
      }
    }
  }

  if (commflag == 2) {
    // forward comm of break_partner for ghosts
    int last = first + n;
    for (int i = first; i < last; i++) {
      break_partner[i] = (tagint) ubuf(buf[m++]).i;
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
------------------------------------------------------------------------- */

void FixBondRupture::bcast_table(Table *tb)
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

/* ----------------------------------------------------------------------
   return a number sampled from the cdf given a unique (i,j) pair key
------------------------------------------------------------------------- */

double FixBondRupture::sample_cdf(long int key, int l) 
{   

  // Draw two random uniform numbers
  double u1 = bond_uniform(key,l);
  double u2 = bond_uniform(key,l + ndata);

  // attempts to draw valid random number
  for (int attempt = 0; attempt < 1000; ++attempt) {
    double x = 0.0;

    if (strcmp(dist_type,"uniform") == 0 ) {
      x = lo[l] + u1 * (hi[l] - lo[l]);
    } else if (strcmp(dist_type,"gauss") == 0) {

      if (u1 < 1.0e-12) u1 = 1.0e-12;

      double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
      x  = mu[l] + sigma[l] * z0;
  
    } else if (strcmp(dist_type,"exponential") == 0) {
      x = -log(1-u1)/lambda[l];
    } else if (strcmp(dist_type,"weibull") == 0) {
        x = alpha[l] * pow(-log(1-u1), 1.0/beta[l]);
    } else {
        error->all(FLERR,"Illegal distribution type in fix bond/rupture");
    }

    // Accept only finite, non-negative values
    if (std::isfinite(x) && x > 0.0) return x;
  }

  error->all(FLERR,"Failed to sample valid bond parameter");
  return 0.0; 
}

/* ----------------------------------------------------------------------
    Psuedo-random uniform number generator given pair key, index, and global seed
  ------------------------------------------------------------------------- */

double FixBondRupture::bond_uniform(long int key, int l)
{
  // 32-bit seed from key and l
  uint32_t x = static_cast<uint32_t>(key);
  uint32_t y = static_cast<uint32_t>(l);

  x ^= static_cast<uint32_t>(seed) * 0x9E3779B9u;

  // Mix key and l (avalanche hash)
  x = x ^ (y * 0x9E3779B9u);
  x = x ^ (x >> 16);
  x = x * 0x7FEB352Du;
  x = x ^ (x >> 15);
  x = x * 0x846CA68Bu;
  x = x ^ (x >> 16);

  // x is now a pseudo-random 32-bit integer.
  // Map to (0,1) excluding endpoints.
  // (x + 1) / (2^32 + 1) is strictly between 0 and 1.
  const double denom = 4294967297.0; // 2^32 + 1
  double u = (static_cast<double>(x) + 1.0) / denom;

  return u;
}