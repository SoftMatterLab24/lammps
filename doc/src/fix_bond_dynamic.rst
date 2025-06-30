.. index:: fix bond/dynamic

fix bond/dynamic command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style Nevery itype jtype bondtype ka kd Rcut keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* Nevery = attempt bond attachment and dettachment every this many steps
* itype,jtype = atoms of itype can bond to atoms of jtype (1-Ntypes or type label)
* bondtype = type of bond modified by this fix
* ka = attachment rate
* kd = dettachment rate
* Rcut = two atoms separated by less than Rcut can attach or dettach (distance units)
* zero or more keyword/value pairs may be appended to args
* keyword = *maxbond* or *seed* or *prob* or *mol* or *critical* or *rouse* or *bell* or *catch*

  .. parsed-literal::

       *maxbond* values = maxbond
         maxbond = max # of bonds of bondtype the itype and jtype atoms can have
       *seed* values = seed
         seed = random number seed (positive integer)
       *prob* values = Pattach Pdettach
         Pattach = create a bond with this proabilitiy if otherwise eligible (fraction)
         Pdettach = remove a bond with this proabilitiy if otherwise eligible (fraction)
       *mol* values = 0 or 1 or 2
         0 = any atom can bond if otherwise eligible (default)
         1 = only atoms on different molecules can bond
         2 = only atoms on same molecules can bond
       *critical* values rcrit
         rcrit = length at which bonds permanently break (distance units)
       *rouse* values b0
         b0 = 
       *bell* values f0 kmax
         f0 =
         kmax =
       *catch* values fs0 fc0 kc0
         fs0 =
         fc0 =
         kc0 =

Examples
""""""""

.. code-block:: LAMMPS

   fix 5 all bond/create 10 1 2 0.8 1
   fix 5 all bond/create 1 3 3 0.8 1 prob 0.5 85784 iparam 2 3
   fix 5 all bond/create 1 3 3 0.8 1 prob 0.5 85784 iparam 2 3 atype 1 dtype 2
   fix 5 all bond/create/angle 10 1 2 1.122 1 aconstrain 120 180 prob 1 4928459 iparam 2 1 jparam 2 2

   labelmap atom 1 c1 2 n2
   labelmap bond 1 c1-n2
   fix 5 all bond/create 10 c1 n2 0.8 c1-n2

Description
"""""""""""

Dynamically attach and/or dettach bonds between pairs of atoms as a
simulation runs according to specified criteria. This can be used to
model the cross-linking of polymers, the formation of a percolation
network, or the continual topological reformation of transient networks, etc. 
In this context, a bond means an interaction between a pair of atoms computed
by the :doc:`bond_style <bond_style>` command. This process is different from 
:doc:`fix bond/create <fix_bond_create>` in that bonds are not permanently
created. The fix first establishes whether any bonds should be removed before 
attempting to create new bonds.

For bond removal, a check for possible bond breaking is performed every *Nevery* time steps.
If two atoms :math:`i` and :math:`j` are within a distance *Rcut* of each
other, atom :math:`i` is of type *itype*, atom :math:`j` is of type *jtype*,
and both :math:`i` and :math:`j` are in the specified fix group, then if a bond is of type *btype*
then the bond is labeled as a "possible" bond break. An eligible bond will break with a 
discrete dettachement probability defined as:

.. math::

   \delta P_d = 1 - \exp{\left( -kd \Delta t \right) } 

where :math:`kd` is the dettachment rate, and :math:`\Delta t` is the timestep. For every 
eligible bond if the probability constraint is satisfied then the bond is removed, otherwise it
remains.

For bond creation, a check for possible new bonds is performed every *Nevery* time steps.
If two atoms :math:`i` and :math:`j` are within a distance *Rcut* of each
other, atom :math:`i` is of type *itype*, atom :math:`j` is of type *jtype*,
and both :math:`i` and :math:`j` are in the specified fix group, then if a bond
does not already exist between atoms :math:`i` and :math:`j`, and if both
:math:`i` and :math:`j` meet their respective *maxbond* requirements (explained
below), then :math:`i` and :math:`j` are labeled as a "possible" bond pair. 
If several atoms are close to an atom, it may have multiple possible bond partners.
An eligible bond will form with a discrete attachement probability defined as:

.. math::

   \delta P_a = 1 - \exp{\left( -ka \Delta t \right) } 

where :math:`kd` is the dettachment rate, and :math:`\Delta t` is the timestep. If the
probability constraint is satisfied, then the bond will be formed. Note that with this
method, each atom may be part of multiple created bonds on a given time step.

It is permissible to have *itype* = *jtype*\ .  *Rcut* must be :math:`\leq` the
pair-wise cutoff distance between *itype* and *jtype* atoms, as defined
by the :doc:`pair_style <pair_style>` command.


#### KEYWORDS

The *maxbond* keyword can be used to limit the number of bonds allowed. 
If either atom :math:`i` of type *itype* or atom :math:`j` of type *jtype* 
has *maxbond* bonds (set by value), then a new bond will not be formed.

.. note::

    To create a new bond, the internal LAMMPS data structures that
    store this information must have space for it.  When LAMMPS is
    initialized from a data file, the list of bonds is scanned and the
    maximum number of bonds per atom is tallied.  If some atom will
    acquire more bonds than this limit as this fix operates, then the
    "extra bond per atom" parameter must be set to allow for it. Therefore, 
    the value of *maxbond* should be less than or equal to what the "extra 
    bond per atom" parameter has been set to. See the :doc:`read_data <read_data>` 
    or :doc:`create_box <create_box>` command for more details. 

The *seed* keyword can be used to specifiy the processor-unique seed 
used to initialized the Marsaglia random number generator. By default the
seed is 12345. The *value* setting must be a positive integer. 

The *prob* keyword can be used to directly set the attachment and
dettachement probabilities. Both the *Pattach* and *Pdettach* settings
must be values between 0.0 and 1.0. For bond creation a uniform random 
number between 0.0 and 1.0 is generated and the eligible bond 
is only created if the random number is less than *Pattach*. Likewise,
for bond deletion a uniform random  number between 0.0 and 1.0 is generated 
and the eligible bond  is only removed if the random number is less than *Pdettach*.
The *prob* keyword cannot be used with keyword *rouse* or *bell* or *catch*.

The *mol* keyword can be used to limit the bonding functionality of
the participating atoms. If the value of *mol* is 1, then in addition 
to the previous constraints, atom :math:`i` will also check to see if 
atom :math:`j` belongs to a different molecule. If this is true, and the
other conditions are met, then :math:`i` and :math:`j` are labeled as a 
"possible" bond pair. If the value of *mol* is 2 then atom :math:`i` will only
label atom :math:`j` as a possible pair if they both belong to the same
molecule. By default *mol* value is 0, in which case atoms do not
check what molecule they belong to.

The *critical* keyword specifies whether bonds permanently rupture after
exceeding a critial legnth set by value *rcrit*. This process is performed
before dettachement and attachment, such that when a bond exceeds its 
critical length, it is immediately labeled for removal. 

The *iparam* and *jparam* keywords can be used to limit the bonding
functionality of the participating atoms.  Each atom keeps track of
how many bonds of *bondtype* it already has.  If atom :math:`i` of type
*itype* already has *maxbond* bonds (as set by the *iparam*
keyword), then it will not form any more, and likewise for atom :math:`j`.
If *maxbond* is set to 0, then there is no limit on the number of bonds
that can be formed with that atom.

The *newtype* value for *iparam* and *jparam* can be used to change
the atom type of atom :math:`i` or :math:`j` when it reaches *maxbond* number
of bonds of type *bondtype*\ .  This means it can now interact in a pair-wise
fashion with other atoms in a different way by specifying different
:doc:`pair_coeff <pair_coeff>` coefficients.  If you do not wish the
atom type to change, simply specify *newtype* as *itype* or *jtype*\ .

The *prob* keyword can also affect whether an eligible bond is
actually created.  The *fraction* setting must be a value between 0.0
and 1.0.  A uniform random number between 0.0 and 1.0 is generated and
the eligible bond is only created if the random number is less than *fraction*.

The *aconstrain* keyword is only available with the fix
bond/create/angle command.  It allows one to specify minimum and maximum
angles *amin* and *amax*, respectively, between the two prospective bonding
partners and a third particle that is already bonded to one of the two
partners. Such a criterion can be important when new angles are defined
together with the formation of a new bond.  Without a restriction on the
permissible angle, and for stiffer angle potentials, very large energies
can arise and lead to unphysical behavior.

Any bond that is created is assigned a bond type of *bondtype*.

When a bond is created, data structures within LAMMPS that store bond
topologies are updated to reflect the creation.  If the bond is part of
new 3-body (angle) or 4-body (dihedral, improper) interactions, you
can choose to create new angles, dihedrals, and impropers as well using
the *atype*, *dtype*, and *itype* keywords.  All of these changes
typically affect pair-wise interactions between atoms that are now part
of new bonds, angles, etc.

.. note::

   One data structure that is not updated when a bond breaks are
   the molecule IDs stored by each atom.  Even though two molecules
   become one molecule due to the created bond, all atoms in the new
   molecule retain their original molecule IDs.

If the *atype* keyword is used and if an angle potential is defined
via the :doc:`angle_style <angle_style>` command, then any new 3-body
interactions inferred by the creation of a bond will create new angles
of type *angletype*, with parameters assigned by the corresponding
:doc:`angle_coeff <angle_coeff>` command.  Likewise, the *dtype* and
*itype* keywords will create new dihedrals and impropers of type
*dihedraltype* and *impropertype*\ .

.. note::

   To create a new bond, the internal LAMMPS data structures that
   store this information must have space for it.  When LAMMPS is
   initialized from a data file, the list of bonds is scanned and the
   maximum number of bonds per atom is tallied.  If some atom will
   acquire more bonds than this limit as this fix operates, then the
   "extra bond per atom" parameter must be set to allow for it.  Ditto
   for "extra angle per atom", "extra dihedral per atom", and "extra
   improper per atom" if angles, dihedrals, or impropers are being added
   when bonds are created.  See the :doc:`read_data <read_data>` or
   :doc:`create_box <create_box>` command for more details.  Note that a
   data file with no atoms can be used if you wish to add non-bonded
   atoms via the :doc:`create atoms <create_atoms>` command (e.g., for a
   percolation simulation).

.. note::

   LAMMPS stores and maintains a data structure with a list of the
   first, second, and third neighbors of each atom (within the bond topology of
   the system) for use in weighting pair-wise interactions for bonded
   atoms.  Note that adding a single bond always adds a new first neighbor
   but may also induce **many** new second and third neighbors, depending on the
   molecular topology of your system.  The "extra special per atom"
   parameter must typically be set to allow for the new maximum total
   size (first + second + third neighbors) of this per-atom list.  There are two
   ways to do this.  See the :doc:`read_data <read_data>` or
   :doc:`create_box <create_box>` commands for details.

.. note::

   Even if you do not use the *atype*, *dtype*, or *itype*
   keywords, the list of topological neighbors is updated for atoms
   affected by the new bond.  This in turn affects which neighbors are
   considered for pair-wise interactions, using the weighting rules set by
   the :doc:`special_bonds <special_bonds>` command.  Consider a new bond
   created between atoms :math:`i` and :math:`j`.  If :math:`j` has a bonded
   neighbor :math:`k`, then :math:`k` becomes a second neighbor of :math:`i`.
   Even if the *atype* keyword is not used to create angle :math:`\angle ijk`,
   the pair-wise interaction between :math:`i` and :math:`k` could potentially
   be turned off or weighted by the 1--3 weighting specified
   by the :doc:`special_bonds <special_bonds>` command.  This is the case
   even if the "angle yes" option was used with that command.  The same
   is true for third neighbors (1--4 interactions), the *dtype* keyword, and
   the "dihedral yes" option used with the
   :doc:`special_bonds <special_bonds>` command.

Note that even if your simulation starts with no bonds, you must
define a :doc:`bond_style <bond_style>` and use the
:doc:`bond_coeff <bond_coeff>` command to specify coefficients for the
*bondtype*\ .  Similarly, if new atom types are specified by the
*iparam* or *jparam* keywords, they must be within the range of atom
types allowed by the simulation and pair-wise coefficients must be
specified for the new types.

Computationally, each time step this fix is invoked, it loops over
neighbor lists and computes distances between pairs of atoms in the
list.  It also communicates between neighboring processors to
coordinate which bonds are created.  Moreover, if any bonds are
created, neighbor lists must be immediately updated on the same
time step.  This is to ensure that any pair-wise interactions that
should be turned "off" due to a bond creation, because they are now
excluded by the presence of the bond and the settings of the
:doc:`special_bonds <special_bonds>` command, will be immediately
recognized.  All of these operations increase the cost of a time step.
Thus, you should be cautious about invoking this fix too frequently.

You can dump out snapshots of the current bond topology via the :doc:`dump local <dump>` command.

.. note::

   Creating a bond typically alters the energy of a system.  You
   should be careful not to choose bond creation criteria that induce a
   dramatic change in energy.  For example, if you define a very stiff
   harmonic bond and create it when two atoms are separated by a distance
   far from the equilibrium bond length, then the two atoms will oscillate
   dramatically when the bond is formed.  More generally, you may need to
   thermostat your system to compensate for energy changes resulting from
   created bonds (and angles, dihedrals, impropers).

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.

This fix computes two statistics which it stores in a global vector of
length 2, which can be accessed by various :doc:`output commands
<Howto_output>`.  The vector values calculated by this fix are
"intensive".

The two quantities in the global vector are

  (1) number of bonds created on the most recent creation time step
  (2) cumulative number of bonds created

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`fix bond/break <fix_bond_break>`, :doc:`fix bond/react <fix_bond_react>`, :doc:`fix bond/swap <fix_bond_swap>`,
:doc:`dump local <dump>`, :doc:`special_bonds <special_bonds>`

Default
"""""""

The option defaults are iparam = (0,itype), jparam = (0,jtype), and
prob = 1.0.
