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
* either ka, kd can be a variable (see below)
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
         b0 = diffusion lengthscale (distance units)
       *bell* values f0
         f0 = force-sensitivity of bond (force units)
       *catch* values fs0 fc0 kc0
         fs0 = force-sensitivity of the bonds slip barrier (force units)
         fc0 = force-sensitivity of the bonds catch barrier (force units)
         kc0 = scale factor for bonds catch barrier (unitless)

Examples
""""""""

.. code-block:: LAMMPS

   fix 5 all bond/dynamic 1 2 2 1 10 0.1 1.7
   fix 5 all bond/dynamic 1 2 2 1 10 0.1 1.7 maxbond 8 prob 0.5 0.7 
   fix 5 all bond/dynamic 1 2 2 1 10 0.1 1.7 maxbond 8 critical 1.5
   fix 5 all bond/dynamic 1 2 2 1 v_ka v_kd 1.7 maxbond 8 critical 1.5 rouse 0.5 catch 2 2 1


Description
"""""""""""

Dynamically attach and/or dettach bonds between pairs of atoms as a
simulation runs according to specified criteria. Bond kinetics
(attachment and dettachment) is treated as a stochastic process, where
each event is independent, thereby treated as a Poisson process. This can be used to
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

   \delta P_d = 1 - \exp{\left( -k_{d} \Delta t \right) } 

where :math:`k_{d}` is the dettachment rate, and :math:`\Delta t` is the timestep. For every 
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

   \delta P_a = 1 - \exp{\left( -k_{a} \Delta t \right) } 

where :math:`k_a` is the attachment rate, and :math:`\Delta t` is the timestep. If the
probability constraint is satisfied, then the bond will be formed. Note that with this
method, each atom may be part of multiple created bonds on a given time step.

Either of the 2 quantities defining the rates can be specified as an equal-style 
:doc:`variable <variable>`, namely *ka*, *kd*. If the value is a variable, it should
be specified as v_name, where name is the variable name. In this case, the variable
will be evaluated each timestep, and its value used to determine the attachment
and dettachment rates respectively.

It is permissible to have *itype* = *jtype*\ .  *Rcut* must be :math:`\leq` the
pair-wise cutoff distance between *itype* and *jtype* atoms, as defined
by the :doc:`pair_style <pair_style>` command.

The *maxbond* keyword can be used to limit the number of bonds allowed. 
If either atom :math:`i` of type *itype* or atom :math:`j` of type *jtype* 
has *maxbond* bonds (set by value), then a new bond will not be formed.
By default *maxbond* is set to the "extra bond per atom" parameter.

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

The *rouse* keyword can be used to modify bonds attachment rate :math:`k_{a}` by
assuming that bonds must 'explore' their surrounding space through a
sub-diffusive Rouse process before attaching. Thus, instead of a fixed 
attachement rate (and thus constant probability) chains rate of attachement 
scales nonlinearly with the distance between two atoms :math:`r` according to:

.. math::

   k_{a}^{rouse} = k_{a} \left( \frac{b0}{r} \right)^4

where :math:`k_{a}` is the nominal or fixed attachment rate and
:math:`b0` is a distance. In the case of polymeric systems, assuming
flexible ergodic chains, :math:`b0` is the molecular distance travelled
in time :math:`1/k_{a}`. The *rouse* keyword cannot be used with *prob*.

The *bell* keyword can be used to modify bonds dettachement rate :math:`k_{a}`
by assuming bonds dettachment kinetics is force-sensitive. In this case, bonds dettachment
rate increases exponentially under increasing force given by Bell's model:

.. math::

   k_d^{bell} = k_{d} \exp{ \left( \frac{f}{f0} \right)}

where :math:`k_{d}` is the nominal or fixed dettachment rate, :math:`f` is
the bonds force, and :math:`f0` characterizes the bonds force-sensitivity. 
The *bell* keyword cannot be used with keyword *prob* or *catch*.

The *catch* keyword can be used to modify bonds dettachement rate :math:`k_{d}`
by assuming bonds dettachment kinetics is force-sensitive. Unlike Bell's model,
bonds dettachment rate at first decreases before subsequently increasing under
increasing force. This is achieved with the two-pathway model:

.. math::

   k_d^{catch} = k_{d} \exp{ \left( \frac{f}{fs0} \right)} +  k_{d} kc0 \exp{ \left( \frac{-f}{fc0} \right)}

where :math:`k_{d}` is the nominal or fixed dettachment rate, :math:`f` is
the bonds force, :math:`fs0` is the force-sensitivity of the slip barrier, 
:math:`fc0` is the force-sensitivity of the catch barrier, and :math:`kc0` is
a scale factor that adjusts the fixed dettachment rate of the catch barrier.
Note that when :math:`kc0` = :math:`0.0` the Bell model is recovered. Also,
when bond forces are small (i.e 0.0) the detachment rate is :math:`k_{d} * (1 + kc0)`.
The *catch* keyword cannot be used with keyword *prob* or *bell*.

Any bond that is created is assigned a bond type of *bondtype*.

When a bond is created, data structures within LAMMPS that store bond
topologies are updated to reflect the creation. 

When a bond is deleted, data structures within LAMMPS that store bond
history variables, such as those used by the BPM package, are updated to
reflect deletion.

.. note::

   One data structure that is not updated when a bond either breaks or is 
   created are the molecule IDs stored by each atom. Even with newly created bonds,
   all atoms in the new "molecule" retain their original molecule IDs.

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

   The list of topological neighbors is updated for atoms
   affected by the new bond.  This in turn affects which neighbors are
   considered for pair-wise interactions, using the weighting rules set by
   the :doc:`special_bonds <special_bonds>` command.  Consider a new bond
   created between atoms :math:`i` and :math:`j`.  If :math:`j` has a bonded
   neighbor :math:`k`, then :math:`k` becomes a second neighbor of :math:`i`. 
   The pair-wise interaction between :math:`i` and :math:`k` could potentially
   be turned off or weighted by the 1--3 weighting specified
   by the :doc:`special_bonds <special_bonds>` command. The same
   is true for third neighbors (1--4 interactions).

Note that even if your simulation starts with no bonds, you must
define a :doc:`bond_style <bond_style>` and use the
:doc:`bond_coeff <bond_coeff>` command to specify coefficients for the
*bondtype*\ .  Similarly, the atom types specified by the
*iype* or *jtype* values must be within the range of atom
types allowed by the simulation.

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

Restrictions
""""""""""""

This fix is part of the TNT package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`fix bond/break <fix_bond_break>`, :doc:`fix bond/react <fix_bond_react>`, :doc:`fix bond/swap <fix_bond_swap>`,
:doc:`dump local <dump>`, :doc:`special_bonds <special_bonds>`

Default
"""""""

The option defaults are seed = 12345, mol = 0, and
maxbond = extra/bond/per/atom value.
