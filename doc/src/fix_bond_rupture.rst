.. index:: fix bond/rupture

fix bond/rupture command
======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID bond/rupture bondtype style args keyword value

* ID, group-ID are documented in :doc:`fix <fix>` command
* bond/rupture = name of this fix command
* bondtype = type of bonds to break (integer or type label)
* style_name = *dist* or *prob/fraction* or *prob/rate* or *prob/slip* or *prob/slip/catch* 

  .. parsed-literal::
       *dist* args = rcrit
         rcrit = bond longer than rcrit can break if otherwise eligible (distance units)
       *prob/fraction* args = fraction seed
         fraction = break a bond with this probability if otherwise eligible
         seed = random number seed (positive integer)
       *prob/rate* args = kr
         kr = bond rutpure rate if otherwise eligible (inverse time units)
       *prob/slip* args = ks0 f0
         ks0 = intrinsic bond rutpure rate if otherwise eligible (inverse time units)
         f0 = force-sensitivity of the bonds slip barrier (force units)
       *prob/slip/catch* args = ks0 kc0 fs0 fc0
         ks0 = slip barrier intrinsic rate (inverse time units)
         kc0 = catch barrier intrinsic rate (inverse time units)
         fs0 = force-sensitivity of the bonds slip barrier (force units)
         fc0 = force-sensitivity of the bonds catch barrier (force units)

* zero or more keyword/value pairs may be appended
* keyword = *bond/table* or *bond/distribution* or *critical*

  .. parsed-literal::

       *bond/table* values = filename key
         filename = file containing the tablulated per/bond style args
         key = section of filename to start reading
       *bond/distribution* values = dist_type style_arg params
         dist_type = *uniform* or *gauss* or *exponential* or *weibull*
            uniform params = lo hi
               lo = lower bound value of style_arg
               hi = upper bound value of style_arg
            gauss params =  mu sigma
               mu = average value of style_arg
               sigma = standard deviation of style_arg
            exponential params = lambda
               lambda = decay lengthscale of style_arg
            weibull params = alpha beta
               alpha = scale factor of style_arg
               beta = shape factor of style_arg
       *critical* value = rcrit
            rcrit = enforce bonds longer than rcrit rupture (distance units)
   
Examples
""""""""

.. code-block:: LAMMPS

   fix 5 all bond/rupture 1 dist 0.3
   fix 5 all bond/rupture 2 prob/fraction 0.3 12345 bond/table rupture.table start critical 0.3
   fix 5 all bond/rupture 1 prob/slip 0.2 1.0 bond/distribution gauss ks0 0.2 2.0 f0 1.0 1.0 

Description
"""""""""""

Break bonds between pairs of atoms as a simulation runs according to
specified criteria.  This can be used to model rupture of chains 
in a polymer network due to stretching of the simulation box or other
deformations.  In this context, a bond means an interaction between a
pair of atoms computed by the :doc:`bond_style <bond_style>` command.
Once the bond is broken it will be permanently deleted, as will all
angle, dihedral, and improper interactions that bond is part of. There
are several possible styles that determine the nature of the rupture criterion.

The *dist* style specifies rupture after bonds exceed a critial legnth set by value *rcrit*.

The *prob/fraction* style specifies bond rupture based on a fixed probability set 
by the value *fraction*, which must be a value between 0 and 1. For rupture,
a uniform random number between 0.0 and 1.0 is generated and the bond is only
broken if the random number is less than *fraction*. The seed can be used to 
specifiy the processor-unique seed used to initialized the Marsaglia random 
number generator. By default the seed is 12345. The value setting must be a positive integer.

The *prob/rate* style specifices bond rupture based on a constant rupture rate set
by the value *kr*, which must be a positive value. A bond will break with a discrete
rupture probability defined as:

.. math::

   \delta P_r = 1 - \exp{\left( -k_{r} \Delta t \right) } 

where :math:`k_{r}` is the rupture rate, and :math:`\Delta t` is the timestep. For every
eligible bond if the probability constraint is satisfied then the bond is destroyed, 
otherwise it remains.

The *prob/slip* style can be used to modify the rupture probability by assuming that
rupture kinetics is force-senstive. In this case, bonds rupture rate increases 
exponentially under increasing force given by Bell's model:

.. math::

   k_r^{slip} = k_{s0} \exp{ \left( \frac{f}{f0} \right)}

where :math:`k_{s0}` is the nominal or fixed rupture rate in the absence of
force, :math:`f` is the bonds force, and :math:`f0` characterizes 
the bonds force-sensitivity.

The *prob/slip/catch* style can be used to modify the rupture probability by
assuming that rupture kinetics is force-senstive. However, unlike *prob/slip*
the rupture rate at first decreases before subsequently increasing under
increasing force. This is achieved with the two-pathway model:

.. math::

 k_r^{slip-catch} =k_{s0} \exp{ \left( \frac{f}{fs0} \right)} + k_{c0} \exp{ \left( -\frac{f}{fc0} \right)}

where :math:`f` is the bonds force, :math:`k_{s0}` and 
:math:`k_{c0}` are the nominal rupture rates of the slip 
and catch pathways respectively. The force-sensitivity of
the slip and catch barriers are given by :math:`fs0 and :math:`fc0.

The *bond/table* keyword allows a unique rupture criterion to be
defined on a per bond basis, by specifying a tabulated file with 
arguments for each bond. The expected number of arguments in the 
table depends on the rupture style. The filename specifies the file containing 
the tablulated arguments, and the keyword specifies a section of 
the file to begin reading from. An example format of the file 
for the *dist* style is provided below. The *bond/table* keyword cannot
be used with *bond/distribution* keyword

The *bond/distribution* keyword allows a unique rupture criterion to be
defined on a per bond basis, by drawing certain style arguments
from a defined distribution type. To distribute and argument, the
argument is specified followed by the parameters that define the
distribution. For example, with distribution *Gauss* and rupture 
style *dist* each bond will have a value *rcrit* drawn from a normal
distribution with mean :math:`\mu` and standard deviation :math:`\sigma`.
One or many arguments can be drawn, each from unique distributions, although
the distribution type cannot be changed. Note that not all style arguments
can be drawn from the distibution, for example *seed* in the *prob/fraction*
style is set globally. The *bond/distribution* keyword cannot
be used with *bond/table* keyword.

The argument drawn, depends on the
rupture style. For instance with distribution *Gauss* and rupture 
style *dist* the value *rcrit* is drawn with mean :math:`\mu` and 
standard deviation :math:`\sigma`. Note that not all style arguments
can be drawn from the distibution, for example *seed* in the *prob/fraction*
style is set globally. The *bond/distribution* keyword cannot
be used with *bond/table* keyword.

The *critical* keyword enforces that bonds rupture after exceeding a
critical length set by value *rcrit*. This can for example be used 
in conjunction with the *prob* rupture styles, to ensure bond rupture. 

When a bond is broken, data structures within LAMMPS that store bond
topologies are updated to reflect the breakage.  Likewise, if the bond
is part of a 3-body (angle) or 4-body (dihedral, improper)
interaction, that interaction is removed as well.  These changes
typically affect pair-wise interactions between atoms that used to be
part of bonds, angles, etc.

.. note::

   One data structure that is not updated when a bond breaks are
   the molecule IDs stored by each atom.  Even though one molecule
   becomes two molecules due to the broken bond, all atoms in both new
   molecules retain their original molecule IDs.

Computationally, each time step this fix is invoked, it loops over all
the bonds in the system and computes distances between pairs of bonded
atoms.  It also communicates between neighboring processors to
coordinate which bonds are broken.  Moreover, if any bonds are broken,
neighbor lists must be immediately updated on the same time step.  This
is to ensure that any pair-wise interactions that should be turned "on"
due to a bond breaking, because they are no longer excluded by the
presence of the bond and the settings of the
:doc:`special_bonds <special_bonds>` command, will be immediately
recognized.  All of these operations increase the cost of a time step.
Thus, you should be cautious about invoking this fix too frequently.

You can dump out snapshots of the current bond topology via the :doc:`dump local <dump>` command.

.. note::

   Breaking a bond typically alters the energy of a system.  You
   should be careful not to choose bond breaking criteria that induce a
   dramatic change in energy.  For example, if you define a very stiff
   harmonic bond and break it when two atoms are separated by a distance
   far from the equilibrium bond length, then the two atoms will be
   dramatically released when the bond is broken.  More generally, you
   may need to thermostat your system to compensate for energy changes
   resulting from broken bonds (as well as angles, dihedrals, and impropers).

See the :doc:`Howto <Howto_broken_bonds>` page on broken bonds for more
information on related features in LAMMPS.

----------

Formatting the table file
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The format of a tabulated file for the *dist* style is as follows (without parenthesized comments):

.. code-block:: LAMMPS

   # Args for bond rupture (dist style)  (one or more comment or blank lines)
   
   DATA                                  (keyword is the first text on line)
   n 5                                   (n bonds)
                                         (blank line)
   1 1 2 0.20                            (index, iatom, jatom, rcrit)
   2 1 3 0.25
   ...
   5 4 5 0.32

The number of bonds *n* defined in the table file must be equal to 
the total number of bonds in the simulation. The first three column 
entries must be the index, iatom, and jatom (in that order). The index 
is a dummy index as LAMMPS uses iatom jatom indexes to properly store data. 
The following column entries depends on the selected style
and must mactch the that styles number of args. 

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the BPM package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`fix bond/create <fix_bond_create>`, :doc:`fix bond/react <fix_bond_react>`, :doc:`fix bond/swap <fix_bond_swap>`,
:doc:`dump local <dump>`, :doc:`special_bonds <special_bonds>`

Default
"""""""

The option defaults are seed = 12345.
