.. index:: bond_style bpm/prony

bond_style bpm/prony command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/prony N keyword value attribute1 attribute2 ...
* N = allocate N history variables for each Mawell element
* optional keyword =  *store/local* or *overlay/pair* or *smooth* or *normalize* or *break* or *plastic* or *nonlinear*

  .. parsed-literal::

       *store/local* values = fix_ID N attributes ...
          * fix_ID = ID of associated internal fix to store data
          * N = prepare data for output every this many timesteps
          * attributes = zero or more of the below attributes may be appended

            *id1, id2* = IDs of two atoms in the bond
            *time* = the timestep the bond broke
            *x, y, z* = the center of mass position of the two atoms when the bond broke (distance units)
            *x/ref, y/ref, z/ref* = the initial center of mass position of the two atoms (distance units)

       *overlay/pair* value = *yes* or *no*
          bonded particles will still interact with pair forces

       *smooth* value = *yes* or *no*
          smooths bond forces near the breaking point

       *normalize* value = *yes* or *no*
          normalizes bond forces by the reference length

       *break* value = *yes* or *no*
          indicates whether bonds break during a run

       *plastic* value = *yes* or *no*
          indicates whether bonds plastically deform

       *nonlinear* value = *yes* or *no*
          indicates whether nonlinear option is used


Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/prony 1
   bond_coeff 1 1.0 0.4 0.1 file.table keyword

   bond_style bpm/prony 1 plastic yes nonlinear yes
   bond_coeff 1 1.0 0.4 0.1 file.table keyword 0.2 2.0

   bond_style bpm/prony 1 myfix 1000 time id1 id2
   dump 1 all local 1000 dump.broken f_myfix[1] f_myfix[2] f_myfix[3]
   dump_modify 1 write_header no

Description
"""""""""""

The *bpm/prony* bond style computes forces based on
deviations from the initial reference state of the two atoms and the strain history. The
reference length :math:`r_0` is stored by each bond when it is first computed in
the setup of a run. Initially, the previous length of the bond :math:`r^{t-1}`
is set equal to :math:`r_0` but evolves during the run. Data is then preserved across
run commands and is written to :doc:`binary restart files <restart>` such that restarting
the system will not reset the reference and previous states of a bond.

This bond style only applies central-body forces which conserve the
translational and rotational degrees of freedom of a bonded set of
particles. The bond force follows a linear viscoelastic formulation based 
on a generalized Maxwell element, as outlined in :ref:`(Groot) <Groot4>`:. 
The bond force is comprised of

.. math::

   F = w (F_{E} + H_D)

where :math:`F_{E}` is the force contribution from the rate-independent
elastic element, and :math:`H_{D}` is the contribution from the rate-dependent 
viscoelastic (Maxwell) elements, and :math:`w` is an optional smoothing factor discussed below.
The elastic bond force has a magnitude of

.. math::

   F_E = k0 (r - r_0)

where :math:`k0` is a stiffness, :math:`r` is the current distance
and :math:`r_0` is the initial distance between the two particles.
The viscoelastic bond force has a magnitude of

.. math::

   H_D = \sum_{j=1}^n h_j^t

where the total viscoelastic bond force is the sum of :math:`j = 1` to :math:`n` 
maxwell elements as defined in the *file.table*. The force contributed by each :math:`j`-th element
at the current timestep is given as

.. math::

   h_j^t = \exp{\frac{-k_j \Delta}{\eta_j}} h_j^{t-1} + \frac{\eta_j}{\Delta t} (1 - \exp{\frac{-k_j \Delta}{\eta_j}} ) [r^{t-1} - r]

where :math:`k_j` is a stiffness, :math:`\eta_j` is a viscosity, :math:`\Delta t` is the timestep,
:math:`r^{t-1}` is the previous bond length, :math:`r` is the current bond length, and :math:`h_j^{t-1}`
is the force contributed on the previous timestep.

Bonds will break at a strain of :math:`\epsilon_c`.  This is done by setting
the bond type to 0 such that forces are no longer computed.

An additional damping force is applied to the bonded
particles.  This forces is proportional to the difference in the
normal velocity of particles using a similar construction as
dissipative particle dynamics :ref:`(Groot) <Groot4>`:

.. math::

   F_D = - \gamma w (\hat{r} \bullet \vec{v})

where :math:`\gamma` is the damping strength, :math:`\hat{r}` is the
radial normal vector, and :math:`\vec{v}` is the velocity difference
between the two particles.

The smoothing factor :math:`w` can be added or removed by setting the
*smooth* keyword to *yes* or *no*, respectively. It is constructed such
that forces smoothly go to zero, avoiding discontinuities, as bonds
approach the critical strain

.. math::

   w = 1.0 - \left( \frac{r - r_0}{r_0 \epsilon_c} \right)^8 .

If the *normalize* keyword is set to *yes*, the elastic bond force will be
normalized by :math:`r_0` such that :math:`k` must be given in force units.

By default, pair forces are not calculated between bonded particles.
Pair forces can alternatively be overlaid on top of bond forces by setting
the *overlay/pair* keyword to *yes*. These settings require specific
:doc:`special_bonds <special_bonds>` settings described in the
restrictions.  Further details can be found in the :doc:`how to <Howto_bpm>`
page on BPMs.

If the *break* keyword is set to *no*, LAMMPS assumes bonds should not break
during a simulation run. This will prevent some unnecessary calculation.
The recommended bond communication distance no longer depends on the value of
:math:`\epsilon_c` (which is ignored) but instead corresponds to the typical
heuristic maximum strain used by typical non-bpm bond styles. Similar behavior
to *break no* can also be attained by setting an arbitrarily high value of
:math:`\epsilon_c`. One cannot use *break no* with *smooth yes*.

The *plastic* keyword toggles whether the elastic element is allowed to plastically
deform as done by :doc:`bpm/spring/plastic <bond_bpm+spring_plastic>`. If set to *yes* the elastic
force has a magnitude of

.. math::
   F_{E} = k0 (r - r_{eq})

where :math:`r_{eq}` is the equlibrium bond length.
If the bond stretches beyond a strain of :math:`\epsilon_p` in compression or extension, 
it will plastically activate and :math:`r_{eq}` will evolve to ensure :math:`|(r-r_{eq})/r_{eq}|`
never exceeds :math:`r_{eq}`. Therefore, if a bond is continually loaded in either tension or compression, 
the force in the elastic element will initially grow elastically before plateauing.

The *nonlinear* keyword toggles whether the force in the elastic element is nonlinear. 
If set to *yes* the elastic force has a magnitude of

.. math::
   F_{E} = k0 (r - r_0)^{\alpha}

where :math:`\alpha` is an exponent chosen to model an arbitrary nonlinear response.
Note that the units of :math:`k_0` will depend on :math:`\alpha` in order for consistent force units.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands:

* :math:`k0`             (force/distance units)
* :math:`\epsilon_c`    (unitless)
* :math:`\gamma`        (force/velocity units)
* filename
* keyword

The filename specifies a file containing the tablulated coefficients for the Maxwell 
elements. The keyword specifies a section of the file. The format of this file is described below.
Additionally, if either *plastic* or *nonlinear* are set to *yes*, a sixth or seventh coefficient
must be provided:

* :math:`\epsilon_p`       (unitless)
* :math:`\alpha`           (unitless)

If the *store/local* keyword is used, an internal fix will track bonds that
break during the simulation. Whenever a bond breaks, data is processed
and transferred to an internal fix labeled *fix_ID*. This allows the
local data to be accessed by other LAMMPS commands. Following this optional
keyword, a list of one or more attributes is specified.  These include the
IDs of the two atoms in the bond. The other attributes for the two atoms
include the timestep during which the bond broke and the current/initial
center of mass position of the two atoms.

Data is continuously accumulated over intervals of *N*
timesteps. At the end of each interval, all of the saved accumulated
data is deleted to make room for new data. Individual datum may
therefore persist anywhere between *1* to *N* timesteps depending on
when they are saved. This data can be accessed using the *fix_ID* and a
:doc:`dump local <dump>` command. To ensure all data is output,
the dump frequency should correspond to the same interval of *N*
timesteps. A dump frequency of an integer multiple of *N* can be used
to regularly output a sample of the accumulated data.

Note that when unbroken bonds are dumped to a file via the
:doc:`dump local <dump>` command, bonds with type 0 (broken bonds)
are not included.
The :doc:`delete_bonds <delete_bonds>` command can also be used to
query the status of broken bonds or permanently delete them, e.g.:

.. code-block:: LAMMPS

   delete_bonds all stats
   delete_bonds all bond 0 remove

----------

The format of a tabulated file is as follows (without parenthesized comments):

.. code-block:: LAMMPS

   # Coefficients for Mawell elements   (one or more comment or blank lines)
   
   MAXWELL                              (keyword is the first text on line)
   n 5                                  (n parameters)
                                        (blank line)
   1 1.0 0.1                            (index, stiffness, viscosity)
   2 2.0 100
   ...
   5 0.5 1.0

The number of parameters *n* defined in the table file must be less than or 
equal to the number of entries *N* allocated via the :doc:`bond_style <bond_style>` command.
Therefore, if each bond type uses a unique tabulated file, *N* 
should be allocated according to the largest number of tabulated entries.

----------

.. versionadded:: 25June2025

Restart and other info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This bond style writes the reference state of each bond to
:doc:`binary restart files <restart>`. Loading a restart
file will properly restore bonds. However, the reference state is NOT
written to data files. Therefore reading a data file will not
restore bonds and will cause their reference states to be redefined.

If the *store/local* option is used, an internal fix will calculate
a local vector or local array depending on the number of input values.
The length of the vector or number of rows in the array is the number
of recorded, broken bonds.  If a single input is specified, a local
vector is produced. If two or more inputs are specified, a local array
is produced where the number of columns = the number of inputs.  The
vector or array can be accessed by any command that uses local values
from a compute as input. See the :doc:`Howto output <Howto_output>` page
for an overview of LAMMPS output options.

The vector or array will be floating point values that correspond to
the specified attribute.

Any settings with the *store/local* option are not saved to a restart
file and must be redefined.

The single() function of this bond style returns 0.0 for the energy of a 
bonded interaction, since energy is not conserved in these dissipative potentials. 
However, the single() function also calculates 4 additional quantities. The first 2 pertain 
to bond lengths, including the reference state :math:`r_0` and equlibrium state :math:`r_{eq}`
if the *plastic* option is utilized. If *plastic* = *no* then the equlibrium state 
:math:`r_{eq}` will equal the reference state :math:`r_0`.
The next 2 quantites (3-4) are the split elastic :math:`F_E`
and viscoelastic :math:`H_D` forces respectively.

These extra quantity can be accessed by the
:doc:`compute bond/local <compute_bond_local>` command as *b1*, *b2*, ..., *b4* \.

Restrictions
""""""""""""

This bond style is part of the BPM package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

By default if pair interactions between bonded atoms are to be disabled,
this bond style requires setting

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

and :doc:`newton <newton>` must be set to bond off.  If the *overlay/pair*
keyword is set to *yes*, this bond style alternatively requires setting

.. code-block:: LAMMPS

   special_bonds lj/coul 1 1 1

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`pair bpm/spring <pair_bpm_spring>`

Default
"""""""

The option defaults are *overlay/pair* = *no*, *smooth* = *yes*, *normalize* = *no*, *break* = *yes*, *plastic* = *no*, and *nonlinear* = *no*

----------

.. _fragment-Clemmer:

**(Clemmer)** Clemmer and Robbins, Phys. Rev. Lett. (2022).

.. _Groot4:

**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).

.. _multibody-Clemmer:

**(Clemmer2)** Clemmer, Monti, Lechman, Soft Matter, 20, 1702 (2024).
