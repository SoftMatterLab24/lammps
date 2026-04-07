.. index:: bond_style poly/dFJC

bond_style poly/dFJC command
============================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style poly/dFJC keyword value ...

* optional keyword = *overlay/pair* or *smooth* or *break* or *temp/shift*

  .. parsed-literal::

       *overlay/pair* value = *yes* or *no*
          bonded particles will still interact with pair forces

       *smooth* value = *yes* or *no*
          accepted by the command; the current implementation does not apply a separate smoothing factor in the main force calculation

       *break* value = *yes* or *no*
          indicates whether bonds break during a run

       *temp/shift* value = *yes* or *no*
          enables the optional shift factor :math:`a_T` in :doc:`bond_coeff <bond_coeff>`


Examples
""""""""

.. code-block:: LAMMPS

   bond_style poly/dFJC break yes overlay/pair yes temp/shift yes
   bond_coeff 1 100.0 50.0 0.20 0.0 1.0 bond.table CHAINS temp 1.0

   bond_style poly/dFJC break yes overlay/pair yes temp/shift yes
   bond_coeff 1 100.0 50.0 0.20 0.0 1.0 bond.table CHAINS temp 1.0 stretch 1.15

   compute 2 all bond/local dist force b1 b2 b3 b4 b5

Description
"""""""""""

The *poly/dFJC* bond style computes a central force for a dissipative,
extensible freely jointed chain. In this model the total bond extension
is split into two parts: a segmental extension of the Kuhn segments and
a conformational extension of the freely jointed chain. The per-bond
values of :math:`N` and :math:`b` are read from a table file, where
:math:`N` is the number of Kuhn segments and :math:`b` is the Kuhn
length.

The code stores internal history variables for each bond and updates
them during the run. The total bond length :math:`r` is decomposed as

.. math::

   r = r_s + r_j

where the segmental contribution is

.. math::

   r_s = \frac{N f}{K_s}

and the conformational FJC contribution is

.. math::

   r_j = N b y

with :math:`f` the bond force, :math:`K_s` the segmental stiffness, and
:math:`y` the internal conformational stretch variable of the FJC part.
The FJC contribution uses a Cohen-Pade approximation to the inverse
Langevin response through the stretch-dependent tangent stiffness

.. math::

   K_j^{\mathrm{eff}}(y) = \frac{K_j}{N b^2}\frac{3-y^2}{1-y^2}

where :math:`K_j` sets the entropic stiffness scale. As
:math:`y \rightarrow 1`, the conformational part stiffens strongly.

Unlike an inextensible FJC model, the Kuhn segments themselves may
stretch in *poly/dFJC*. The segmental stretch can therefore exceed 1,
with

.. math::

   \lambda_s = 1 + \frac{f}{K_s b}

This is the quantity used by the optional *stretch* rupture criterion.

The dissipative part of the model is controlled by the segmental
friction coefficient :math:`\zeta`. For each bond, the implementation
uses

.. math::

   \eta = 2 \zeta N

and, when the *temp/shift* keyword is enabled,

.. math::

   \eta_{\mathrm{eff}} = a_T \eta

where :math:`a_T` is a type-specific shift factor. The code updates the
bond force by a fast predictor-corrector solve for modest
conformational stretch and switches to an iterative solve when the
internal conformational stretch becomes larger.

An additional damping force is applied in the bond direction,
proportional to the relative normal velocity of the bonded atoms:

.. math::

   F_d = - \gamma (\hat{r} \bullet \vec{v})

where :math:`\gamma` is the damping strength, :math:`\hat{r}` is the
bond unit vector, and :math:`\vec{v}` is the relative velocity.

By default, rupture is controlled by the critical segmental extension
:math:`r_c`. The code checks

.. math::

   \left|\frac{f}{K_s}\right| > r_c

and breaks the bond when that condition is met.

If the optional *stretch* argument is added to the
:doc:`bond_coeff <bond_coeff>` command, rupture becomes stretch-based
instead. In that case the segmental stretch :math:`\lambda_s` defined
above is compared to a critical value :math:`\lambda_c`, and the bond
breaks when :math:`\lambda_s > \lambda_c`. Because the Kuhn segments
are extensible in this model, the code requires :math:`\lambda_c > 1`.

By default, pair forces are not calculated between bonded particles.
Pair forces can alternatively be overlaid on top of bond forces by
setting the *overlay/pair* keyword to *yes*. These settings require
specific :doc:`special_bonds <special_bonds>` settings as described in
the restrictions.

If the *break* keyword is set to *no*, LAMMPS assumes the bonds do not
rupture during the run. This removes the break checks and also changes
the recommended bond communication distance so that it no longer depends
on the bond rupture limit. As in the code, *smooth yes* cannot be used
together with *break no*.

The *smooth* keyword is accepted for consistency with related *poly*
bond styles, but the current *poly/dFJC* implementation does not apply
a separate smoothing factor in the main compute force path.

If *temp/shift yes* is used, the optional *temp* value in
:doc:`bond_coeff <bond_coeff>` sets the current shift factor
:math:`a_T`. The style also exposes :math:`a_T` through its bond-style
extract interface as *aT*, so it can be modified by commands such as
:doc:`fix adapt <fix_adapt>`.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands:

* :math:`K_s`        (segmental stiffness, force/distance units)
* :math:`K_j`        (entropic stiffness scale, energy or force-distance units)
* :math:`r_c`        (critical segmental extension for the default rupture criterion)
* :math:`\gamma`     (force/velocity units)
* :math:`\zeta`      (segmental friction coefficient, force-time units)
* filename
* keyword

Optional trailing arguments are:

* *temp* :math:`a_T`
* *stretch* :math:`\lambda_c`

The *temp* argument is only recognized when *temp/shift yes* is set in
the :doc:`bond_style <bond_style>` command. If *stretch* is specified,
the rupture criterion switches from the default critical extension
:math:`r_c` to the extensible-segment stretch criterion using
:math:`\lambda_c > 1`.

The filename specifies a table containing the per-bond values of
:math:`N` and :math:`b`. The keyword selects a named section in that
file. Each entry is matched against a bonded atom pair by atom ID.

----------

Formatting the table file
"""""""""""""""""""""""""""""""""

The table file assigns :math:`N` and :math:`b` to specific bonded atom
pairs. Its format is as follows (without parenthesized comments):

.. code-block:: LAMMPS

   # dFJC chain parameters for specific bonds

   CHAINS                              (keyword is the first text on line)
   N 3                                 (number of table entries)

   1 10 11 40 0.50                     (index, atom I, atom J, N, b)
   2 11 12 40 0.50
   3 12 13 60 0.40

The integer atom IDs on each line must match the IDs of a bonded atom
pair. The order of the two IDs is irrelevant.

----------

Restart and other info
"""""""""""""""""""""""""""""""""""""""""""""

This bond style writes its current bond history, type coefficients, and
table data to :doc:`binary restart files <restart>`. Loading a restart
file restores the saved bond state for continued simulation. However,
the internal history is not written to data files, so reading a data
file reinitializes the bond state.

This is a dissipative bond style. The main force computation tallies
zero bond energy, and the single() function also returns 0.0 for the
energy.

The single() function provides five extra per-bond quantities through
:doc:`compute bond/local <compute_bond_local>` as *b1* through *b5*.
For this style they are:

* *b1* = stored :math:`N`
* *b2* = stored :math:`b`
* *b3* = current segmental extension :math:`r_s`
* *b4* = current internal conformational stretch variable :math:`y`
* *b5* = current shift factor :math:`a_T`

For example:

.. code-block:: LAMMPS

   compute 2 all bond/local dist force b1 b2 b3 b4 b5

Here *dist* and *force* are the standard bond-local outputs, while the
five *b* values are the extra quantities supplied by *poly/dFJC*.

Restrictions
""""""""""""

This bond style is part of the POLY-NET package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

By default, if pair interactions between bonded atoms are to be
disabled, this bond style requires setting

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

and :doc:`newton <newton>` must be set to bond off. If the
*overlay/pair* keyword is set to *yes*, this bond style instead
requires

.. code-block:: LAMMPS

   special_bonds lj/coul 1 1 1

As with the other *poly* bond styles, breaking bonds is not compatible
with simulations that also use angles, dihedrals, impropers, or atom
style template.

This bond style requires ghost atom velocities.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`compute bond/local <compute_bond_local>`,
:doc:`fix adapt <fix_adapt>`, :doc:`bond poly/FJC <bond_poly_FJC>`,
:doc:`bond poly/uFJC <bond_poly_uFJC>`

Default
"""""""

The keyword defaults are *break* = *yes*, *overlay/pair* = *no*,
*smooth* = *yes*, and *temp/shift* = *no*. If the optional *stretch*
argument is not given in :doc:`bond_coeff <bond_coeff>`, rupture uses
the default critical segmental extension :math:`r_c` criterion.
