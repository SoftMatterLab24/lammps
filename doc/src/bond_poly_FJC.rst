.. index:: bond_style poly/FJC

bond_style poly/FJC command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style poly/FJC keyword value ...

* optional keyword = *overlay/pair* or *smooth* or *break*

  .. parsed-literal::

       *overlay/pair* value = *yes* or *no*
          bonded particles will still interact with pair forces

       *smooth* value = *yes* or *no*
          toggles the smoothing option used with bond breaking

       *break* value = *yes* or *no*
          indicates whether bonds break during a run

Examples
""""""""

.. code-block:: LAMMPS

   bond_style poly/FJC
   bond_coeff 1 1.0 10.0 0.1 bond.table KEY

   bond_style poly/FJC break yes overlay/pair yes smooth no
   bond_coeff 1 $(v_kBT) $(100*v_kBT) 0 bond.table KEY

   bond_style poly/FJC break yes overlay/pair yes smooth no
   bond_coeff 1 $(v_kBT) $(100*v_kBT) 0 bond.table KEY stretch $(v_lamc)
   compute 2 all bond/local dist force b1 b2

Description
"""""""""""

The *poly/FJC* bond style computes a central force for a finitely
extensible freely jointed polymer chain using a Cohen-Pade
approximation to the inverse Langevin function.  The initial bond
length :math:`r_0` is stored when the bond is first initialized during
setup of a run and is then preserved across subsequent run commands and
written to :doc:`binary restart files <restart>`.

In addition to the stored reference length, each bond also stores two
chain parameters, :math:`N` and :math:`b`, that are read from a table
file for each bonded atom pair.  Here :math:`N` is the number of Kuhn
segments and :math:`b` is the Kuhn length.  The contour length of the
chain is then :math:`N b`.

The scalar bond force implemented by this style is

.. math::

   F = -\frac{k_0}{b}\frac{\lambda (3 - \lambda^2)}{1 - \lambda^2}

where

.. math::

   \lambda = \frac{r}{N b}

is the chain stretch ratio and :math:`r` is the current bond length.
As :math:`\lambda \rightarrow 1`, the force diverges, reflecting the
finite extensibility of the chain.

An additional dissipative contribution is added in the bond direction,
proportional to the relative normal velocity of the bonded atoms:

.. math::

   F_d = - \gamma (\hat{r} \bullet \vec{v})

where :math:`\gamma` is the damping strength, :math:`\hat{r}` is the
bond unit vector, and :math:`\vec{v}` is the relative velocity.

By default, bonds break when the magnitude of the FJC bond force exceeds
the critical force :math:`f_c`.  This is done by setting the bond type
to 0 so that the bond no longer contributes forces.

Alternatively, a stretch-based rupture criterion can be selected by
adding the optional *stretch* argument to the :doc:`bond_coeff
<bond_coeff>` command.  If *stretch* is present, rupture switches from
the default force-based criterion using :math:`f_c` to a stretch-based
criterion, and the bond breaks when the chain stretch ratio exceeds the
specified critical value :math:`\lambda_c` with
:math:`0 < \lambda_c < 1`.

The *smooth* keyword controls the smoothing option associated with bond
failure.  Since smoothing is only meaningful when bonds may fail,
*smooth yes* cannot be used together with *break no*.

By default, pair forces are not calculated between bonded particles.
Pair forces can alternatively be overlaid on top of bond forces by
setting the *overlay/pair* keyword to *yes*.  These settings require
specific :doc:`special_bonds <special_bonds>` settings as described in
the restrictions.

If the *break* keyword is set to *no*, LAMMPS assumes the bonds do not
rupture during the run.  This removes the break checks and also changes
the recommended bond communication distance so that it no longer depends
on a maximum chain stretch.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands:

* :math:`k_0`        (energy or force-distance units, usually presented in units of :math:`k_bT`)
* :math:`f_c`        (critical bond force for the default rupture criterion)
* :math:`\gamma`     (force/velocity units)
* filename
* keyword

An optional rupture criterion may be added after those required values:

* *stretch* :math:`\lambda_c`

If *stretch* is specified, :math:`f_c` is still read by the command but
bond breaking is determined by :math:`\lambda_c` instead of the force
threshold.

The filename specifies a table containing the per-bond values of
:math:`N` and :math:`b`.  The keyword selects a named section in that
file.  Each entry is matched against a bonded atom pair by atom ID.

----------

Formatting the table file
"""""""""""""""""""""""""""""""""

The table file assigns :math:`N` and :math:`b` to specific bonded atom
pairs.  Its format is as follows (without parenthesized comments):

.. code-block:: LAMMPS

   # FJC chain parameters for specific bonds

   CHAINS                              (keyword is the first text on line)
   N 3                                 (number of table entries)

   1 10 11 40 0.5                      (index, atom I, atom J, N, b)
   2 11 12 40 0.5
   3 12 13 60 0.4

The integer atom IDs on each line must match the IDs of a bonded atom
pair.  The order of the two IDs is irrelevant.

----------

Restart and other info
"""""""""""""""""""""""""""""""""""""""""""""

This bond style writes the reference state of each bond, the stored
table data, and the style settings to :doc:`binary restart files
<restart>`.  Loading a restart file will therefore properly restore the
bond model.  However, the per-bond reference state is not written to
data files, so reading a data file will reinitialize the bond state.

The single() function returns the free-energy expression used by this
style,

.. math::

   E = N\left(\frac{\lambda^2}{2} - \ln(1 - \lambda^2)\right)

It also provides two extra per-bond quantities through
:doc:`compute bond/local <compute_bond_local>` as *b1* and *b2*.  For
this style, these are the two extra values returned by the bond's
single() method and correspond to the stored chain parameters
:math:`N` and :math:`b`, respectively.  For example:

.. code-block:: LAMMPS

   compute 2 all bond/local dist force b1 b2

Here *dist* and *force* are the standard bond-local outputs, while *b1*
and *b2* are the extra values supplied by *poly/FJC*.

This bond style does not support dynamic bond creation after the bond
state has been initialized.

Restrictions
""""""""""""

This bond style is part of the POLY-NET package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

By default, if pair interactions between bonded atoms are to be
disabled, this bond style requires setting

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

and :doc:`newton <newton>` must be set to bond off.  If the
*overlay/pair* keyword is set to *yes*, this bond style instead
requires

.. code-block:: LAMMPS

   special_bonds lj/coul 1 1 1

As with the other *poly* bond styles, breaking bonds is not compatible
with simulations that also use angles, dihedrals, impropers, or atom
style template.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`compute bond/local <compute_bond_local>`,
:doc:`bond poly/uFJC <bond_poly_uFJC>`, :doc:`bond poly/dFJC <bond_poly_dFJC>`,
:doc:`fix bond/rupture <fix_bond_rupture>`

Default
"""""""

The keyword defaults are *break* = *yes*, *overlay/pair* = *no*, and
*smooth* = *yes*.  Bond rupture is force-based by default through
:math:`f_c`; if the optional *stretch* argument is used in
:doc:`bond_coeff <bond_coeff>`, rupture becomes stretch-based instead.
