.. index:: bond_style poly/uFJC

bond_style poly/uFJC command
============================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style poly/uFJC keyword value ...

* optional keyword = *overlay/pair* or *smooth* or *break*

  .. parsed-literal::

       *overlay/pair* value = *yes* or *no*
          bonded particles will still interact with pair forces

       *smooth* value = *yes* or *no*
          accepted by the command; the current implementation only applies it in the bond single() path, not in the main compute force calculation

       *break* value = *yes* or *no*
          indicates whether bonds break during a run

Examples
""""""""

.. code-block:: LAMMPS

   bond_style poly/uFJC
   bond_coeff 1 1.0 100.0 10.0 0.1 bond.table CHAINS

   bond_style poly/uFJC break yes overlay/pair yes smooth no
   bond_coeff 1 1.0 100.0 10.0 0.1 bond.table CHAINS stretch 1.05

   bond_style poly/uFJC break yes overlay/pair yes smooth no
   bond_coeff 1 1.0 100.0 10.0 0.1 bond.table CHAINS bond/distribution 1.6 exponential 30 27 35

   compute 2 all bond/local dist force b1 b2 b3 b4

Description
"""""""""""

The *poly/uFJC* bond style computes a central force for an extensible
freely jointed chain, which considers both entropic elasticity
due to chain conformational changes, and the enthalpic stretching of 
the Kuhn segments. For each bonded atom pair, the style reads the
chain parameters :math:`N` and :math:`b` from a table file, where
:math:`N` is the number of Kuhn segments and :math:`b` is the Kuhn
length. The contour length of the chain is then :math:`N b`.

The implementation also stores the initial bond length :math:`r_0` when
the bond is first initialized. That stored reference state is preserved
across subsequent run commands and written to :doc:`binary restart files
<restart>`. The constitutive force law itself is based on the current
bond length :math:`r` together with the per-bond chain parameters
:math:`N` and :math:`b`.

The chain stretch is defined as

.. math::

   \lambda = \frac{r}{N b}

with the segmental stretch computed using the Bergstrom approximant and 
following quadratic root analysis

.. math::

   \lambda_v = \frac{\lambda + 1 + \sqrt{(\lambda - 1)^2 + 4/\kappa}}{2}

where :math:`\kappa` is the nondimensional segment stiffness for the Kuhn
segments. The entropic stretch variable used by the force law is then

.. math::

   y = \lambda - \lambda_v + 1

and the corresponding inverse-Langevin factor is

.. math::

   \xi = \frac{y(3-y^2)}{1-y^2}

The scalar bond force implemented by this style is

.. math::

   F = -\frac{k_0}{b} \xi

so that explicitly

.. math::

   F = -\frac{k_0}{b} \frac{y(3-y^2)}{1-y^2}

An additional dissipative contribution is added in the bond direction,
proportional to the relative normal velocity of the bonded atoms:

.. math::

   F_d = - \gamma (\hat{r} \bullet \vec{v})

where :math:`\gamma` is the damping strength, :math:`\hat{r}` is the
bond unit vector, and :math:`\vec{v}` is the relative velocity.

By default, rupture is controlled by the force threshold
:math:`f_c`. The code breaks the bond when the magnitude of the uFJC
bond force exceeds :math:`f_c`. This is done by setting the bond type
to 0 so that the bond no longer contributes forces.

Alternatively, a stretch-based rupture criterion can be selected by
adding the optional *stretch* argument to the :doc:`bond_coeff
<bond_coeff>` command. If *stretch* is present, rupture switches from
the default force-based criterion to a segmental-stretch criterion, and
the bond breaks when :math:`\lambda_v > \lambda_c`. In the current
implementation the command only requires :math:`\lambda_c > 0`. Since
:math:`\lambda_v` is the segmental stretch, choosing
:math:`\lambda_c > 1` corresponds to rupture after the Kuhn segments
have stretched beyond their undeformed length.

The *smooth* keyword is accepted for consistency with related *poly*
bond styles, but the current *poly/uFJC* implementation does not apply
a separate smoothing factor in the main compute force path. The current
code only applies smoothing in the bond single() path used for
single-bond evaluation.

By default, pair forces are not calculated between bonded particles.
Pair forces can alternatively be overlaid on top of bond forces by
setting the *overlay/pair* keyword to *yes*. These settings require
specific :doc:`special_bonds <special_bonds>` settings as described in
the restrictions.

If the *break* keyword is set to *no*, LAMMPS assumes the bonds do not
rupture during the run. This removes the break checks and also changes
the recommended bond communication distance so that it no longer depends
on the rupture criterion. As in the code, *smooth yes* cannot be used
together with *break no*.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands:

* :math:`k_0`        (energy or force-distance units)
* :math:`\kappa`     (nondimensional segment stiffness)
* :math:`f_c`        (critical bond force for the default rupture criterion)
* :math:`\gamma`     (force/velocity units)
* filename
* keyword

Zero or more additional settings/keywords may be appended after those required values:

* keyword = *stretch* or *bond/distribution* or *seed*

  .. parsed-literal::

      *stretch* values = lambda_c
            lambda_c = enforce bonds with stretch greater than lambda_c rupture
      *bond/distribution* values = b dist_type params
         b =  Kuhn segment length
         dist_type = *uniform* or *gauss* or *exponential* or *weibull*
            uniform params = lo hi
               lo = lower bound value of style_arg
               hi = upper bound value of style_arg
            gauss params =  mu sigma
               mu = average value of style_arg
               sigma = standard deviation of style_arg
            exponential params = lambda
               lambda = mean of style_arg
               l_min = lower bound value of style_arg
               l_max = upper bound value of style_arg
            weibull params = alpha beta
               alpha = scale factor of style_arg
               beta = shape factor of style_arg
      *seed* value = 
            seed = random number seed (positive integer)


If *stretch* is specified, bond breaking is determined by 
:math:`\lambda_c` through the segmental stretch :math:`\lambda_v` 
instead of the force threshold.

If *bond/distribution* keyword is specified newly created bonds
state variable :math:`N` is assigned by drawing the value from a 
defined distribution type. The state variable :math:`b` is currently not 
supported and is instead set as a fixed value. Note that this setting 
is required if bonds are added dynamically during a simulation run
as achieved by :doc:`fix bond/create <fix_bond_dynamic>` and
:doc:`fix bond/create <fix_bond_create>` for example. To distribute 
:math:`N`, the distribution type is specified followed by the 
parameters that define the distribution. For example, 
with distribution *Gauss* each new bond will have a value :math:`N` 
drawn from a normal distribution with mean 
:math:`\mu` and standard deviation :math:`\sigma`.

The *seed* keyword  can be used to specifiy the processor-unique seed 
used to initialized the Marsaglia random number generator. By default 
the seed is 12345. The value setting must be a positive integer.

The filename specifies a table containing the per-bond values of
:math:`N` and :math:`b`. The keyword selects a named section in that
file. Each entry is matched against a bonded atom pair by atom ID.

----------

Formatting the table file
""""""""""""""""""""""""""""""

The table file assigns :math:`N` and :math:`b` to specific bonded atom
pairs. Its format is as follows (without parenthesized comments):

.. code-block:: LAMMPS

   # uFJC chain parameters for specific bonds

   CHAINS                              (keyword is the first text on line)
   N 3                                 (number of table entries)

   1 10 11 40 0.50                     (index, atom I, atom J, N, b)
   2 11 12 40 0.50
   3 12 13 60 0.40

The integer atom IDs on each line must match the IDs of a bonded atom
pair. The order of the two IDs is irrelevant.

----------

Restart and other info
""""""""""""""""""""""""""""""""

This bond style writes its stored reference state, table data, and
style settings to :doc:`binary restart files <restart>`. Loading a
restart file will therefore restore the saved bond model. However, the
reference state is not written to data files, so reading a data file
reinitializes the bond state.

The single() function returns 0.0 for the energy. It also provides four
extra per-bond quantities through :doc:`compute bond/local
<compute_bond_local>` as *b1* through *b4*. For this style they are:

* *b1* = stored :math:`N`
* *b2* = stored :math:`b`
* *b3* = current segmental stretch :math:`\lambda_v`
* *b4* = current inverse-Langevin factor :math:`\xi`

For example:

.. code-block:: LAMMPS

   compute 2 all bond/local dist force b1 b2 b3 b4

Here *dist* and *force* are the standard bond-local outputs, while the
four *b* values are the extra quantities supplied by *poly/uFJC*.

This bond style does not support dynamic bond creation after the bond
state has been initialized.

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
:doc:`bond poly/FJC <bond_poly_FJC>`, :doc:`bond poly/dFJC <bond_poly_dFJC>`

Default
"""""""

The keyword defaults are *break* = *yes*, *overlay/pair* = *no*, and
*smooth* = *yes*. If the optional *stretch* argument is not given in
:doc:`bond_coeff <bond_coeff>`, rupture uses the default force-based
criterion through :math:`f_c`.
