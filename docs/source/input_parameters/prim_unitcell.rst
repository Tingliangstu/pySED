prim_unitcell
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   prim_unitcell = 3.11905 0 0 -1.55952 2.70117 0 0 0 12.6227

**Meaning**
   Primitive-cell lattice vectors in Angstrom. Write the three lattice vectors
   as nine ordinary decimal numbers, and pySED reshapes them into a 3 by 3
   matrix.

**Default**
   No useful default. Set it explicitly.

**Notes**
   The rows are lattice vectors. Use the primitive cell consistent with
   ``basis.in`` and ``supercell_dim``.

   This parameter does not support fractional strings such as ``1/3`` or
   ``1/2``. Use decimal values instead, for example ``0.3333333333`` or
   ``0.5``. Fraction input is supported for ``q_path``, but not for
   ``prim_unitcell``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`prim_axis <prim_axis>`
- :doc:`supercell_dim <supercell_dim>`
- :doc:`rescale_prim <rescale_prim>`
