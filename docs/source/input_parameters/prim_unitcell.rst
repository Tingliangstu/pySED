prim_unitcell
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   prim_unitcell = 3.11905 0 0 -1.55952 2.70117 0 0 0 12.6227

**Meaning**
   Primitive-cell lattice vectors in Angstrom, written as nine numbers and
   reshaped into a 3 by 3 matrix.

**Default**
   No useful default. Set it explicitly.

**Notes**
   The rows are lattice vectors. Use the primitive cell consistent with
   ``basis.in`` and ``supercell_dim``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`prim_axis <prim_axis>`
- :doc:`supercell_dim <supercell_dim>`
- :doc:`rescale_prim <rescale_prim>`
