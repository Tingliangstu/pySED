prim_axis
~~~~~~~~~

**Syntax**

.. code-block:: text

   prim_axis = 0.0 0.5 0.5 0.5 0.0 0.5 0.5 0.5 0.0

**Meaning**
   Optional transformation from the input cell to a primitive cell, following
   the same idea as phonopy primitive axes. Use it only when the structure file
   is not already expressed in the primitive-cell basis expected by pySED.

**Default**
   ``None``.

**Notes**
   Use only when the input structure cell differs from the primitive cell needed
   for q-point construction. For this workflow, keep the MD cell consistent and
   verify the transformed primitive cell printed by pySED.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`prim_unitcell <prim_unitcell>`
- :doc:`supercell_dim <supercell_dim>`
- :doc:`rescale_prim <rescale_prim>`
