rescale_prim
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   rescale_prim = 1

**Meaning**
   Controls whether pySED reconstructs the primitive cell from the actual
   trajectory cell. This is useful after NPT relaxation, where the final cell
   may differ slightly from the ideal ``prim_unitcell * supercell_dim``.

**Default**
   ``1``.

**Notes**
   This is useful after NPT relaxation. For strict NVT workflows, the trajectory
   cell should normally match the expected cell.

   Set this parameter to ``0`` with
   :doc:`dump_group_supercell <dump_group_supercell>`. In that workflow the trajectory
   may contain the full simulation cell while ``basis.in`` maps only a local
   group, so reconstructing a primitive cell from the trajectory cell would be
   physically inconsistent.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`prim_unitcell <prim_unitcell>`
- :doc:`prim_axis <prim_axis>`
- :doc:`supercell_dim <supercell_dim>`
- :doc:`dump_group_supercell <dump_group_supercell>`
