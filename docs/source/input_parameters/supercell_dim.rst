supercell_dim
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   supercell_dim = 12 12 16

**Meaning**
   Number of primitive-cell repeats used to build the MD supercell. This
   should match the actual MD trajectory cell, because pySED uses it to decide
   which q-points are commensurate with the simulation.

**Default**
   ``1 1 1``.

**Notes**
   This is one of the most important parameters in ``input_SED.in`` because it
   determines the repetition matrix used for commensurate q-point selection.
   pySED does not use a separate ``num_qpoints`` parameter; the available
   q-points are determined by ``supercell_dim``, the primitive cell, the
   trajectory cell, and ``q_path``.

   For a simple diagonal supercell with repetition matrix
   :math:`P=\mathrm{diag}(N_x,N_y,N_z)`, a path from Gamma to the zone boundary
   along one reduced reciprocal direction, for example
   ``q_path = 0 0 0  0 0 0.5``, gives allowed points
   :math:`0, 1/N_z, 2/N_z, \ldots, \lfloor N_z/2 \rfloor/N_z`. The number of
   q-points is therefore :math:`\lfloor N_z/2 \rfloor + 1`, or
   :math:`N_z/2 + 1` when :math:`N_z` is even.

   A non-orthogonal real-space cell does not by itself change this counting if
   the supercell is still a diagonal repetition of the primitive lattice in
   reduced coordinates. What matters is the integer repetition matrix
   :math:`P` satisfying ``supercell = P @ primitive``.

   For a non-diagonal or transformed supercell, the number of q-points depends
   on the selected path direction. pySED keeps only fractional positions
   :math:`f` on the line
   :math:`q(f)=q_{\mathrm{start}}+f(q_{\mathrm{end}}-q_{\mathrm{start}})`
   that satisfy :math:`q(f)P^T \in \mathbb{Z}^3`. In this case, inspect the
   printed ``Number of q-points generated`` message or the output ``.Qpts``
   file. **Increase the supercell repeats in the lattice directions that
   project onto the desired q-path to improve resolution.** If the calculation
   can still fit in memory, use as large a supercell as practical; a larger
   supercell gives denser q-point sampling.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`prim_unitcell <prim_unitcell>`
- :doc:`basis_lattice_file <basis_lattice_file>`
- :doc:`num_qpaths <num_qpaths>`
- :doc:`q_path_name <q_path_name>`
- :doc:`q_path <q_path>`
- :doc:`prim_axis <prim_axis>`
