q_path
~~~~~~

**Syntax**

.. code-block:: text

   q_path = 0.0 0.0 0.0  0.5 0.0 0.0  0.3333333 0.3333333 0.0  0.0 0.0 0.0

**Meaning**
   High-symmetry q-points in reduced reciprocal coordinates.

**Default**
   No useful default. Set it explicitly.

**Notes**
   Provide ``num_qpaths + 1`` triples. Fractions such as ``1/3`` are accepted.
   pySED keeps only q-points commensurate with the supercell.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`num_qpaths <num_qpaths>`
- :doc:`q_path_name <q_path_name>`
