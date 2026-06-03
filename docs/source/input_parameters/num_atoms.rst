num_atoms
~~~~~~~~~

**Syntax**

.. code-block:: text

   num_atoms = 13824

**Meaning**
   Number of atoms in the MD supercell used for the trajectory. It must match
   both the trajectory file and the atom mapping in ``basis.in``.

**Default**
   ``0``. This must be set for real calculations.

**Example**
   In the MoS\ :sub:`2` GPUMD example, ``num_atoms = 13824``.

**Notes**
   The value must match the maximum atom id in ``basis.in``. pySED stops if
   these two values are inconsistent.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`total_num_steps <total_num_steps>`
- :doc:`time_step <time_step>`
- :doc:`output_data_stride <output_data_stride>`
