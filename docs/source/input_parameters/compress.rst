compress
~~~~~~~~

**Syntax**

.. code-block:: text

   compress = 1

**Meaning**
   Controls whether pySED first stores the trajectory coordinates and
   velocities in an HDF5 file before calculating SED. Keeping this enabled is
   usually recommended because it avoids repeatedly parsing large trajectory
   files.

**Default**
   ``1``.

**Notes**
   Keep this enabled for normal workflows. If the HDF5 file already exists,
   pySED reuses it.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`num_splits <num_splits>`
