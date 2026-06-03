compress
~~~~~~~~

**Syntax**

.. code-block:: text

   compress = 1

**Meaning**
   Compress trajectory coordinates and velocities into an HDF5 file before SED
   computation.

**Default**
   ``1``.

**Notes**
   Keep this enabled for normal workflows. If the HDF5 file already exists,
   pySED reuses it.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`num_splits <num_splits>`
