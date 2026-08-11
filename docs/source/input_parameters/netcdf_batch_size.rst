netcdf_batch_size
~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   netcdf_batch_size = 64

**Meaning**
   Number of text-trajectory frames parsed at a time while converting a GPUMD
   XYZ or LAMMPS dump to NetCDF, or while supplying direct text input to the
   common block interface. It is not the NetCDF/HDF5 on-disk chunk size and it
   does not change the SED averaging block length.

**Default**
   ``64``.

**Notes**
   Smaller values reduce temporary parser memory. Larger values may reduce
   Python batching overhead, but do not guarantee faster trajectory reads.
   Native NetCDF input is unaffected.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`num_blocks <num_blocks>`
