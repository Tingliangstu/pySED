trajectory_read_mode
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_read_mode = cache

**Meaning**
   Selects how a GPUMD XYZ or one-file LAMMPS dump is supplied to SED
   computation:

   - ``cache`` converts the text once to
     ``<source-name>.pysed.nc`` beside the input file and reuses that binary
     cache.
   - ``direct`` parses the requested text frames sequentially during SED
     computation and creates no intermediate NetCDF file.

   Direct text reading uses one common block interface for GPUMD and LAMMPS
   and never restarts the text scan at each averaging block. A native NetCDF
   ``trajectory_file`` is always read directly and is unaffected by this
   parameter.

**Default**
   ``cache``.

**Choosing a mode**
   ``direct`` avoids conversion time and cache storage for a one-off
   calculation. ``cache`` is usually better when the same trajectory will be
   recomputed with different Q paths or settings. To retain a cache without
   compression overhead, combine ``cache`` with
   ``netcdf_compression_level = 0``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_prefetch <trajectory_prefetch>`
- :doc:`output_data_stride <output_data_stride>`
