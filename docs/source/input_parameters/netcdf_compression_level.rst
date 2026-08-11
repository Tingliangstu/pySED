netcdf_compression_level
~~~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   netcdf_compression_level = 1

**Meaning**
   Sets the compression level used when a GPUMD XYZ or LAMMPS dump is
   converted to a reusable ``.pysed.nc`` cache. Accepted values are ``0``
   through ``9``:

   - ``0`` disables compression;
   - ``1`` enables light compression and is the default;
   - ``2`` through ``9`` request progressively stronger compression.

   Stronger compression can reduce file size but generally makes conversion
   slower and does not guarantee faster SED reads. Changing this value causes
   pySED to rebuild its converted cache automatically.

   This parameter has no effect with ``trajectory_read_mode = direct`` or
   when ``trajectory_file`` already points to a native NetCDF file.

**Default**
   ``1``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`trajectory_prefetch <trajectory_prefetch>`
