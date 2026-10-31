trajectory_file
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_file = ../gpumd_run/dump.xyz

**Meaning**
   Path to one supported trajectory:

   - GPUMD extended XYZ with positions, velocities, and lattice information;
   - one LAMMPS custom dump with atom IDs, positions, velocities, and box
     information;
   - compatible GPUMD, LAMMPS, or pySED NetCDF.

   Relative paths are resolved from the directory containing the input file
   (``input_SED.in`` by default).
   pySED reads exactly the named file; it does not choose a neighboring
   ``.xyz`` or ``.nc`` file automatically.

**Default**
   ``dump.xyz``.

**Notes**
   For text input, ``trajectory_read_mode`` selects direct parsing or a
   reusable ``<source-name>.pysed.nc`` cache beside the source trajectory. A named
   native NetCDF file is always read directly.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`trajectory_prefetch <trajectory_prefetch>`
- :doc:`lammps_unit <lammps_unit>`
