trajectory_prefetch
~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_prefetch = 1

**Meaning**
   This is a Boolean switch and accepts only ``0`` or ``1``. ``0`` disables
   prefetch. When set to ``1``, one background thread prepares the next
   trajectory block while the current block is used for SED computation.
   Exactly one block is prefetched at most; the value is not a configurable
   prefetch depth, and block order is unchanged.

   The same prefetch layer is used by direct GPUMD/LAMMPS text, converted
   NetCDF caches, and native NetCDF trajectories. It can be most helpful when
   text parsing or NetCDF reads would otherwise leave the compute backend
   waiting.

**Default**
   ``1``.

**Memory tradeoff**
   Prefetch retains one additional raw block of positions and velocities in
   CPU memory. Set it to ``0`` when memory is limited. Its speedup depends on
   the filesystem, block size, and CPU/GPU backend, so compare elapsed times on
   the actual trajectory.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`output_data_stride <output_data_stride>`
