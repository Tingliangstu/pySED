output_data_stride
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   output_data_stride = 50

**Meaning**
   Number of MD simulation steps between two saved trajectory frames. This is
   the same output interval used in the MD run, so it controls how often pySED
   receives atomic positions and velocities.

**Default**
   ``0``. This must be set.

**Example**
   If GPUMD uses ``dump_exyz 50 1``, set ``output_data_stride = 50``.

   With ``time_step = 1`` fs and ``output_data_stride = 50``, pySED samples the
   trajectory every 50 fs. The Nyquist frequency is therefore 10 THz.

**Notes**
   ``output_data_stride`` is a sampling parameter. During the production run,
   atomic velocities and positions are saved every ``output_data_stride`` MD
   steps. According to the Shannon sampling theorem, the maximum frequency that
   can be resolved by SED is

   .. math::

      f_{\max}(\mathrm{THz})
      =
      \frac{500}
      {\mathrm{time\_step(fs)} \times \mathrm{output\_data\_stride}}.

   This is the usual Nyquist limit,
   ``f_max = 1 / (2 * time_step * output_data_stride)``. The factor ``500``
   comes only from unit conversion: when ``time_step`` is given in fs, ``1/fs``
   equals ``1000`` THz, and the Nyquist factor ``1/2`` gives ``1000 / 2 = 500``.

   For example, if ``time_step = 1`` fs and ``output_data_stride = 50``, the
   maximum resolvable frequency is ``10`` THz. Modes above this frequency
   cannot be recovered from the saved trajectory.

   A longer production trajectory gives more saved frames and improves the
   frequency resolution of the SED map. The trade-off is memory and file size:
   saving frames more often or running longer gives better spectral information,
   but requires more storage and memory during analysis.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`time_step <time_step>`
- :doc:`total_num_steps <total_num_steps>`
- :doc:`num_splits <num_splits>`
- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`num_atoms <num_atoms>`
