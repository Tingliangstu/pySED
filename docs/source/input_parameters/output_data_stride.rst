output_data_stride
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   output_data_stride = 50

**Meaning**
   Number of MD steps between saved trajectory frames.

**Default**
   ``0``. This must be set.

**Example**
   If GPUMD uses ``dump_exyz 50 1``, set ``output_data_stride = 50``.

   With ``time_step = 1`` fs and ``output_data_stride = 50``, pySED samples the
   trajectory every 50 fs. The Nyquist frequency is therefore 10 THz.

**Notes**
   This is a sampling parameter, not only a file-reading parameter. pySED
   performs a Fourier transform of the sampled velocities, following the SED
   formulation in which kinetic energy is represented in the frequency domain.
   By the Shannon-Nyquist sampling theorem, the saved trajectory interval must
   be short enough to represent the fastest phonons you want to analyze.

   Define the saved-frame interval as

   .. math::

      \Delta t_{\mathrm{save}}
      =
      \mathrm{time\_step}
      \times
      \mathrm{output\_data\_stride}.

   If ``time_step`` is in fs, the sampling frequency and Nyquist frequency are

   .. math::

      f_s(\mathrm{THz})
      =
      \frac{1000}{\Delta t_{\mathrm{save}}(\mathrm{fs})},
      \qquad
      f_{\mathrm{Nyquist}}(\mathrm{THz})
      =
      \frac{500}{\Delta t_{\mathrm{save}}(\mathrm{fs})}.

   pySED keeps the positive half of the FFT spectrum. Frequencies close to the
   Nyquist limit are the highest frequencies available from the saved
   trajectory; modes above this limit cannot be recovered and may alias into
   lower frequencies.

   The frequency spacing is controlled mainly by the time length of each
   averaging block:

   .. math::

      \Delta f(\mathrm{THz})
      =
      \frac{1000 \times \mathrm{num\_splits}}
      {\mathrm{total\_num\_steps} \times \mathrm{time\_step(fs)}}.

   Therefore, choose ``output_data_stride`` small enough for the maximum phonon
   frequency of interest, then use a sufficiently long production trajectory to
   improve frequency resolution. As a quick rule, for a target maximum
   frequency ``f_max`` in THz,

   .. math::

      \mathrm{output\_data\_stride}
      \le
      \frac{500}
      {\mathrm{time\_step(fs)} \times f_{\max}(\mathrm{THz})}.

   For example, if ``time_step = 1`` fs and modes up to 20 THz are needed,
   choose ``output_data_stride <= 25``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`time_step <time_step>`
- :doc:`total_num_steps <total_num_steps>`
- :doc:`num_splits <num_splits>`
- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`num_atoms <num_atoms>`
