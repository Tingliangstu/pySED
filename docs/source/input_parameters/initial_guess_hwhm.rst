initial_guess_hwhm
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   initial_guess_hwhm = 0.001

**Meaning**
   Initial half-width at half maximum (HWHM) supplied to each Lorentz or DHO
   nonlinear fit. If fits fail or converge to unreasonable widths, this value
   can help guide the optimizer.

**Default**
   ``0.001``.

**Notes**
   This value has units of THz. The optimizer clips it inside the valid width
   bounds before fitting.

Frequency resolution
--------------------

   ``initial_guess_hwhm`` does not set the frequency resolution. Because
   pySED divides the requested trajectory into ``num_blocks`` equal blocks
   before averaging their spectra, the approximate FFT frequency spacing is

   .. math::

      \Delta f\,[\mathrm{THz}]
      =
      \frac{1}{T_\mathrm{block}\,[\mathrm{ps}]}
      =
      \frac{1000\,\mathtt{num\_blocks}}
      {\mathtt{total\_num\_steps}\,
       \mathtt{time\_step}\,[\mathrm{fs}]}.

   Here

   .. math::

      T_\mathrm{block}\,[\mathrm{ps}]
      =
      \frac{\mathtt{total\_num\_steps}\,
             \mathtt{time\_step}\,[\mathrm{fs}]}
           {1000\,\mathtt{num\_blocks}}.

   ``output_data_stride`` does not appear in the final expression because it
   decreases the number of saved frames and increases the time between frames
   by the same factor. It still controls the Nyquist frequency and must match
   the trajectory output interval.

   For example,

   .. code-block:: text

      total_num_steps    = 300000
      time_step          = 1
      num_blocks         = 5

   gives :math:`T_\mathrm{block}=60\,\mathrm{ps}` and
   :math:`\Delta f\approx 0.0167\,\mathrm{THz}`. Therefore,
   ``initial_guess_hwhm = 0.01`` is of the same order as one frequency bin and
   is a reasonable optimizer starting value for this example. It cannot make a
   linewidth narrower than the available resolution physically resolvable.

   In the SED theory used by pySED, the fitted peak is written as a Lorentzian

   .. math::

      \Phi(\mathbf{q},\omega)
      =
      \frac{I}
      {1+\left[(\omega-\omega_c)/\gamma\right]^2},

   where :math:`\gamma` is the HWHM. The PYSED paper defines the lifetime from
   this linewidth as :math:`\tau = 1/(2\gamma)`. In pySED output, frequencies
   are handled in THz and lifetimes are written in ps using the code convention
   described in the theory page.

   A reasonable initial HWHM helps the nonlinear fit converge to the physical
   linewidth. If ``initial_guess_hwhm`` is much too small, the fit may lock onto
   a very narrow spike or fail for broadened peaks. If it is much too large, the
   fit can over-broaden nearby peaks or converge slowly. For sharp crystalline
   peaks, values such as ``0.0005`` to ``0.005`` THz are often a useful starting
   range, but the best value depends on the material, temperature, trajectory
   length, and frequency resolution.

   Tune this parameter only after ``peak_min_significance`` detects the
   intended peaks. Then check the fitted curve visually before using
   ``lorentz_fit_all_qpoint = 1``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
- :doc:`modulate_factor <modulate_factor>`
