lorentz_fit_freq_max
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   lorentz_fit_freq_max = 20

**Meaning**
   Maximum frequency included in spectral peak fitting, in THz. Use it to fit
   only the frequency range where the SED peaks are meaningful and avoid
   fitting noisy high-frequency regions. Combine it with
   ``lorentz_fit_freq_min`` to fit a finite interval.

**Default**
   ``None``. When both ``lorentz_fit_freq_min`` and
   ``lorentz_fit_freq_max`` are ``None`` or omitted, pySED fits the complete
   available frequency interval.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`modulate_factor <modulate_factor>`
