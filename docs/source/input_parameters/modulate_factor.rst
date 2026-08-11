modulate_factor
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   modulate_factor = 2

**Meaning**
   Removes this many samples from each end of a peak's local valley-bounded
   fitting range before the independent Lorentz/DHO fit. That range is the
   complete basin between the smoothed valleys separating neighboring
   detected peaks. The upper boundary
   itself is excluded by the legacy NumPy slice convention.

**Default**
   ``0``.

**Notes**
   Use sparingly. A nonzero value discards part of the measured peak tail and
   can change the fitted HWHM. It does not introduce a fitted background or
   a multi-peak decomposition.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
