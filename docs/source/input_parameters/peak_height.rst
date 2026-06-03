peak_height
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   peak_height = 8.0e-7

**Meaning**
   Minimum absolute SED intensity passed to ``scipy.signal.find_peaks`` through
   the ``height`` argument.

**Default**
   No parser default.

**Notes**
   This is the first gate for automatic peak detection. A peak must be higher
   than this value before pySED attempts a Lorentzian fit.

   Set this parameter from a single-q-point slice first. Use
   ``plot_slice = 1`` and inspect both the plotted spectrum and the frequencies
   printed by pySED. If ``peak_height`` is too high, weak but real branches are
   missed. If it is too low, random noise and small shoulders may be fitted as
   artificial phonon modes.

   For noisy spectra, tune ``peak_height`` together with
   ``peak_prominence``. Height controls the absolute intensity threshold, while
   prominence controls how clearly a peak stands out from its local background.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`peak_prominence <peak_prominence>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz <lorentz>`
- :doc:`lorentz_fit_cutoff <lorentz_fit_cutoff>`
- :doc:`modulate_factor <modulate_factor>`
