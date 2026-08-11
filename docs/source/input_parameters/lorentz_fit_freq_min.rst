lorentz_fit_freq_min
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   lorentz_fit_freq_min = 2.0

**Meaning**
   Minimum frequency included in spectral peak detection and fitting, in
   THz. Combine it with ``lorentz_fit_freq_max`` to select a finite fitting
   interval:

   .. code-block:: text

      lorentz_fit_freq_min = 2.0
      lorentz_fit_freq_max = 5.0

   This detects and fits peaks from 2 to 5 THz, including both endpoints. The
   fitted single-Q figure displays the same frequency range.

**Default**
   ``None``. Fitting starts from the lowest available frequency, normally
   0 THz.

**Notes**
   The value must be non-negative and smaller than
   ``lorentz_fit_freq_max`` when both are specified. The requested interval
   must contain at least three sampled frequency points. When both frequency
   bounds are ``None`` or omitted, pySED fits the complete available
   frequency interval.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
- :doc:`peak_min_significance <peak_min_significance>`
