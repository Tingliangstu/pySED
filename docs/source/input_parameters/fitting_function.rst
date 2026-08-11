fitting_function
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   fitting_function = auto

**Meaning**
   Selects the spectral peak line shape. ``lorentz`` uses a Lorentzian with
   fitted HWHM. ``dho`` uses the velocity-spectrum damped harmonic oscillator
   and reports half of its FWHM as the equivalent HWHM. ``auto`` fits both
   candidates on the same local valley-bounded peak range and compares AICc.
   If candidates are within two AICc units, or if the fitted DHO is in
   the weak-damping limit ``hwhm / center < 0.05``, it selects Lorentz.

**Default**
   ``auto``.

**DHO interpretation**
   Overdamped DHO values are retained in the SED lifetime data for qualitative
   comparison, but a linewidth-derived ``tau_SED`` is not a strict phonon
   lifetime in that regime. A peak is marked incomplete in the terminal when both fitted
   half-maximum points are not contained in its actual fitting range. Its
   reported HWHM is still the width of the
   complete analytic model, estimated by extrapolating the fitted line shape
   beyond the available side of the peak; it is not a complete FWHM directly
   measured from two observed half-maximum crossings. Its linewidth-derived
   lifetime should therefore be compared qualitatively only. In strongly
   anharmonic systems, do not interpret an overdamped result as a strict phonon
   lifetime.

**Related parameters**

- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`peak_min_significance <peak_min_significance>`

:doc:`Back to Peak Detection and Line-Shape Fitting <../peak_fitting>`
