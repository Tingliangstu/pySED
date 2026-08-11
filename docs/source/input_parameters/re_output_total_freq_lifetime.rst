re_output_total_freq_lifetime
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   re_output_total_freq_lifetime = 1

**Meaning**
   Rewrites ``Lifetime/Fitting-All-Qpoints.Fre_lifetime`` and regenerates
   ``Fitting-Frequency-Lifetime.png`` after you re-fit a selected
   q-point. Use it when one q-point fit was poor and you want to adjust
   ``peak_min_significance`` without re-fitting all q-points.

**Default**
   ``0``.

**Notes**
   Use this with ``lorentz_fit_all_qpoint = 0`` after improving a single
   q-point fit. It is active only in single-Q fitting mode.

   When ``lorentz_fit_all_qpoint = 1``, pySED always rebuilds
   ``Lifetime/Fitting-All-Qpoints.Fre_lifetime`` and redraws
   ``Fitting-Frequency-Lifetime.png`` after all Q points have been fitted.
   In that mode, ``re_output_total_freq_lifetime`` is ignored and should
   normally be left at ``0``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`qpoint_slice_index <qpoint_slice_index>`
- :doc:`lorentz_fit_all_qpoint <lorentz_fit_all_qpoint>`
- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
