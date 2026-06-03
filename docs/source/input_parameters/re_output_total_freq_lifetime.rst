re_output_total_freq_lifetime
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   re_output_total_freq_lifetime = 1

**Meaning**
   Rewrites ``TOTAL-LORENTZ-Qpoints.Fre_lifetime`` after you re-fit a selected
   q-point. Use it when one q-point fit was poor and you want to adjust
   ``peak_height`` or ``peak_prominence`` without re-fitting all q-points.

**Default**
   ``0``.

**Notes**
   Use this with ``lorentz_fit_all_qpoint = 0`` after improving a single
   q-point fit.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`qpoint_slice_index <qpoint_slice_index>`
- :doc:`lorentz_fit_all_qpoint <lorentz_fit_all_qpoint>`
- :doc:`lorentz <lorentz>`
- :doc:`peak_height <peak_height>`
- :doc:`peak_prominence <peak_prominence>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
