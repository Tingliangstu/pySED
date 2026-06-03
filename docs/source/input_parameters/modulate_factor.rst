modulate_factor
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   modulate_factor = 2

**Meaning**
   Narrows the fitting window around each detected peak before the Lorentzian
   fit. This can help when neighboring peaks or background intensity disturb
   the fit.

**Default**
   ``0``.

**Notes**
   This can help when neighboring peaks or noisy baselines distort the fit.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`lorentz <lorentz>`
- :doc:`peak_height <peak_height>`
- :doc:`peak_prominence <peak_prominence>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_cutoff <lorentz_fit_cutoff>`
