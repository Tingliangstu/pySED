peak_max_hwhm
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   peak_max_hwhm = 0.1

**Meaning**
   Maximum allowed fitted half-width at half maximum (HWHM). It prevents very
   broad or failed Lorentzian fits from being accepted as physical peaks.

**Default**
   ``1e6``.

**Notes**
   Set a smaller value when broad or noisy peaks cause unreasonable fits.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_height <peak_height>`
- :doc:`peak_prominence <peak_prominence>`
- :doc:`lorentz <lorentz>`
- :doc:`lorentz_fit_cutoff <lorentz_fit_cutoff>`
- :doc:`modulate_factor <modulate_factor>`
