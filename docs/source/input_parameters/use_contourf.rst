use_contourf
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   use_contourf = 1

**Meaning**
   Controls whether pySED uses contour-style plotting for the SED map. Enable
   it when contour filling gives a cleaner multi-path figure.

**Default**
   ``0``.

**Notes**
   pySED uses contour plotting automatically for multi-segment q-paths. For a
   single path, set this to ``1`` when the contour style gives a cleaner figure.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`plot_SED <plot_SED>`
- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_min <colorbar_min>`
- :doc:`colorbar_max <colorbar_max>`
