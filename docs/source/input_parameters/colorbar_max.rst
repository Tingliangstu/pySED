colorbar_max
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   colorbar_max = -10

**Meaning**
   Upper limit of the natural-log color scale
   :math:`\ln\left(\Phi/(1\,\mathrm{eV\,THz^{-1}})\right)`. The ratio
   inside the logarithm makes this parameter dimensionless. Increase it when
   the strongest branches saturate the plot; decrease it when weak branches
   are too hard to see.

**Default**
   ``None``.

**Notes**
   The value must be greater than ``colorbar_min``. Negative values are
   valid. Automatic limits are obtained from the extrema of the natural-log
   SED values. See :doc:`../sed_units` for the unit convention.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_min <colorbar_min>`
- :doc:`use_contourf <use_contourf>`
