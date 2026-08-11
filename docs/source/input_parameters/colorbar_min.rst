colorbar_min
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   colorbar_min = -22

**Meaning**
   Lower limit of the natural-log color scale
   :math:`\ln\left(\Phi/(1\,\mathrm{eV\,THz^{-1}})\right)`. The ratio
   inside the logarithm makes this parameter dimensionless. Increase it to
   hide weak background noise; decrease it to reveal low-intensity branches.

**Default**
   ``None``.

**Notes**
   Run an initial plot with the default ``None`` first, then tune this value
   for contrast. Negative values are valid. Automatic limits are obtained
   from the extrema of the natural-log SED values. See
   :doc:`../sed_units` for the unit convention.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_max <colorbar_max>`
- :doc:`use_contourf <use_contourf>`
