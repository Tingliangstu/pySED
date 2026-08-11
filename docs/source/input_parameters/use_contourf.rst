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
   pySED uses ``contourf`` automatically for multiple q-paths. A single
   q-path uses ``imshow`` by default; set this parameter to ``1`` to use
   ``contourf`` instead. A single-q spectrum requested by ``plot_slice`` is a
   separate line plot with physical SED values in
   :math:`\mathrm{eV\,THz^{-1}}` on a logarithmic y-axis.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_min <colorbar_min>`
- :doc:`colorbar_max <colorbar_max>`
