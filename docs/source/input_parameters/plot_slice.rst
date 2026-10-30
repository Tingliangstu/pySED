plot_slice
~~~~~~~~~~

**Syntax**

.. code-block:: text

   plot_slice = 1

**Meaning**
   Controls whether pySED plots a one-dimensional SED spectrum at a selected
   q-point. This is useful for checking ``peak_min_significance`` before
   Lorentz/DHO fitting. The plotted quantity is
   :math:`\Phi(\mathbf q,\omega)` in
   :math:`\mathrm{eV\,THz^{-1}}`; only the y-axis scaling is logarithmic.

**Default**
   ``0``.

**Notes**
   Use this before all-q-point spectral fitting to tune peak detection settings.
   Unlike the dispersion heatmap, this plot does not replace the physical SED
   values with their natural logarithms. See :doc:`../sed_units`.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`action <action>`
- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_min <colorbar_min>`
- :doc:`colorbar_max <colorbar_max>`
