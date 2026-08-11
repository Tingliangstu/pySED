qpoint_slice_index
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   qpoint_slice_index = 2

**Meaning**
   Zero-based index of the q-point used for single-q-point plotting and
   fitting. Use the ``.Qpts`` output file to identify which index corresponds
   to the q-point you want to inspect.

   For fitting, this selection is used when
   ``lorentz_fit_all_qpoint = 0`` or when that parameter is omitted. With
   ``action = thinking``, pySED fits this Q point only if its corresponding
   ``Lifetime/Fitting-<q-index>-qpoint.Fre_lifetime`` file is missing.

**Default**
   ``0``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_min <colorbar_min>`
- :doc:`colorbar_max <colorbar_max>`
