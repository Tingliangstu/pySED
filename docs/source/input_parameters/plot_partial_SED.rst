plot_partial_SED
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   plot_partial_SED = 3
   plot_partial_SED = 3 x

**Meaning**
   Plot the partial SED for a one-based atom type index. With only the type
   index, pySED sums x, y, and z directions. With a direction, pySED plots only
   that component.

**Default**
   ``0``.

**Allowed directions**
   ``x``, ``y``, or ``z``.

**Notes**
   Requires partial files generated earlier with ``output_partial = 1``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`output_partial <output_partial>`
- :doc:`plot_SED <plot_SED>`
- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_min <colorbar_min>`
