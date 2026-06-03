if_show_figures
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   if_show_figures = 0

**Meaning**
   Controls whether Matplotlib figures are shown on screen after pySED creates
   them. Use ``0`` for batch jobs or remote runs, and ``1`` when working
   interactively.

**Default**
   ``0``.

**Notes**
   Use ``0`` on clusters or headless environments.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`plot_SED <plot_SED>`
- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_min <colorbar_min>`
- :doc:`colorbar_max <colorbar_max>`
