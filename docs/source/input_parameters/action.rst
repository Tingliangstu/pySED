action
~~~~~~

**Syntax**

.. code-block:: text

   action = thinking

**Meaning**
   Selects the SED workflow stage:

   - ``thinking`` runs the next missing stage based on existing outputs;
   - ``compute`` recalculates the numerical SED and overwrites its main data;
   - ``plot`` redraws figures from existing SED data;
   - ``fit`` detects and fits peaks in existing SED data.

   ``thinking`` does not compare the parameters used to create an existing
   output with the current input file. Use an explicit action to force a
   recalculation, redraw, or refit after changing settings.

**Default**
   ``thinking``.

:doc:`Back to Parameter Index <../input_parameters>`
