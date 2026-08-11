backend
~~~~~~~

**Syntax**

.. code-block:: text

   backend = numpy

**Meaning**
   Selects the SED compute backend. ``numpy`` uses one process when
   ``max_cores = 1`` and persistent CPU workers when ``max_cores > 1``.
   ``cupy`` uses one CUDA GPU per pySED process. Plotting, peak detection,
   and spectral fitting remain CPU operations.

**Allowed values**
   ``numpy`` or ``cupy``.

**Default**
   ``numpy``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`max_cores <max_cores>`
- :doc:`trajectory_prefetch <trajectory_prefetch>`
