max_cores
~~~~~~~~~

**Syntax**

.. code-block:: text

   max_cores = 4

**Meaning**
   Maximum number of NumPy CPU worker processes. ``max_cores = 1`` runs the SED
   kernel serially; larger values enable multiprocessing.

**Default**
   ``4``.

**Notes**
   More cores can reduce compute time but increase memory use because each worker
   handles trajectory data. This setting is ignored by the single-GPU CuPy SED
   kernel.

:doc:`Back to Parameter Index <../input_parameters>`
