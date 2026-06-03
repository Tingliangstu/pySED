use_parallel
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   use_parallel = 1

**Meaning**
   Enables multiprocessing over q-points. It can speed up SED calculation, but
   each worker needs memory, so serial mode can be safer on memory-limited
   machines.

**Default**
   ``1``.

**Notes**
   On Windows, or when memory is limited, set ``use_parallel = 0`` or use a small
   ``max_cores``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`max_cores <max_cores>`
