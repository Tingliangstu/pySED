max_cores
~~~~~~~~~

**Syntax**

.. code-block:: text

   max_cores = 4

**Meaning**
   Maximum number of worker processes used in parallel mode.

**Default**
   ``4``.

**Notes**
   More cores can reduce compute time but increase memory use because each worker
   handles trajectory data.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`use_parallel <use_parallel>`
