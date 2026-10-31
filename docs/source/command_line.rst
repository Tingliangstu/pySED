Command Line
============

Invocation
----------

.. code-block:: bash

   pySED [input-file]
   pySED -h
   pysed [input-file]
   PYSED [input-file]

``pySED`` is the canonical spelling. When no path is supplied, pySED reads
``input_SED.in`` in the current directory. To use another file name, pass it
explicitly, for example ``pysed input.in``. If ``input_SED.in`` is missing,
running without a path reports an error; it does not select another file.

Workflow actions
----------------

The ``action`` parameter controls work performed by one invocation:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Value
     - Behavior
   * - ``thinking``
     - Inspect existing outputs and run the next missing SED stage.
   * - ``compute``
     - Read the trajectory and overwrite the main numerical SED outputs.
   * - ``plot``
     - Read existing numerical SED files and regenerate requested figures.
   * - ``fit``
     - Read existing SED files, detect peaks, fit line shapes, and write
       lifetime data.

Automatic mode does not prove that old outputs were generated with the current
input file. Use an explicit action when reproducibility matters, and archive or
remove stale outputs before recalculating.

Exit status
-----------

Successful work exits with status zero. Invalid parameters, missing input
files, unreadable trajectories, missing plot/fit prerequisites, and numerical
fit-stage failures produce an error and a nonzero status. This makes the CLI
safe to use in batch schedulers and shell pipelines.
