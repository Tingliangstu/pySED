Command Line
============

Invocation
----------

.. code-block:: bash

   pySED [input-file]
   pySED -h
   pysed [input-file]
   PYSED [input-file]

Create an input file
--------------------

Run ``pysed init`` in your calculation directory to create a basic
``input_SED.in`` template:

.. code-block:: bash

   pysed init

Edit the generated file for your trajectory, structure, and MD sampling, then
run ``pysed`` to start the calculation. If ``input_SED.in`` already exists,
the command preserves it and returns an error.

``pysed init -h`` explains the command without creating a file.

Parameter help
--------------

``pysed -h`` prints the program introduction, references, online manual link,
and a short guide to looking up parameters. Detailed parameter descriptions
are available with ``pysed <parameter> -h``:

.. code-block:: bash

   pysed action -h
   pysed total_num_steps -h
   pysed device --help
   pysed trajectory_file -h

Replace the name with any supported input parameter. Each help page explains
the syntax, meaning, default or required status, examples, and practical notes.
The help is bundled with pySED and works offline, without ``input_SED.in`` or
trajectory data. It does not run a calculation. Unknown parameter names return
an error and, when possible, suggest the closest valid name.

The `online manual <https://pysed.readthedocs.io/>`_ provides the full parameter
guide, workflow examples, and theory notes.

Input file selection
--------------------

``pySED`` is the canonical spelling. When no path is supplied, pySED reads
``input_SED.in`` in the current directory. To use another file name, pass it
explicitly, for example ``pysed input.in``. If ``input_SED.in`` is missing,
running without a path reports an error; it does not select another file.
Except for help requests and the ``init`` command, a single argument is an
input-file path, even if it is named like a parameter. For example,
``pysed action`` reads a file named ``action``; it does not select the workflow
stage. To read an input file literally named ``init``, use ``pysed ./init``.

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
