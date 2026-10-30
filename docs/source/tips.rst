Tips
====

This page collects small workflow tips for running pySED more efficiently.

Use a Run Script
----------------

For repeated calculations, it is convenient to use a small shell script to
switch ``input.in`` between compute mode and plot/fitting mode. This is
especially useful on clusters, where the same folder may be submitted several
times while tuning plotting and Lorentzian fitting parameters.

The older silicon example contains a simple script:
`run_SED.sh <https://github.com/Tingliangstu/pySED/blob/codex/pysed-3.0-architecture/example/For_old_version_example/Silicon/SED/run_SED.sh>`_.
The idea is:

1. set ``action = compute`` and run pySED once to compute ``.SED``, ``.Qpts``, and
   ``.THz`` files;
2. set ``action = plot`` and run pySED again to plot the existing SED data;
3. use ``action = fit`` while tuning ``peak_min_significance``,
   ``lorentz_fit_freq_max``, and related parameters.

A more robust Bash version is:

.. code-block:: bash

   #!/usr/bin/env bash
   set -euo pipefail

   input_file="input.in"

   set_param() {
       local key="$1"
       local value="$2"
       sed -i -E "s|^(${key}[[:space:]]*=[[:space:]]*).*|\\1${value}|" "${input_file}"
   }

   # First run: compute SED from the trajectory.
   set_param action compute
   pysed "${input_file}"

   # Second run: read existing SED data and plot.
   set_param action plot
   pysed "${input_file}"

This avoids changing a fixed line number. It is safer when comments or new
parameters are added to ``input.in``.

On Windows PowerShell, the same idea can be written as:

.. code-block:: powershell

   $inputFile = "input.in"

   function Set-SedParam($key, $value) {
       $content = Get-Content $inputFile
       $content = $content -replace "^($key\s*=\s*).*", "`${1}$value"
       Set-Content -Path $inputFile -Value $content
   }

   Set-SedParam "action" "compute"
   pysed $inputFile

   Set-SedParam "action" "plot"
   pysed $inputFile

Keep Compute and Plot Modes Separate
------------------------------------

Use ``action = compute`` only when the trajectory, ``basis.in``, q-path, and cell
settings are ready. After SED files are written, use ``action = plot`` for
plotting; use ``action = fit`` for fitting. This avoids rereading trajectories each time
you only want to tune figure or fitting parameters.

Sync Documentation to ReadTheDocs
---------------------------------

The online manual is built from the files under ``docs/source``. The repository
also contains ``.readthedocs.yml``, so ReadTheDocs rebuilds the manual after the
documentation changes are pushed to GitHub.

Recommended update workflow:

.. code-block:: bash

   python -m sphinx -b html docs/source docs/_build/html
   git add README.md docs/source example publications .readthedocs.yml
   git commit -m "Update pySED documentation"
   git push

For the publication list, edit only ``publications/readme.md``. During each
Sphinx build, pySED copies that file into the manual as the ``Publications``
page, so the online page stays synchronized after ReadTheDocs rebuilds.
