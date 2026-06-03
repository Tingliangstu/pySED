Tips
====

This page collects small workflow tips for running pySED more efficiently.

Use a Run Script
----------------

For repeated calculations, it is convenient to use a small shell script to
switch ``input_SED.in`` between compute mode and plot/fitting mode. This is
especially useful on clusters, where the same folder may be submitted several
times while tuning plotting and Lorentzian fitting parameters.

The older silicon example contains a simple script:
`run_SED.sh <https://github.com/Tingliangstu/pySED/blob/main/example/For_old_version_example/Silicon/SED/run_SED.sh>`_.
The idea is:

1. set ``plot_SED = 0`` and run pySED once to compute ``.SED``, ``.Qpts``, and
   ``.THz`` files;
2. set ``plot_SED = 1`` and run pySED again to plot the existing SED data;
3. repeat plot/fitting mode while tuning ``peak_height``, ``peak_prominence``,
   ``lorentz_fit_cutoff``, and related parameters.

A more robust Bash version is:

.. code-block:: bash

   #!/usr/bin/env bash
   set -euo pipefail

   input_file="input_SED.in"

   set_param() {
       local key="$1"
       local value="$2"
       sed -i -E "s|^(${key}[[:space:]]*=[[:space:]]*).*|\\1${value}|" "${input_file}"
   }

   # First run: compute SED from the trajectory.
   set_param plot_SED 0
   pysed "${input_file}"

   # Second run: read existing SED data and plot or fit.
   set_param plot_SED 1
   pysed "${input_file}"

This avoids changing a fixed line number. It is safer when comments or new
parameters are added to ``input_SED.in``.

On Windows PowerShell, the same idea can be written as:

.. code-block:: powershell

   $inputFile = "input_SED.in"

   function Set-SedParam($key, $value) {
       $content = Get-Content $inputFile
       $content = $content -replace "^($key\s*=\s*).*", "`${1}$value"
       Set-Content -Path $inputFile -Value $content
   }

   Set-SedParam "plot_SED" "0"
   pysed $inputFile

   Set-SedParam "plot_SED" "1"
   pysed $inputFile

Keep Compute and Plot Modes Separate
------------------------------------

Use ``plot_SED = 0`` only when the trajectory, ``basis.in``, q-path, and cell
settings are ready. After SED files are written, use ``plot_SED = 1`` for
plotting and fitting. This avoids recompressing large trajectories each time
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
