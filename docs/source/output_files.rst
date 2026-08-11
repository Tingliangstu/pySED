Output Files
============

Numerical and plotting outputs
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Output
     - Contents
   * - ``<prefix>.SED``
     - One-sided SED in eV/THz; rows are positive-frequency bins and columns
       follow ``<prefix>.Qpts``.
   * - ``<prefix>.Qpts``
     - Retained reduced Q coordinates.
   * - ``<prefix>.THz``
     - Positive ordinary-frequency grid in THz.
   * - ``<prefix>.Q_distances_and_labels``
     - Cumulative path distances and high-symmetry labels.
   * - ``<prefix>-SED.png``
     - Total SED dispersion heatmap.
   * - ``SED-<q-index>-qpoint.png``
     - Requested unfitted single-Q spectrum.

With ``output_partial = 1``, directional matrices are written beneath
``<prefix>_partial_SED/`` as
``<basename>.SED_<Element>_<x|y|z>``.

Fitting outputs
---------------

.. code-block:: text

   Fitting-Qpoint/
     Fitting-<q-index>-qpoint.png

   Lifetime/
     Fitting-<q-index>-qpoint.Fre_lifetime
     Fitting-All-Qpoints.Fre_lifetime

   Fitting-Frequency-Lifetime.png

Lifetime files contain fitted center frequency in THz and linewidth-derived
lifetime in ps. They do not encode branch assignment; use the numbered figures
and terminal diagnostics to assess each fit.

Trajectory cache
----------------

For text input with ``trajectory_read_mode = cache``, pySED creates or reuses
``<trajectory-name>.pysed.nc`` beside the source. Native NetCDF and direct
text reading do not create this cache.
