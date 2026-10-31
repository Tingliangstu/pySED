Code Architecture and Reading Guide
===================================

pySED is an SED-only package. The ``sed/`` package contains the calculation,
and ``postprocess/`` contains result reading, writing, fitting, and plotting.
The package root holds the command entry point, input configuration, and
workflow. Existing module names are retained inside their new directories.

Module responsibilities
-----------------------

.. code-block:: text

   pySED/
     __init__.py         Package version export
     version.py          Version number
     main.py             CLI entry, overview help, and error reporting
     templates/          Commented input template for pysed init
     parser.py           Read the input file
     parameters.py       One parameter table and its validation
     pipeline.py         Select compute, plot, fit, or thinking
     sed/
       construct_BZ.py   Crystal geometry, basis mapping, and Q paths
       Phonon.py         Prepare topology; calculate and average time blocks
       sed_kernel.py     Spatial projection, FFT, power, and normalization
     postprocess/
       FileIO.py         Read/write analysis results
       Lorentz.py        Detect and fit spectral peaks; no file/plot output
       Plot_SED.py       Plot spectra and fitted peaks
       output_files.py   Names of all analysis output files
     backends/
       __init__.py       Select CPU or GPU execution
       cpu.py            Serial batches, worker pool, and shared memory
       gpu.py            CuPy arrays, GPU batches, memory checks, and timing
       common.py         Backend memory estimates and progress messages
     io/                 Trajectory formats, cache, and block iteration
     structure/          Input structures, supercells, and basis.in generation

The calculation direction is ``pipeline`` -> ``sed.Phonon`` -> ``backends``
-> ``sed.sed_kernel``. The pipeline also prepares the trajectory through
``io`` and the geometry through ``sed.construct_BZ``. After computation, it
writes the spectrum through ``postprocess.FileIO``; plot and fit actions
operate on saved spectra through ``postprocess``.
The kernel does not import the pipeline, trajectory readers, or backends.
CuPy is loaded only when the GPU backend is entered. Both new package
initializers are lightweight: importing the numerical kernel does not load
the other calculation modules, plotting, or a GPU backend.

Q-path construction is contained in ``sed.construct_BZ.BZPathHelper``.
``build_commensurate_path`` normalizes the endpoints, obtains exact path
parameters through the GCD/Bézout solver ``_commensurate_parameters``, and
checks the generated coordinates against the supercell. Input-file fractions
remain exact. Python float inputs such as ``7 / 15`` can recover their rational
value when it reproduces the same float exactly. See :doc:`theory` for the
integer-lattice derivation.

Start with the serial calculation
---------------------------------

Use these settings to follow one path without multiprocessing or prefetch:

.. code-block:: text

   action = compute
   device = cpu
   max_cores = 1
   trajectory_prefetch = 0

Read these functions in order:

1. ``pipeline.step_compute_sed`` prepares the trajectory and Q-point geometry.
2. ``sed.Phonon.spectral_energy_density.compute_sed`` prepares the basis metadata
   and opens the selected backend for the lifetime of the calculation.
3. ``sed.Phonon._get_simulation_data`` takes velocities and averages reference
   positions within the current block.
4. ``sed.Phonon._loop_over_blocks`` requests one spectrum from the backend,
   applies block normalization, and adds it to the running sum.
5. ``sed.sed_kernel.calculate_q_batch`` contains the spatial projection and FFT.
6. ``sed.sed_kernel.finalize_spectrum`` averages blocks, folds frequencies, and
   converts the output to eV/THz unless ``legacy_unit`` is enabled.

Names prefixed with ``sed.Phonon.`` in steps 3 and 4 refer to methods of
``spectral_energy_density``. There is no need to read the GPU implementation
to understand the SED formula.

The numerical core
------------------

Let ``T`` be frames per block, ``C`` primitive cells, ``B`` basis atoms,
and ``Q`` Q points in a batch. For total SED the calculation is:

.. code-block:: python

   # velocities: (T, C * B, 3); cell_vectors: (C, 3); qpoints: (Q, 3)
   phases = exp(1j * (cell_vectors @ qpoints.T))  # (C, Q)
   for basis in basis_atoms:
       projected = tensordot(
           velocities[:, basis.indices, :], phases, axes=([1], [0])
       )                                      # (T, 3, Q)
       transformed = fft(projected, axis=0)    # (frequency, 3, Q)
       power = abs(transformed)**2 * basis.mass_in_kg
       spectrum += power.sum(axis=1)          # (frequency, Q)

Partial SED retains the type and Cartesian-direction axes instead of summing
them. The backend result is unnormalized. Block scaling uses
``dt**2 / (4*pi*block_duration*num_unit_cells)`` by default; finalization retains
DC and positive frequencies, excludes an even-length Nyquist bin, and converts
to eV/THz. The optional :doc:`legacy_unit <input_parameters/legacy_unit>` mode
uses the same computation path with 2.x normalization and pair averaging.

Version 2.4 performed the time FFT before summing over cells. The batched
implementation sums cells first, then performs the FFT. The operations
commute because the block's phase reference is time independent and the FFT
is linear. The normalization changes between 2.4 and 3.0 are separate from
this architectural refactor; see :doc:`sed_units` and :doc:`migration_v3`.

Backend contract and resource ownership
---------------------------------------

There are two concrete backend classes, without a plugin registry or an
inheritance hierarchy. Both provide:

* A context manager to open and reliably release their resources.
* ``prepare`` for the first block's memory estimates and GPU batch selection.
* ``compute(velocities, cell_vectors, result)`` to fill one NumPy result array.
* ``report`` for optional timing information after the calculation.

The CPU backend owns its worker pool and shared velocity buffers. The GPU
backend owns its device arrays, CUDA events, and memory pools. Neither
backend changes physical normalization, frequency selection, or output files.
Timing wrappers are supplied to the kernel as ordinary projection/FFT
functions; no CUDA-event logic is embedded in the formula.

Fitting is separate from output
-------------------------------

The pipeline makes these three operations explicit:

.. code-block:: python

   from pySED.postprocess import FileIO, Lorentz, Plot_SED

   fit_result = Lorentz.lorentz(data, params)
   FileIO.write_phonon_lifetime(fit_result, params)
   Plot_SED.plot_slice(data, params, fit_result=fit_result)

The fitter returns its parameters and per-peak curves in
``fit_result.peak_fits``. It does not write files, draw figures, or attach
results to ``params``. Calling ``plot_slice(data, params)`` without a fit
result always gives a plain spectrum. All-Q fitting reads the SED data once
and uses a separate parameter copy for each Q point.

This is a deliberate change for Python callers of the earlier development
version: constructing ``Lorentz.lorentz`` no longer produces output files
as a side effect. The CLI still performs all three operations.

Where to make a change
----------------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Intended change
     - Starting point
   * - SED formula or normalization
     - ``sed/sed_kernel.py`` and numerical regression tests
   * - Time-block averaging or basis metadata
     - ``sed/Phonon.py``
   * - CPU performance or shared-memory handling
     - ``backends/cpu.py``
   * - GPU performance, timing, or memory
     - ``backends/gpu.py``
   * - Trajectory format or unit conversion on input
     - ``io/text.py`` or ``io/netcdf.py``
   * - Q-point geometry
     - ``sed/construct_BZ.py``
   * - Peak detection or line-shape fitting
     - ``postprocess/Lorentz.py``
   * - Output filenames or directories
     - ``postprocess/output_files.py``
   * - Spectrum files or plotting
     - ``postprocess/FileIO.py`` or ``postprocess/Plot_SED.py``
   * - Parameters or CLI workflow
     - ``parameters.py``, ``parser.py``, or ``pipeline.py``

Rules for future additions
--------------------------

* Keep pySED SED-only. Add input options to the single ``PARAMETERS`` table;
  do not introduce method registries or a plugin framework.
* Keep mathematics separate from file output and plotting. Pass results
  explicitly instead of storing fit results or plotting flags on ``params``.
* Let each backend own its execution resources. A numerical feature should
  have one shared implementation unless the hardware genuinely requires
  different algorithms.
* Put output naming in ``postprocess/output_files.py``. Writers, readers, and workflow
  readiness checks must use the same names.
* Create a module only for a distinct responsibility. Do not split every
  small helper into its own file or add a base class for two simple backends.
* Add a serial numerical reference test before optimizing a new calculation.

The old-parameter error messages and the fallback input filename remain
useful migration guidance; they do not execute a second legacy workflow.
Computing new results does not automatically delete old-format result files.

The ``thinking`` action checks whether expected files exist, including all
four SED data files. It does not compare calculation settings or file contents.
After changing the trajectory or calculation parameters, use ``action = compute``
explicitly.

Python imports and maintenance
------------------------------

Programmatic imports follow the calculation and postprocessing packages:

.. code-block:: python

   from pySED.sed.Phonon import spectral_energy_density
   from pySED.sed.construct_BZ import BZ_methods
   from pySED.postprocess import FileIO, Plot_SED

Scripts using the previous top-level analysis imports must use the paths
above. Trajectory readers remain under ``pySED.io``, CPU/GPU execution under
``pySED.backends``, and structure generation under ``pySED.structure``.
The directory reorganization leaves input parameters, CLI commands, SED units,
fitting models, and output filenames unchanged.

For synchronization from MDtrace, port a mathematical change into the kernel
and a resource-management change into the corresponding backend. Do not copy
the upstream ``Phonon.py`` wholesale over the split implementation.

Run the regression suite from the repository root:

.. code-block:: bash

   python -m pytest -q

Tests cover reference-kernel agreement, partial spectra, normalization,
serial/multiprocess agreement, GPU batching through an adapter, cleanup after
errors, and commensurate q-path construction. A real CuPy numerical test is
conditional on CuPy and a usable CUDA device.
