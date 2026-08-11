Tips
====

Use explicit actions in reproducible runs
-----------------------------------------

Use ``action = compute``, ``plot``, or ``fit`` in scripts.
``thinking`` is convenient interactively, but it only detects whether files
exist; it cannot prove that an old file matches the current input.

Tune fitting locally first
--------------------------

Start with ``lorentz_fit_all_qpoint = 0``, select several representative
Q points, and inspect every fit figure. Only then fit all Q points. Adjust
``peak_min_significance`` before optimizer-specific settings.

Choose cache or direct mode deliberately
----------------------------------------

Use ``trajectory_read_mode = cache`` when the same text trajectory will be
analyzed repeatedly. Use ``direct`` for a one-off run or when an intermediate
NetCDF file is undesirable. Native NetCDF is read directly in either case.

Keep sampling divisible
-----------------------

Check that ``total_num_steps / output_data_stride`` is an integer multiple
of ``num_blocks``. More blocks improve averaging but shorten each FFT and
therefore reduce frequency resolution.
