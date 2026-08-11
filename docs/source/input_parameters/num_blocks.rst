num_blocks
~~~~~~~~~~

**Syntax**

.. code-block:: text

   num_blocks = 5

**Meaning**
   Divides the requested trajectory frames into equal time blocks. pySED
   calculates one spectrum per block and averages the spectra. The requested
   number of saved frames must be at least ``num_blocks`` and divisible by it.

**Default**
   ``5``.

**Tradeoff**
   More blocks reduce the number of frames held and processed in each block
   and provide more spectral averaging, but shorten each block and coarsen the
   FFT frequency spacing. Fewer blocks improve frequency resolution but use
   larger trajectory and working arrays. Do not change this parameter only to
   tune I/O performance without considering the resulting spectrum.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`total_num_steps <total_num_steps>`
- :doc:`output_data_stride <output_data_stride>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
