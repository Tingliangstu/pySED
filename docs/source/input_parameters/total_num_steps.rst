total_num_steps
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   total_num_steps = 500000

**Meaning**
   Total number of MD steps in the production trajectory used for SED.

**Default**
   ``0``. This must be set.

**Notes**
   The number of trajectory frames used by pySED is
   ``total_num_steps / output_data_stride``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`num_atoms <num_atoms>`
- :doc:`time_step <time_step>`
- :doc:`output_data_stride <output_data_stride>`
