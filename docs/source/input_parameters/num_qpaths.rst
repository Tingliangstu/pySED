num_qpaths
~~~~~~~~~~

**Syntax**

.. code-block:: text

   num_qpaths = 3

**Meaning**
   Number of q-path segments to plot. For example, a path ``G-M-K-G`` has three
   segments, so ``num_qpaths = 3`` and ``q_path`` must contain four coordinate
   triples.

**Default**
   ``None``. Set it explicitly.

**Notes**
   A path with labels ``GMKG`` has three segments: Gamma to M, M to K, and K to
   Gamma.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`q_path_name <q_path_name>`
- :doc:`q_path <q_path>`
