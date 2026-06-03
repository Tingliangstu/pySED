num_splits
~~~~~~~~~~

**Syntax**

.. code-block:: text

   num_splits = 5

**Meaning**
   Number of blocks used to split the trajectory before averaging SED. More
   blocks can make the SED map smoother, but each block is shorter, so the
   frequency resolution becomes coarser.

**Default**
   ``1``.

**Notes**
   Values such as ``5`` or ``10`` often give smoother SED maps. Each split uses
   fewer frames, so very large values reduce frequency resolution.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`compress <compress>`
