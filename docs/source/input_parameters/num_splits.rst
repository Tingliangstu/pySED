num_splits
~~~~~~~~~~

**Syntax**

.. code-block:: text

   num_splits = 5

**Meaning**
   Number of trajectory blocks used for block averaging.

**Default**
   ``1``.

**Notes**
   Values such as ``5`` or ``10`` often give smoother SED maps. Each split uses
   fewer frames, so very large values reduce frequency resolution.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`compress <compress>`
