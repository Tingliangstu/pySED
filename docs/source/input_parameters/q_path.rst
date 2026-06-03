q_path
~~~~~~

**Syntax**

.. code-block:: text

   q_path = 0 0 0  1/2 0 0  1/3 1/3 0  0 0 0

**Meaning**
   High-symmetry q-points in reduced reciprocal coordinates. These are the
   endpoint coordinates of the path that pySED samples for the SED plot.

**Default**
   No useful default. Set it explicitly.

**Notes**
   Provide ``num_qpaths + 1`` triples. Fractions such as ``1/3`` are accepted.
   If a high-symmetry point is a rational value, prefer fraction input such as
   ``1/3`` or ``1/2`` instead of rounded decimals such as ``0.3333333``. This
   avoids unnecessary rounding when pySED checks which q-points are
   commensurate with the supercell.

   For conventional high-symmetry paths and labels, users may refer to the
   Setyawan-Curtarolo band-structure path paper [Setyawan2010]_.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`num_qpaths <num_qpaths>`
- :doc:`q_path_name <q_path_name>`

References
----------

.. [Setyawan2010] W. Setyawan and S. Curtarolo, "High-throughput electronic
   band structure calculations: Challenges and tools," *Computational
   Materials Science* **49**, 299-312 (2010).
   https://doi.org/10.1016/j.commatsci.2010.05.010
