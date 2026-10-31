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
   Provide ``num_qpaths + 1`` triples. Explicit fractions such as ``1/3`` or
   ``7/15`` are preserved exactly and are preferred for rational endpoints.
   Decimal inputs within ``1e-5`` of a fraction with denominator at most 12
   are normalized to that fraction; for example, ``0.33333`` becomes ``1/3``.
   Other decimal input is retained exactly as written. To preserve a value
   near a simple fraction without normalization, write it explicitly as a
   fraction, for example ``33333/100000``.

   pySED uses the same normalized endpoints for path selection and coordinate
   output. Only points satisfying the supercell periodicity condition are
   included; a supplied endpoint is not automatically an allowed q-point.
   A repeated endpoint contributes one point only if it satisfies that
   condition. See :doc:`../theory` for the GCD/Bézout construction.

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
