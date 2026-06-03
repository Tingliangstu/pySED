initial_guess_hwhm
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   initial_guess_hwhm = 0.001

**Meaning**
   Initial guess for the Lorentzian HWHM.

**Default**
   ``0.001``.

**Notes**
   This value is the starting guess for the Lorentzian half-width at
   half-maximum, passed to ``scipy.optimize.curve_fit``. It has units of THz in
   the current pySED fitting workflow.

   In the SED theory used by pySED, the fitted peak is written as a Lorentzian

   .. math::

      \Phi(\mathbf{q},\omega)
      =
      \frac{I}
      {1+\left[(\omega-\omega_c)/\gamma\right]^2},

   where :math:`\gamma` is the HWHM. The PYSED paper defines the lifetime from
   this linewidth as :math:`\tau = 1/(2\gamma)`. In pySED output, frequencies
   are handled in THz and lifetimes are written in ps using the code convention
   described in the theory page.

   A reasonable initial HWHM helps the nonlinear fit converge to the physical
   linewidth. If ``initial_guess_hwhm`` is much too small, the fit may lock onto
   a very narrow spike or fail for broadened peaks. If it is much too large, the
   fit can over-broaden nearby peaks or converge slowly. For sharp crystalline
   peaks, values such as ``0.0005`` to ``0.005`` THz are often a useful starting
   range, but the best value depends on the material, temperature, trajectory
   length, and frequency resolution.

   Tune this parameter only after ``peak_height`` and ``peak_prominence`` detect
   the correct peaks. Then check the fitted curve visually before using
   ``lorentz_fit_all_qpoint = 1``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`peak_height <peak_height>`
- :doc:`peak_prominence <peak_prominence>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_cutoff <lorentz_fit_cutoff>`
- :doc:`lorentz <lorentz>`
- :doc:`modulate_factor <modulate_factor>`
