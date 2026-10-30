Units and Normalization
=======================

Public SED quantity
-------------------

pySED writes a one-sided density per ordinary frequency,

.. math::

   \Phi_f(\mathbf q,f),

in eV/THz. The corresponding ``<prefix>.THz`` values are ordinary frequency
in THz, so integrating a spectrum over THz gives energy in eV.

The discrete Fourier transform includes the squared sample spacing
``dt**2``. Each retained positive-frequency bin includes its negative-frequency
partner, while DC is counted once. The current frequency grid retains FFT
indices from zero through ``N // 2 - 1`` for a block of ``N`` frames; for
even-length blocks, the Nyquist bin is excluded. This convention is unchanged
by the architectural refactor.

Display conventions
-------------------

The dispersion heatmap displays

.. math::

   \ln\!\left[
   \frac{\Phi_f(\mathbf q,f)}
        {1\,\mathrm{eV\,THz^{-1}}}
   \right],

so ``colorbar_min`` and ``colorbar_max`` are dimensionless natural-log
values. Single-Q plots and fitted amplitudes remain in physical eV/THz, usually
shown on a logarithmic vertical axis.

``peak_min_significance`` is dimensionless. Fitted center frequencies and
HWHM are in THz, and linewidth-derived lifetimes are in ps.
