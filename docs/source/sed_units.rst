Units and Normalization
=======================

Public SED quantity
-------------------

pySED writes a one-sided density per ordinary frequency,

.. math::

   \Phi_f(\mathbf q,f),

in eV/THz by default (``legacy_unit = False``). The corresponding
``<prefix>.THz`` values are ordinary frequency in THz, so integrating a spectrum
over THz gives energy in eV.

The discrete Fourier transform includes the squared sample spacing
``dt**2``. Each retained positive-frequency bin includes its negative-frequency
partner, while DC is counted once. For a block of ``N`` frames, the frequency
grid retains DC and every positive-frequency bin. Even-length blocks exclude
the Nyquist bin; odd-length blocks retain the bin at index ``N // 2``.

Historical output scale
-----------------------

Set :doc:`legacy_unit <input_parameters/legacy_unit>` to ``True`` to compare
with pySED 2.x intensities. This mode omits ``dt**2``, averages each
positive/negative frequency pair, and skips the eV/THz conversion. Its output
is labelled ``legacy 2.x scale``; the historical ``J*s`` label does not denote
the same normalized density as the default output. Frequency is still in THz.

This compatibility option will be removed in a future major release.
Use ``legacy_unit = False`` for new calculations. Changing the intensity
scale affects peak amplitudes and spectral integrals; a uniform scaling of
a nonzero-frequency peak does not change its center or linewidth.

SED file headers record the output scale and control subsequent plot labels.
For files without scale metadata, the input ``legacy_unit`` value specifies
the scale; loading does not convert the stored values.

Display conventions
-------------------

With the default physical scale, the dispersion heatmap displays

.. math::

   \ln\!\left[
   \frac{\Phi_f(\mathbf q,f)}
        {1\,\mathrm{eV\,THz^{-1}}}
   \right],

so ``colorbar_min`` and ``colorbar_max`` are dimensionless natural-log
values. Single-Q plots and fitted amplitudes use eV/THz by default, usually
shown on a logarithmic vertical axis. Legacy results use the same plot style
with the explicit ``legacy 2.x scale`` label and their stored intensities.

``peak_min_significance`` is dimensionless. Fitted center frequencies and
HWHM are in THz, and linewidth-derived lifetimes are in ps.
