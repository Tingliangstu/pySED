Scattering Extensions
=====================

This page records the theoretical and implementation basis for the new
experiment-facing scattering modules in pySED.  The design goal is not to copy
dynasor, pynamic-structure-factor, or PySlice, but to provide a pySED-native
route from MD trajectories to experiment-comparable phonon scattering maps with
mode-level interpretation.

The modules currently added are:

- :mod:`pySED.scattering_kernel`: shared phase, density-amplitude, time-FFT,
  and block-statistics kernels.
- :mod:`pySED.probe_weights`: neutron/X-ray probe weights, including
  q-dependent Cromer-Mann X-ray form factors.
- :mod:`pySED.quantum_correction`: optional classical-to-quantum
  detailed-balance corrections and energy-axis conversion.
- :mod:`pySED.q_advisor`: commensurate q-point validation, nearest-q mapping,
  finite-size warnings, and supercell suggestions.
- :mod:`pySED.dsf`: coherent and incoherent neutron dynamic structure factors,
  X-ray weighted DSF, and longitudinal/transverse current spectra.
- :mod:`pySED.mode_visibility`: coherent one-phonon neutron/X-ray
  mode-visibility diagnostics from eigenvectors.
- :mod:`pySED.visibility_diagnostics`: branch-level reason labels for
  polarization selection, basis interference, finite-size q mapping, and weak
  experimental visibility.
- :mod:`pySED.eigen_sed`: phonopy-eigenvector projection for branch-resolved
  SED.
- :mod:`pySED.electron_kinematics`: relativistic electron wavelength and
  small-angle q-EELS momentum-axis calibration.
- :mod:`pySED.eels`: first kinematic q-EELS visibility and extended-zone maps.
- :mod:`pySED.compare_experiment`: broadening, smoothing, normalization,
  residual maps, line cuts, and peak tracking for experimental comparison.
- :mod:`pySED.scattering_plot`: publication-style map, residual, and line-cut
  plotting helpers.
- :mod:`pySED.scattering_workflow`: high-level workflow wrappers that connect
  commensurate q validation, eigenvector SED, extended-zone EELS visibility,
  and optional experiment-map comparison.
- :mod:`pySED.scattering_export`: reproducibility bundles for q-plans,
  simulated maps, visibility diagnostics, and comparison tables.

The core innovation is a unified scattering kernel:

.. math::

   \mathrm{MD\ trajectory}
   \rightarrow
   \left\{
   \Phi(\mathbf{q},\omega),
   S(\mathbf{Q},\omega),
   I_{\mathrm{EELS}}(\mathbf{Q},\omega)
   \right\},

where the same finite-supercell q-point rules are used for SED, DSF, EELS, and
mode projection.

The low-level implementation is centralized in
:mod:`pySED.scattering_kernel`.  This module contains the shared CPU/CuPy-ready
kernels for

.. math::

   e^{i\mathbf{Q}\cdot\mathbf{r}_i(t)},\quad
   \rho_w(\mathbf{Q},t),\quad
   |\mathrm{FFT}_t[\rho_w]|^2,

as well as block slicing and standard-error estimates.  Keeping these
operations in one module is intentional: optional GPU implementations can
replace this kernel layer without changing the DSF/EELS theory or high-level
user API.

The default backend is CPU/NumPy:

.. code-block:: python

   result = compute_dsf(positions, qpoints, dt, backend="cpu")

Users with a compatible CUDA environment can request the optional CuPy backend:

.. code-block:: python

   result = compute_dsf(positions, qpoints, dt, backend="cupy")

The eigenvector-projected SED path uses the same backend mechanism:

.. code-block:: python

   eigen_sed = compute_eigen_sed(..., backend="cupy")
   eels = compute_eels_workflow(..., sed_backend="cupy")
   ixs = compute_one_phonon_workflow(..., sed_backend="cupy")

CuPy is not a pySED hard dependency because CUDA wheels are platform- and
runtime-specific.  Users should install the CuPy package that matches their own
CUDA runtime, for example ``pip install ".[gpu-cuda12x]"`` or
``pip install ".[gpu-cuda11x]"`` from the pySED source tree.  If CuPy is not
installed, pySED raises a clear error and the CPU path remains available.  The
current GPU coverage targets the numerically heavy selected-Q density FFTs,
current correlations, and eigenvector SED projection/FFT.  The lightweight
visibility weighting, q-advisor, and experimental comparison steps remain CPU
NumPy operations.  This keeps GPU acceleration opt-in and avoids forcing
non-GPU users to compile or install GPU packages.

Commensurate Q Points
---------------------

For a primitive cell :math:`\mathbf{p}` and a simulation supercell
:math:`\mathbf{S}`, both written with lattice vectors as rows,

.. math::

   \mathbf{S} = \mathbf{P}\mathbf{p},

where :math:`\mathbf{P}` is an integer repetition matrix.  The reciprocal basis
of the primitive cell is

.. math::

   \mathbf{B}_p = 2\pi(\mathbf{p}^{-1})^T.

A reduced wave vector :math:`\mathbf{q}_{\mathrm{red}}` is compatible with the
periodic MD trajectory only when

.. math::

   \mathbf{q}_{\mathrm{red}}\mathbf{P}^{T} \in \mathbb{Z}^3.

For a diagonal supercell
:math:`\mathbf{P}=\mathrm{diag}(N_x,N_y,N_z)`, this reduces to the familiar
grid

.. math::

   \mathbf{q}_{\mathrm{red}}
   =
   \left(
   \frac{n_x}{N_x},
   \frac{n_y}{N_y},
   \frac{n_z}{N_z}
   \right).

Non-commensurate wave vectors are not simply "slightly interpolated" phonons.
They violate the periodic boundary condition of the finite trajectory.  In a
direct Fourier transform, this appears as finite-window leakage or spurious
intensity.  Therefore, pySED should reject or explicitly label
non-commensurate q points by default.

The :mod:`pySED.q_advisor` module exposes this condition directly.  It can:

- validate requested q points,
- enumerate the finite-supercell commensurate grid,
- construct only the allowed points on a high-symmetry path,
- map an experimental :math:`\mathbf{Q}` to the nearest commensurate point,
  either from reduced coordinates or directly from Cartesian reciprocal-space
  coordinates,
- report the reduced and Cartesian error, and
- suggest diagonal supercell repeats for a target set of experimental q
  points.

This differs from pynamic-structure-factor in an important way.  pynamic uses
direct user-supplied Q points and parallelizes over them, which is efficient
for selected experimental paths and maps.  Its documentation also warns that
the user must choose Q points commensurate with the simulation cell.  Therefore
pynamic's advantage is the direct selected-Q estimator, not a more rigorous
q-point algorithm.  pySED keeps the stricter q-advisor layer and uses the
direct selected-Q estimator where appropriate.

Pynamic-style Selected-Q Grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The selected-Q idea is useful and should be kept.  The point is not that the
Q-point generation itself is mathematically superior; the point is that the
calculation is restricted to the experimental vectors that will actually be
plotted or compared.  For a selected set
:math:`\mathcal{Q}_{\mathrm{sel}}=\{\mathbf{Q}_{\ell}\}_{\ell=1}^{N_Q}` and a
block of :math:`N_t` trajectory frames, pySED evaluates

.. math::

   \rho_w(\mathbf{Q}_{\ell},t_n)
   =
   \sum_{i=1}^{N}
   w_i(\mathbf{Q}_{\ell})
   \exp[
   i\mathbf{Q}_{\ell}\cdot\mathbf{r}_i(t_n)]

and then

.. math::

   S_w(\mathbf{Q}_{\ell},\omega_k)
   =
   \frac{1}{N N_t}
   \left|
   \sum_{n=0}^{N_t-1}
   \rho_w(\mathbf{Q}_{\ell},t_n)
   \exp(-i2\pi kn/N_t)
   \right|^2 .

This is the same finite-time estimator as the pynamic manual's direct
expression for :math:`S(\mathbf{Q},\omega)` [ScatPynamic2025]_, and it is the
efficient production form used by :func:`pySED.dsf.compute_dsf`.

For a selected rectangular mesh in reciprocal-lattice units, pySED uses

.. math::

   H_a = H_{\min}
   +
   \frac{a}{N_H-1}(H_{\max}-H_{\min}),
   \quad
   K_b = K_{\min}
   +
   \frac{b}{N_K-1}(K_{\max}-K_{\min}),
   \quad
   L_c = L_{\min}
   +
   \frac{c}{N_L-1}(L_{\max}-L_{\min}),

and forms

.. math::

   \mathbf{Q}_{abc}
   =
   H_a\mathbf{b}_1
   +
   K_b\mathbf{b}_2
   +
   L_c\mathbf{b}_3 .

For selected paths, segment :math:`s` is sampled as

.. math::

   \mathbf{Q}_{s,m}
   =
   \mathbf{Q}^{\mathrm{start}}_s
   +
   \frac{m}{M_s}
   \left(
   \mathbf{Q}^{\mathrm{end}}_s
   -
   \mathbf{Q}^{\mathrm{start}}_s
   \right),
   \qquad
   m=0,\ldots,M_s-1,

so consecutive path segments do not duplicate shared vertices.  These selected
points are generated by :func:`pySED.q_advisor.qpoints_from_mesh` and
:func:`pySED.q_advisor.qpoints_from_path_segments`:

.. code-block:: python

   from pySED.q_advisor import QAdvisor, qpoints_from_mesh, qpoints_from_path_segments

   mesh = qpoints_from_mesh([-1, 1, 41], [-1, 1, 41], 0.0)
   path = qpoints_from_path_segments(
       starts=[[0, 0, 0], [1, 0, 0]],
       ends=[[1, 0, 0], [1, 1, 0]],
       steps=[40, 40],
   )

   advisor = QAdvisor(primitive_cell, md_supercell)
   report = advisor.advise_experimental_path(mesh.qpoints_reduced, coordinates="reduced")

Only after this generation step does pySED apply the finite-supercell
commensurability policy.  The workflow is therefore:

.. math::

   \mathrm{selected\ experimental\ } \mathbf{Q}
   \rightarrow
   \mathrm{commensurability\ report}
   \rightarrow
   \rho(\mathbf{Q},t)
   \rightarrow
   S(\mathbf{Q},\omega).

The computational cost of the selected-Q direct route is approximately

.. math::

   \mathcal{O}(N_Q N_t N)
   +
   \mathcal{O}(N_Q N_t\log N_t),

where the first term is the phase/density construction and the second term is
the time FFT.  This is more efficient than evaluating every finite-supercell
commensurate point when
:math:`N_Q \ll |\det \mathbf{P}|`, or every point on a dense visualization
mesh when only a line cut or experimental detector window is needed.  It is
not a license to use arbitrary incommensurate Q points: the same periodicity
condition

.. math::

   \exp(i\mathbf{Q}\cdot\mathbf{R}_s)=1

for every supercell lattice vector :math:`\mathbf{R}_s` still applies, giving
the reduced-coordinate rule
:math:`\mathbf{q}_{\mathrm{red}}\mathbf{P}^T\in\mathbb{Z}^3`.

The q-advisor can report this expected saving before a calculation:

.. code-block:: python

   efficiency = advisor.estimate_selected_q_efficiency(
       mesh.qpoints_reduced,
       num_frames=trajectory_num_frames,
       num_atoms=num_atoms,
   )
   print(efficiency.to_table())

The reported ``qpoint_reduction_factor`` is

.. math::

   R_Q
   =
   \frac{|\det \mathbf{P}|}
        {N_Q^{\mathrm{selected}}},

where :math:`|\det\mathbf{P}|` is the number of distinct commensurate q points
in the finite supercell.  This is an upper-bound style Q-count estimate: the
actual wall-time speedup also depends on memory layout, CPU/GPU backend, form
factor evaluation, block size, and I/O.  It is still useful for planning
because the dominant phase construction scales linearly with
:math:`N_Q`.

For neutron, X-ray, and EELS comparison the experimental vector may be an
extended-zone vector,

.. math::

   \mathbf{Q}
   =
   \mathbf{q}
   +
   \mathbf{G}.

The commensurability check is applied to the full reduced
:math:`\mathbf{Q}_{\mathrm{red}}`, not only to the folded first-BZ
:math:`\mathbf{q}`.  This matters because :math:`\mathbf{G}` can change the
atomic form factor and basis-interference phase.  Therefore, when pySED maps a
non-commensurate experimental point to the nearest allowed point, the high-level
DSF workflow preserves the reciprocal-lattice image:

.. math::

   \mathbf{Q}_{\mathrm{red}}=1.26
   \rightarrow
   1.25,
   \qquad
   \mathrm{not}
   \quad
   0.25.

The recommended user entry point is
:func:`pySED.scattering_workflow.compute_dsf_workflow`:

.. code-block:: python

   from pySED.scattering_workflow import compute_dsf_workflow

   result = compute_dsf_workflow(
       positions,
       qpoints=experimental_Q_cartesian,
       primitive_cell=primitive_cell,
       supercell_cell=md_supercell,
       dt=1.0e-12,
       qpoint_coordinates="cartesian",
       q_policy="nearest",       # or "strict"
       max_error_cartesian=0.02,
       atom_types=atom_types,
       experiment="xray",
       map_component="total",    # "coherent" or "incoherent" are also allowed
       output_energy_unit="meV",
       experimental_map=ixs_map,
       sigma_energy=1.5,
   )

``q_policy="strict"`` raises an error for any non-commensurate point.
``q_policy="nearest"`` computes at the nearest commensurate point and stores
the requested-to-used mapping in ``result.q_advice``.  This makes finite-size
approximations visible in the output rather than hiding them inside the DSF.
The same workflow converts the selected DSF component into
``result.scattering_map`` and, when an experimental map is supplied, stores a
resolution-broadened residual comparison in ``result.comparison``.  This puts
DSF on the same paper-ready map path as the EELS and one-phonon visibility
workflows.

Before running an expensive trajectory transform, users can create a standalone
Experiment-Q Advisor report:

.. code-block:: python

   from pySED.q_advisor import QAdvisor

   advisor = QAdvisor(primitive_cell, md_supercell)
   report = advisor.advise_experimental_path(
       experimental_Q_cartesian,
       coordinates="cartesian",
       preserve_image=True,
   )

   print(report.to_table())
   report.write_csv("q_advisor_report.csv")

The report records which experimental points are already commensurate, the
nearest commensurate reduced and Cartesian vectors, the error introduced by
nearest-point mapping, the finite-supercell integer index
:math:`\mathbf{q}_{\mathrm{red}}\mathbf{P}^T`, and a diagonal-supercell
recommendation that would make all requested points exactly commensurate.  This
is the practical guardrail behind pySED's claim that finite-size q sampling is
explicit rather than hidden.

When the experimental data come from a calibrated q-EELS/IXS/INS image, the
image object already contains a scalar q axis but not necessarily the full
three-dimensional reciprocal vector at each pixel column.  pySED handles this
by asking for the vector endpoints of the plotted path and interpolating the
vector :math:`\mathbf{Q}` along the calibrated axis:

.. math::

   \mathbf{Q}_i
   =
   \mathbf{Q}_{\mathrm{start}}
   +
   \frac{s_i-s_0}{s_N-s_0}
   \left(
   \mathbf{Q}_{\mathrm{end}}-\mathbf{Q}_{\mathrm{start}}
   \right),

where :math:`s_i` is the scalar q-axis value stored in the experimental map.
Use :func:`pySED.q_advisor.qpoints_from_path_axis` for this interpolation, or
use the :class:`pySED.q_advisor.QAdvisor` convenience methods directly:

.. code-block:: python

   experimental_map = load_scattering_map_image(
       "experimental_qeels.png",
       q_range=(0.0, 2.0),       # scalar path coordinate on the image
       energy_range=(0.0, 180.0),
   )

   advisor = QAdvisor(primitive_cell, md_supercell)
   q_plan = advisor.plan_map_q_axis(
       experimental_map,
       start_qpoint=Q_start_cartesian,
       end_qpoint=Q_end_cartesian,
       coordinates="cartesian",
       q_policy="nearest",
       max_error_cartesian=0.02,
   )

``q_plan.qpoints_reduced`` are the folded commensurate q points to pass to
phonopy/eigen-SED, while ``q_plan.g_vectors_reduced`` retains the extended-zone
image for EELS, INS, and IXS visibility.

Dynamic Structure Factor
------------------------

The orthodox correlation-function route starts from the Van Hove
density-density correlation function [ScatVanHove1954]_.  For atoms of type
:math:`A`, define the spatial Fourier component of the microscopic density as

.. math::

   \rho_A(\mathbf{Q},t)
   =
   \sum_{i\in A}
   \exp[i\mathbf{Q}\cdot\mathbf{r}_i(t)] .

The partial intermediate scattering function is

.. math::

   F_{AB}(\mathbf{Q},t)
   =
   \frac{1}{\sqrt{N_A N_B}}
   \left<
   \rho_A(\mathbf{Q},t)
   \rho_B(-\mathbf{Q},0)
   \right>,

and the partial dynamic structure factor is its time Fourier transform,

.. math::

   S_{AB}(\mathbf{Q},\omega)
   =
   \frac{1}{2\pi}
   \int_{-\infty}^{\infty}
   F_{AB}(\mathbf{Q},t)
   e^{-i\omega t}\,dt .

For a probe with complex scattering weights :math:`w_A(\mathbf{Q})`, the
coherent probe-specific DSF is obtained by the weighted sum

.. math::

   S_{\mathrm{coh}}(\mathbf{Q},\omega)
   =
   \sum_{A,B}
   w_A(\mathbf{Q})
   w_B^*(\mathbf{Q})
   S_{AB}(\mathbf{Q},\omega).

This is the correlation-function language used by dynasor 2 [ScatDynasor2025]_.
It is the safest theoretical basis for publication because it connects directly
to Van Hove functions, partial correlations, and probe-specific scattering
weights.

Direct Estimator Used in pySED
------------------------------

For efficient computation on selected experimental Q points, pySED also uses
the equivalent direct density-amplitude estimator.  Define

.. math::

   \rho_w(\mathbf{Q},t)
   =
   \sum_i
   w_i(\mathbf{Q})
   \exp[i\mathbf{Q}\cdot\mathbf{r}_i(t)] .

For a finite trajectory block of duration :math:`T`, the coherent spectrum is
estimated as

.. math::

   S_{\mathrm{coh}}(\mathbf{Q},\omega)
   \approx
   \frac{1}{N T}
   \left|
   \int_0^T
   \rho_w(\mathbf{Q},t)e^{-i\omega t}\,dt
   \right|^2,

with block averaging used to estimate uncertainty and reduce noise.  This is
equivalent to Fourier transforming the density-density correlation function
under the usual finite-time and ergodic approximations, but it avoids storing
:math:`F(\mathbf{Q},t)` explicitly.  This is the practical algorithmic idea
used by pynamic-structure-factor [ScatPynamic2025]_.

If the trajectory is split into :math:`N_b` blocks, pySED stores the block
estimate

.. math::

   S_b(\mathbf{Q},\omega)
   =
   \frac{1}{N T_b}
   \left|
   \int_{b}
   \rho_w(\mathbf{Q},t)e^{-i\omega t}\,dt
   \right|^2,

and returns both the block mean

.. math::

   \bar{S}(\mathbf{Q},\omega)
   =
   \frac{1}{N_b}
   \sum_{b=1}^{N_b}
   S_b(\mathbf{Q},\omega)

and the standard error

.. math::

   \mathrm{SEM}[S(\mathbf{Q},\omega)]
   =
   \frac{\mathrm{std}_b[S_b(\mathbf{Q},\omega)]}
   {\sqrt{N_b}}.

The same convention is used for the coherent, incoherent, and total spectra.

The current implementation in :func:`pySED.dsf.compute_dsf` evaluates

.. math::

   \rho_w(\mathbf{Q},t_n)
   =
   \sum_i
   w_i(\mathbf{Q})
   \exp[i\mathbf{Q}\cdot\mathbf{r}_i(t_n)]

for each requested :math:`\mathbf{Q}` and applies a time FFT:

.. math::

   S_{\mathrm{coh}}(\mathbf{Q},\omega_k)
   \propto
   \left|
   \mathrm{FFT}_t[\rho_w(\mathbf{Q},t)]
   \right|^2 .

For validation and small reference calculations,
:func:`pySED.dsf.compute_coherent_dsf_via_correlation` explicitly constructs
the finite-block intermediate scattering function

.. math::

   F_w(\mathbf{Q},m\Delta t)
   =
   \frac{1}{N_t}
   \sum_{n=0}^{N_t-1}
   \rho_w(\mathbf{Q},t_n+m\Delta t)
   \rho_w^*(\mathbf{Q},t_n),

with circular indexing inside each block, and then computes

.. math::

   S_w(\mathbf{Q},\omega_k)
   \propto
   \mathcal{F}_m[
   F_w(\mathbf{Q},m\Delta t)].

This reference path is slower, but it is useful for testing the equivalence
between the correlation-function formula used in the paper and the direct
selected-Q estimator used in production.

For species-resolved interpretation,
:func:`pySED.dsf.compute_partial_dsf` computes the coherent partial spectra
directly from the type-resolved density amplitudes.  With

.. math::

   \rho_A(\mathbf{Q},t)
   =
   \sum_{i\in A}
   \exp[i\mathbf{Q}\cdot\mathbf{r}_i(t)],

the finite-block direct partial estimator is

.. math::

   S_{AB}(\mathbf{Q},\omega_k)
   \propto
   \mathcal{F}_t[\rho_A(\mathbf{Q},t)]_{\omega_k}
   \mathcal{F}_t[\rho_B(\mathbf{Q},t)]_{\omega_k}^{*}.

The probe-weighted coherent total is reconstructed as

.. math::

   S_{\mathrm{coh}}(\mathbf{Q},\omega)
   =
   \sum_{A,B}
   w_A(\mathbf{Q})w_B^*(\mathbf{Q})
   S_{AB}(\mathbf{Q},\omega).

This partial matrix helps identify which atom-type pairs dominate a neutron or
X-ray feature before applying experimental form factors or scattering lengths.
For map-level analysis, the partial result can be sent through the same
paper-ready map path used by total DSF, EELS, and one-phonon visibility:

.. code-block:: python

   from pySED.scattering_workflow import compute_partial_dsf_workflow

   result = compute_partial_dsf_workflow(
       positions,
       qpoints=experimental_qpoints,
       primitive_cell=primitive,
       supercell_cell=supercell,
       dt=dt,
       atom_types=atom_types,
       qpoint_coordinates="cartesian",
       q_policy="nearest",
       species_weights={"Si": 4.1491, "O": 5.803},
       map_component="weighted_total",
   )

To inspect one cross term rather than the weighted total, use
``map_component="partial"` together with ``species_pair=("Si", "O")``.
Because off-diagonal partial spectra are complex, the map converter requires a
choice of ``complex_part="real"``, ``"imag"``, or ``"abs"``.  The default
paper-facing weighted total remains real.

The incoherent neutron contribution is estimated from the self terms,

.. math::

   S_{\mathrm{inc}}(\mathbf{Q},\omega)
   \approx
   \frac{1}{N T}
   \sum_i
   \sigma_i^{\mathrm{inc}}
   \left|
   \int_0^T
   e^{i\mathbf{Q}\cdot\mathbf{r}_i(t)}
   e^{-i\omega t}\,dt
   \right|^2 .

For neutron calculations, :math:`w_i` is the coherent scattering length and
:math:`\sigma_i^{\mathrm{inc}}` is the incoherent cross section.  For X-ray
calculations, pySED evaluates the neutral-atom Cromer-Mann form factor
[ScatCromerMann1968]_

.. math::

   f_i^0(s)
   =
   \sum_{m=1}^{4}
   a_{im}\exp[-b_{im}s^2] + c_i,
   \qquad
   s=\frac{|\mathbf{Q}|}{4\pi},

so the coherent X-ray weight is

.. math::

   w_i(\mathbf{Q}) = f_i^0(|\mathbf{Q}|).

The built-in table covers the same common elements as pySED's default neutron
tables.  Users can override the coefficients for ions, custom tabulations, or
anomalous terms through ``xray_form_factor_table`` or by passing fully
user-defined ``coherent_weights``.  If an element is missing from the
Cromer-Mann table, pySED can fall back to atomic-number weights for exploratory
work, but such fallback should be reported explicitly in publication-quality
X-ray comparisons.

The implemented DSF is classical.  As in classical MD-based DSF workflows, the
raw trajectory spectrum does not automatically satisfy quantum detailed
balance.  pySED therefore provides an optional harmonic
classical-to-quantum correction in :mod:`pySED.quantum_correction`.  With
signed energy transfer :math:`E=\hbar\omega` and
:math:`x=E/(k_B T)`, the correction factor is

.. math::

   C_{\mathrm{q}}(E,T)
   =
   \frac{x}{1-\exp(-x)}.

The experiment-facing estimate is then

.. math::

   S_{\mathrm{q}}(\mathbf{Q},E)
   \approx
   C_{\mathrm{q}}(E,T)
   S_{\mathrm{cl}}(\mathbf{Q},E).

This signed form obeys detailed balance,

.. math::

   \frac{C_{\mathrm{q}}(-E,T)}
   {C_{\mathrm{q}}(E,T)}
   =
   \exp[-E/(k_B T)],

and has the well-defined zero-energy limit
:math:`C_{\mathrm{q}}(0,T)=1`.  When the pySED axis is a frequency in THz,
:mod:`pySED.quantum_correction` uses :math:`E=h\nu`; it can also convert axes
among THz, meV, eV, and :math:`\mathrm{cm}^{-1}`.  This correction should be
reported in any experimental comparison because it changes relative intensities
across energy, not peak positions [ScatSquires2012]_.

Current Correlations
--------------------

The microscopic current density is

.. math::

   \mathbf{j}_w(\mathbf{Q},t)
   =
   \sum_i
   w_i
   \mathbf{v}_i(t)
   \exp[i\mathbf{Q}\cdot\mathbf{r}_i(t)].

With :math:`\hat{\mathbf{Q}}=\mathbf{Q}/|\mathbf{Q}|`, the longitudinal current
is

.. math::

   j_L(\mathbf{Q},t)
   =
   \hat{\mathbf{Q}}\cdot\mathbf{j}_w(\mathbf{Q},t),

and the transverse current is

.. math::

   \mathbf{j}_T(\mathbf{Q},t)
   =
   \mathbf{j}_w(\mathbf{Q},t)
   -
   j_L(\mathbf{Q},t)\hat{\mathbf{Q}}.

The frequency-domain longitudinal and transverse current spectra are useful for
separating longitudinal acoustic/optic features from transverse motion.  In the
continuum relation,

.. math::

   \omega^2 S(\mathbf{Q},\omega)
   =
   Q^2 C_L(\mathbf{Q},\omega),

so current spectra can make dispersive modes clearer when density spectra have
strong elastic backgrounds [ScatDynasor2025]_.

The production current estimator follows the same block-average convention as
the DSF estimator.  :func:`pySED.dsf.compute_current_correlations` returns
``longitudinal_sem`` and ``transverse_sem`` when multiple blocks are used.  For
experiment-facing maps, use
:func:`pySED.scattering_workflow.compute_current_correlation_workflow`:

.. code-block:: python

   from pySED.scattering_workflow import compute_current_correlation_workflow

   result = compute_current_correlation_workflow(
       positions,
       velocities,
       qpoints=experimental_qpoints,
       primitive_cell=primitive,
       supercell_cell=supercell,
       dt=dt,
       qpoint_coordinates="cartesian",
       q_policy="nearest",
       map_component="longitudinal",
       num_blocks=8,
   )

The returned ``result.scattering_map`` can be broadened, normalized, and
compared to an experimental map using the same
:mod:`pySED.compare_experiment` path as total DSF and q-EELS.  This is useful
when the density DSF has a strong elastic line but the longitudinal current
spectrum exposes the inelastic branch more clearly.

Mode-Resolved INS/IXS Visibility
--------------------------------

The correlation-function DSF gives the trajectory spectrum
:math:`S(\mathbf{Q},\omega)`.  To explain which phonon branches are visible in
coherent neutron or X-ray scattering, pySED also provides the harmonic
one-phonon mode visibility in :mod:`pySED.mode_visibility`.  This is an
interpretation layer built from eigenvectors; it does not replace the
finite-temperature MD correlation spectrum.

For a reduced first-BZ wave vector :math:`\mathbf{q}`, reciprocal vector
:math:`\mathbf{G}`, and experimental vector
:math:`\mathbf{Q}=\mathbf{q}+\mathbf{G}`, the elementary contribution from
basis atom :math:`b` and Cartesian direction :math:`\alpha` is

.. math::

   A_{b\alpha,\nu}^{p}(\mathbf{Q})
   =
   w_b^{p}(\mathbf{Q})
   \frac{Q_\alpha e_{b\alpha,\nu}(\mathbf{q})}
   {\sqrt{m_b}}
   \exp[i\mathbf{Q}\cdot\boldsymbol{\tau}_b],

where :math:`p` labels the probe.  For coherent neutron scattering
:math:`w_b^{p}` is the coherent scattering length.  For X-ray scattering it is
the q-dependent atomic form factor :math:`f_b^0(\mathbf{Q})`.  pySED then forms

.. math::

   I_{\nu}^{p}(\mathbf{Q})
   =
   \left|
   \sum_{b,\alpha}
   A_{b\alpha,\nu}^{p}(\mathbf{Q})
   \right|^2.

When desired, the common harmonic one-phonon prefactor can be included as
:math:`I_\nu \rightarrow I_\nu/\omega_{\mathbf{q}\nu}` by using
``frequency_power=-1``.  The default is ``frequency_power=0`` so users first see
the pure polarization, mass, form-factor, and basis-interference visibility.

As in the EELS diagnostic, pySED stores the complex amplitude before the final
squared magnitude.  It reports atom-grouped amplitudes
:math:`A_{b,\nu}=\sum_\alpha A_{b\alpha,\nu}` and direction-grouped amplitudes
:math:`A_{\alpha,\nu}=\sum_b A_{b\alpha,\nu}`.  The atom interference residual

.. math::

   I_\nu^{\mathrm{interf}}
   =
   |A_\nu|^2
   -
   \sum_b |A_{b,\nu}|^2

identifies constructive interference when positive and extinction when
negative.  This is the branch-level explanation layer for deciding whether a
mode that is present in SED is hidden in INS/IXS by
:math:`\mathbf{Q}\cdot\mathbf{e}`, a weak probe weight, mass weighting, or
basis interference.

For paper-ready branch labels, :mod:`pySED.visibility_diagnostics` evaluates
the same quantities as normalized diagnostic scores.  In addition to
:math:`I_\nu^{\mathrm{interf}}`, it records the Cartesian-direction residual

.. math::

   I_\nu^{\mathrm{dir\ interf}}
   =
   |A_\nu|^2
   -
   \sum_\alpha |A_{\alpha,\nu}|^2,

and the normalized fractions

.. math::

   \eta_\nu^{\mathrm{basis}}
   =
   \frac{I_\nu^{\mathrm{interf}}}
        {\sum_b |A_{b,\nu}|^2},
   \qquad
   \eta_\nu^{\mathrm{dir}}
   =
   \frac{I_\nu^{\mathrm{dir\ interf}}}
        {\sum_\alpha |A_{\alpha,\nu}|^2}.

Values near :math:`-1` indicate near-complete extinction by destructive
interference.  The diagnostic table also stores the dominant atom, dominant
Cartesian direction, relative branch visibility at each :math:`(q,G)`, and any
finite-size q-mapping error from :mod:`pySED.q_advisor`.

.. code-block:: python

   from pySED.visibility_diagnostics import diagnose_mode_visibility

   diagnostics = diagnose_mode_visibility(
       result.visibility,
       q_advice=result.q_advice,
       atom_labels=atom_types,
   )
   diagnostics.write_csv("one_phonon_visibility_diagnostics.csv", include_all=False)

The energy-resolved coherent one-phonon map is obtained by combining this
visibility with a branch-resolved spectral function,

.. math::

   I_{p}(\mathbf{Q},\omega)
   =
   \sum_\nu
   I_{\nu}^{p}(\mathbf{Q})
   A_\nu(\mathbf{q},\omega).

As for EELS, :math:`A_\nu` can be a harmonic broadened line shape around the
phonon frequency, or the finite-temperature branch SED
:math:`\Phi_\nu(\mathbf{q},\omega)` from :mod:`pySED.eigen_sed`.  The latter is
implemented by
:func:`pySED.mode_visibility.build_one_phonon_map_from_mode_spectra` and the
high-level workflow
:func:`pySED.scattering_workflow.compute_one_phonon_workflow`:

.. code-block:: python

   from pySED.scattering_workflow import compute_one_phonon_workflow

   result = compute_one_phonon_workflow(
       velocities,
       unitcell_vectors,
       basis_index,
       masses,
       primitive_cell,
       supercell_cell,
       phonopy_eigenvectors,
       basis_positions_reduced,
       g_vectors_reduced,
       dt,
       atom_types=atom_types,
       experiment="xray",      # or "neutron"
       frequency_power=-1,
   )

The returned ``result.visibility`` gives the mode-level explanation, while
``result.scattering_map`` can be passed directly to the experimental comparison
and plotting utilities.

For an experimental extended-zone INS/IXS path, use the q-plan wrapper:

.. code-block:: python

   from pySED.scattering_workflow import compute_one_phonon_workflow_from_q_path

   result = compute_one_phonon_workflow_from_q_path(
       velocities=velocities,
       unitcell_vectors=unitcell_vectors,
       basis_index=basis_index,
       masses=masses,
       primitive_cell=primitive_cell,
       supercell_cell=md_supercell,
       eigenvectors=phonopy_eigenvectors,
       basis_positions_reduced=basis_positions_reduced,
       experimental_qpoints=experimental_Q_cartesian,
       dt=dt,
       qpoint_coordinates="cartesian",
       q_policy="nearest",
       max_error_cartesian=0.02,
       atom_types=atom_types,
       experiment="xray",      # or "neutron"
       frequency_power=-1,
   )

This follows the same :math:`Q=q+G` rule as q-EELS: phonopy eigenvectors and
eigen-SED are evaluated at ``result.q_plan.qpoints_reduced``, while
``result.q_plan.g_vectors_reduced`` is retained in the one-phonon scattering
visibility.

Eigenvector-Resolved SED
------------------------

The existing pySED production workflow uses the eigenvector-free SED formula
described in the main theory chapter and in the SED literature
[ScatThomas2010]_.  For branch-level interpretation, :mod:`pySED.eigen_sed`
adds a phonopy eigenvector projection.  For primitive-cell basis atom
:math:`b`, repeated cell :math:`l`, and Cartesian direction :math:`\alpha`,

.. math::

   \dot{Q}_{\mathbf{q}\nu}(t)
   =
   \sum_{l,b,\alpha}
   \sqrt{m_b}\,
   v_{lb\alpha}(t)
   e^*_{b\alpha}(\mathbf{q},\nu)
   \exp[i\mathbf{q}\cdot\mathbf{R}_l].

The branch-resolved SED is

.. math::

   \Phi_{\nu}(\mathbf{q},\omega)
   =
   \left|
   \mathcal{F}_t[
   \dot{Q}_{\mathbf{q}\nu}(t)]
   \right|^2 .

Even when eigenvectors are read from phonopy, the q points must still satisfy
the finite MD supercell condition

.. math::

   \mathbf{q}_{\mathrm{red}}\mathbf{P}^{T} \in \mathbb{Z}^3.

Phonopy can compute eigenvectors at arbitrary reciprocal-space points because
it evaluates a force-constant dynamical matrix.  MD SED cannot use arbitrary
q points in the same way because the trajectory is finite and periodic.  The
reasonable workflow is therefore:

1. choose experimental or high-symmetry q points,
2. use :mod:`pySED.q_advisor` to keep or map them to commensurate q points,
3. ask phonopy for eigenvectors at exactly those commensurate q points,
4. project the MD velocities with :mod:`pySED.eigen_sed`, and
5. compare branch intensities with DSF/EELS visibility.

This resolves the practical question: eigenvectors do not remove the
commensurability requirement for an MD trajectory.  They add branch resolution
after a valid finite-supercell q point has been selected.

Kinematic q-EELS Visibility
---------------------------

Momentum-resolved vibrational EELS measures an electron-scattering-weighted
phonon spectrum.  The first pySED layer is a kinematic one-phonon visibility
model in the extended Brillouin zone.

Experimental q-EELS axes are often calibrated in detector angle rather than
direct reciprocal-space units.  :mod:`pySED.electron_kinematics` converts this
geometry into Cartesian momentum transfer before the q-advisor step.  For an
incident kinetic energy :math:`E_0`, the relativistic electron momentum is

.. math::

   pc
   =
   \sqrt{E_0(E_0+2m_ec^2)},
   \qquad
   \lambda
   =
   \frac{hc}{pc},
   \qquad
   k_0
   =
   \frac{2\pi}{\lambda}.

For small scattering angles :math:`(\theta_x,\theta_y)` and energy loss
:math:`\Delta E`, pySED uses the q-EELS calibration [ScatEgerton2011]_

.. math::

   Q_x \simeq k_0\theta_x,
   \qquad
   Q_y \simeq k_0\theta_y,
   \qquad
   Q_z \simeq \frac{\Delta E}{\hbar v}
   =
   \frac{\Delta E}{\hbar c\beta},

where :math:`\beta=v/c`.  For vibrational losses the longitudinal term is
usually much smaller than the transverse detector momentum, so the
line-scan helper leaves it off by default:

.. code-block:: python

   from pySED.electron_kinematics import eels_qpoints_from_angle_axis

   experimental_Q_cartesian = eels_qpoints_from_angle_axis(
       angle_axis_mrad,
       beam_energy_ev=200_000.0,
       direction=(1.0, 0.0),
       angle_unit="mrad",
       include_longitudinal=False,
   )

The resulting ``experimental_Q_cartesian`` array can be passed directly to
:class:`pySED.q_advisor.QAdvisor` or to
``compute_eels_workflow_from_q_path(..., qpoint_coordinates="cartesian")``.
The same small-angle conversion can translate an experimental angular
resolution into the momentum-resolution width used for map comparison:

.. code-block:: python

   from pySED.electron_kinematics import angular_resolution_to_q_sigma

   sigma_q = angular_resolution_to_q_sigma(
       sigma_angle=2.0,
       beam_energy_ev=200_000.0,
       angle_unit="mrad",
   )

   comparison = prepare_map_comparison(
       simulated_map,
       experimental_map,
       sigma_q=sigma_q,
       sigma_energy=3.0,
       energy_unit="meV",
   )

For a complete MD-to-q-EELS calculation, pySED also provides the one-call
workflow:

.. code-block:: python

   from pySED.scattering_workflow import compute_eels_workflow_from_angle_axis

   result = compute_eels_workflow_from_angle_axis(
       velocities=velocities,
       unitcell_vectors=unitcell_vectors,
       basis_index=basis_index,
       masses=masses,
       primitive_cell=primitive_cell,
       supercell_cell=md_supercell,
       eigenvectors=phonopy_eigenvectors,
       basis_positions_reduced=basis_positions_reduced,
       angle_axis=angle_axis_mrad,
       beam_energy_ev=200_000.0,
       direction=(1.0, 0.0),
       angle_unit="mrad",
       dt=dt,
       q_policy="nearest",
       max_error_cartesian=0.02,
       atom_types=atom_types,
       electron_form_factor_model="mott-bethe",
   )

For the scattering visibility itself, given

.. math::

   \mathbf{Q} = \mathbf{q} + \mathbf{G},

where :math:`\mathbf{q}` is folded into the first Brillouin zone and
:math:`\mathbf{G}` is a reciprocal lattice vector, pySED estimates the
branch-resolved visibility as

.. math::

   I_{\nu}(\mathbf{Q})
   \propto
   \left|
   \sum_b
   f_b^{e}(\mathbf{Q})
   \frac{\mathbf{Q}\cdot
   \mathbf{e}_{b\nu}(\mathbf{q})}
   {\sqrt{m_b}}
   \exp[i\mathbf{Q}\cdot\boldsymbol{\tau}_b]
   \right|^2 .

This expression contains the main terms needed for interpreting
momentum-resolved EELS selection and interference:

- :math:`\mathbf{Q}\cdot\mathbf{e}_{b\nu}` gives polarization selection,
- :math:`\exp[i\mathbf{Q}\cdot\boldsymbol{\tau}_b]` gives basis interference,
- :math:`f_b^e(\mathbf{Q})` gives electron probe weighting, and
- :math:`\mathbf{Q}=\mathbf{q}+\mathbf{G}` exposes higher-Brillouin-zone
  extinction and enhancement.

By default pySED keeps :math:`f_b^e(\mathbf{Q})=1` to preserve the earlier
unit-weight visibility.  For a first quantitative electron-probe estimate,
``electron_form_factor_model="mott-bethe"`` uses the neutral-atom Mott-Bethe
relation [ScatMottBethe]_ [ScatPeng1996]_

.. math::

   f_b^e(\mathbf{Q})
   =
   C_e
   \frac{Z_b-f_b^X(|\mathbf{Q}|)}
        {|\mathbf{Q}|^2},

where :math:`f_b^X` is the Cromer-Mann X-ray form factor
[ScatCromerMann1968]_ and :math:`Z_b` is the atomic number.  The current
kinematic EELS implementation drops the constant :math:`C_e` because relative
intensities are sufficient for branch visibility, extinction, and map
comparison after normalization.  For :math:`|\mathbf{Q}|\rightarrow0`, pySED
uses the analytic derivative of the Cromer-Mann expansion,

.. math::

   \lim_{Q\rightarrow0}
   \frac{Z-f^X(Q)}{Q^2}
   =
   \frac{1}{16\pi^2}
   \sum_i a_i b_i,

provided the neutral-atom coefficients satisfy
:math:`f^X(0)=\sum_i a_i+c\simeq Z`.  Users can still pass
``electron_form_factors`` explicitly when they have a dedicated electron
scattering table, ionic form factors, or a later multislice/TACAW model.

For a real q-EELS line scan, the measured points are usually a pairwise path
:math:`\mathbf{Q}_i`, not the Cartesian product of all folded q points and all
reciprocal-lattice vectors.  pySED therefore provides an explicit extended-zone
planner:

.. code-block:: python

   from pySED.q_advisor import QAdvisor

   advisor = QAdvisor(primitive_cell, md_supercell)
   q_plan = advisor.plan_extended_zone_q_path(
       experimental_Q_cartesian,
       coordinates="cartesian",
       q_policy="nearest",
       max_error_cartesian=0.02,
   )

   q_plan.write_phonopy_qpoints("phonopy_qpoints.dat")

The plan stores

.. math::

   \mathbf{Q}_i^{\mathrm{used}}
   =
   \mathbf{q}_i^{\mathrm{folded}}
   +
   \mathbf{G}_i,

where :math:`\mathbf{Q}_i^{\mathrm{used}}` is commensurate with the MD
supercell, :math:`\mathbf{q}_i^{\mathrm{folded}}` is the q point that should be
used for phonopy eigenvectors and eigen-SED, and :math:`\mathbf{G}_i` is the
integer extended-zone vector.  Use
:func:`pySED.eels.compute_mode_visibility_for_q_path`, or
``compute_eels_workflow(..., pairwise_g_vectors=True)``, when
``qpoints_reduced`` and ``g_vectors_reduced`` have this one-to-one path
meaning.  The default :func:`pySED.eels.compute_mode_visibility` behavior still
computes the full q-by-G product, which is useful for making extended-zone maps
over several Brillouin zones.

For production scripts, the recommended single entry point is
:func:`pySED.scattering_workflow.compute_eels_workflow_from_q_path`:

.. code-block:: python

   from pySED.scattering_workflow import compute_eels_workflow_from_q_path

   result = compute_eels_workflow_from_q_path(
       velocities=velocities,
       unitcell_vectors=unitcell_vectors,
       basis_index=basis_index,
       masses=masses,
       primitive_cell=primitive_cell,
       supercell_cell=md_supercell,
       eigenvectors=phonopy_eigenvectors,
       basis_positions_reduced=basis_positions_reduced,
       experimental_qpoints=experimental_Q_cartesian,
       dt=dt,
       qpoint_coordinates="cartesian",
       q_policy="nearest",
       max_error_cartesian=0.02,
       atom_types=atom_types,
       electron_form_factor_model="mott-bethe",
   )

This wrapper returns ``result.q_plan`` together with ``result.eigen_sed``,
``result.visibility``, ``result.eels_map``, and ``result.scattering_map``.  It
also validates that the phonopy eigenvectors are supplied at the folded
commensurate q points from ``result.q_plan.qpoints_reduced``.  In practice this
prevents a common mistake: using the extended-zone experimental
:math:`\mathbf{Q}` as the phonopy q point instead of using the folded
:math:`\mathbf{q}` and retaining :math:`\mathbf{G}` only in the scattering
visibility.

For mode-level interpretation,
:func:`pySED.eels.compute_mode_visibility_decomposition` keeps the complex
amplitude before the final squared magnitude.  The elementary contribution is

.. math::

   A_{b\alpha,\nu}(\mathbf{Q})
   =
   f_b^e(\mathbf{Q})
   \frac{Q_\alpha e_{b\alpha,\nu}(\mathbf{q})}
   {\sqrt{m_b}}
   \exp[i\mathbf{Q}\cdot\boldsymbol{\tau}_b].

The total mode amplitude and visibility are

.. math::

   A_\nu(\mathbf{Q})
   =
   \sum_{b,\alpha}
   A_{b\alpha,\nu}(\mathbf{Q}),
   \qquad
   I_\nu(\mathbf{Q})
   =
   |A_\nu(\mathbf{Q})|^2.

pySED also reports coherent diagnostic amplitudes grouped by atom and by
Cartesian direction,

.. math::

   A_{b,\nu}=\sum_\alpha A_{b\alpha,\nu},
   \qquad
   A_{\alpha,\nu}=\sum_b A_{b\alpha,\nu}.

The diagnostic quantities :math:`|A_{b,\nu}|^2` and
:math:`|A_{\alpha,\nu}|^2` identify which basis atoms and displacement
directions are visible before interference.  They are not additive intensities:

.. math::

   |A_\nu|^2
   \ne
   \sum_b |A_{b,\nu}|^2

when basis interference is present.  pySED therefore records the atom
interference residual

.. math::

   I_\nu^{\mathrm{interf}}
   =
   |A_\nu|^2
   -
   \sum_b |A_{b,\nu}|^2,

which is positive for constructive interference and negative for extinction.
This is the practical diagnostic for deciding whether a branch disappears
because of polarization selection, a weak form factor, or destructive basis
interference.

The same :mod:`pySED.visibility_diagnostics` helper can be applied to the EELS
decomposition object.  For an EELS workflow, use
``compute_mode_visibility_decomposition`` when calling the low-level API, or
pass ``result.visibility`` from a workflow when it already contains atom and
direction amplitudes.  The resulting CSV table gives a compact branch-by-branch
classification: ``visible``, ``basis_interference``,
``polarization_or_zero_weight``, ``cartesian_interference``,
``finite_size_q_mapping``, or ``weak_visibility``.  These labels are diagnostic
metadata; the simulated intensity still comes from the coherent visibility
:math:`I_\nu(\mathbf{Q})`.

The energy-resolved kinematic map is then

.. math::

   I_{\mathrm{EELS}}(\mathbf{Q},\omega)
   =
   \sum_{\nu}
   I_{\nu}(\mathbf{Q})
   A_{\nu}(\mathbf{q},\omega),

where :math:`A_{\nu}` can be supplied in two ways.  For a harmonic reference,
:func:`pySED.eels.build_eels_map` represents :math:`A_{\nu}` by a Gaussian or
Lorentzian broadening around the phonon frequency.  For a finite-temperature MD
calculation, :func:`pySED.eels.build_eels_map_from_mode_spectra` can consume the
branch-resolved spectral function from :mod:`pySED.eigen_sed`,

.. math::

   A_{\nu}(\mathbf{q},\omega)
   \equiv
   \Phi_{\nu}(\mathbf{q},\omega),

and evaluate the same visibility-weighted sum directly.  This is the preferred
pySED-native route for a MD-to-EELS comparison because anharmonic peak shifts
and linewidths are inherited from the trajectory rather than imposed by an
external broadening model.

The same branch spectra can be used to make energy-resolved diagnostic maps:

.. math::

   I_b(\mathbf{Q},\omega)
   =
   \sum_\nu
   |A_{b,\nu}(\mathbf{Q})|^2
   A_\nu(\mathbf{q},\omega),
   \qquad
   I_{\mathrm{interf}}(\mathbf{Q},\omega)
   =
   \sum_\nu
   I_{\nu}^{\mathrm{interf}}(\mathbf{Q})
   A_\nu(\mathbf{q},\omega).

They obey

.. math::

   I_{\mathrm{EELS}}(\mathbf{Q},\omega)
   =
   \sum_b I_b(\mathbf{Q},\omega)
   +
   I_{\mathrm{interf}}(\mathbf{Q},\omega),

while the Cartesian-direction maps remain diagnostic projections rather than an
additive partition of the total intensity.  Use
:func:`pySED.eels.build_eels_decomposition_map_from_mode_spectra` for EELS and
:func:`pySED.mode_visibility.build_one_phonon_decomposition_map_from_mode_spectra`
for coherent neutron/X-ray one-phonon maps.

A later PySlice-like layer can replace this kinematic intensity with a
pySED-native multislice or TACAW electron propagation model [ScatPySlice2026]_.
The current layer is still valuable because it already explains many missing
or enhanced branches through polarization, basis interference, form factors,
and finite-size q sampling.

Experiment Comparison
---------------------

The :mod:`pySED.compare_experiment` module provides a paper-ready comparison
path:

1. load an experimental intensity map with calibrated q and energy axes,
2. optionally convert simulation energy axes to the experimental unit,
3. optionally apply the harmonic quantum detailed-balance correction,
4. broaden simulated spectra with the experimental energy resolution,
5. convolve or smooth along q to model momentum resolution,
6. normalize simulated and experimental maps consistently,
7. compute residual maps and scalar residual metrics,
8. extract line cuts at selected q or energy values, and
9. track simulated and experimental peak positions.

Experimental figures can be imported directly when only an image is available.
For a calibrated q-EELS image, use
:func:`pySED.compare_experiment.load_scattering_map_image` to convert pixels
into a :class:`pySED.compare_experiment.ScatteringMap`:

.. code-block:: python

   from pySED.compare_experiment import load_scattering_map_image

   experimental_map = load_scattering_map_image(
       "experimental_qeels.png",
       q_range=(-0.6, 0.6),       # left and right image limits
       energy_range=(0.0, 180.0), # bottom and top image limits
       crop=(40, 420, 30, 730),   # row_start, row_stop, col_start, col_stop
       origin="upper",
       invert=False,
       normalization="max",
   )

The image columns are mapped to the q axis and the image rows are mapped to the
energy axis.  ``origin="upper"`` is appropriate for ordinary image files where
row zero is the top of the image; pySED flips the pixel array so the returned
energy axis increases from ``energy_range[0]`` to ``energy_range[1]``.  If the
experimental map uses dark intensity on a bright background, set ``invert=True``.

For direct figure export, :mod:`pySED.scattering_plot` can render individual
maps, overlaid line cuts, or a combined simulated/experimental/residual figure.
The combined figure uses the normalized arrays stored in
:class:`pySED.compare_experiment.MapComparison`, so the plotted residual is
exactly the quantity used for RMSE, MAE, and correlation:

.. code-block:: python

   from pySED.compare_experiment import prepare_map_comparison
   from pySED.scattering_plot import plot_map_comparison

   comparison = prepare_map_comparison(
       simulated_map,
       experimental_map,
       sigma_q=0.02,
       sigma_energy=3.0,
       normalization="max",
       quantum_temperature=300.0,
       energy_unit="meV",
   )
   fig = plot_map_comparison(
       comparison,
       q_label="Q path (1/Angstrom)",
       energy_label="Energy loss (meV)",
       save_path="pysed_eels_comparison.png",
   )

The same comparison object can be exported as manuscript-ready tables:

.. code-block:: python

   from pySED.compare_experiment import write_map_comparison_report

   paths = write_map_comparison_report(
       comparison,
       output_dir="paper_eels_comparison",
       prefix="figure_3",
       linecut_q_indices=[10, 25],
       linecut_energy_indices=[40],
   )

This writes ``figure_3_summary.json`` with RMSE, MAE, correlation, axis ranges,
and metadata; ``figure_3_residual.csv`` with flattened
:math:`R(q,\omega)` values; and ``figure_3_peaks.csv`` with the simulated and
experimental peak energy at each q point.  When line-cut indices are supplied,
it also writes ``figure_3_linecuts.csv``.  Fixed-q rows report intensity as a
function of energy, while fixed-energy rows report intensity along the q path.
These exported tables are intended to accompany the plotted residual map and
line cuts in a reproducible paper workflow.

For a complete workflow result, :mod:`pySED.scattering_export` can write a
single reproducibility bundle:

.. code-block:: python

   from pySED.scattering_export import write_scattering_workflow_bundle

   written = write_scattering_workflow_bundle(
       result,
       output_dir="figure_3_bundle",
       prefix="qeels",
       atom_labels=atom_types,
       linecut_q_indices=[10, 25],
       linecut_energy_indices=[40],
   )

For q-EELS and one-phonon workflows this bundle includes the extended-zone
``q_plan`` CSV, folded ``phonopy_qpoints.dat``, branch visibility diagnostics,
the simulated scattering map, and any experimental comparison summary,
residual, peak, and line-cut tables.  For DSF workflows it exports the q-advisor
table and scattering map.  This is the recommended archive format for
paper-figure provenance because the finite-size q mapping and the residual
metrics are stored beside the simulated intensity.

For a two-dimensional q-energy map, the high-level workflow is represented by
:class:`pySED.compare_experiment.ScatteringMap` and
:func:`pySED.compare_experiment.prepare_map_comparison`.  If the simulated map
has axes :math:`(q_s,\omega_s)` and the experimental map has axes
:math:`(q_e,\omega_e)`, pySED first interpolates

.. math::

   I_s(q_s,\omega_s)
   \rightarrow
   I_s(q_e,\omega_e),

then can apply the optional detailed-balance correction
:math:`C_{\mathrm{q}}(E,T)`.  The experimental resolution is then modeled as a
Gaussian convolution,

.. math::

   \tilde{I}_s(q,\omega)
   =
   G_q(\sigma_q) * G_\omega(\sigma_\omega) *
   I_s(q,\omega),

where :math:`\sigma_q` and :math:`\sigma_\omega` are supplied in the physical
axis units rather than pixel units.  The residual map is

.. math::

   R(q,\omega)
   =
   N[\tilde{I}_s(q,\omega)]
   -
   N[I_{\mathrm{exp}}(q,\omega)],

where :math:`N[\cdot]` is the selected normalization, for example max, area, or
z-score normalization.  pySED reports RMSE, MAE, and the Pearson correlation of
the aligned maps.

The intended manuscript-level workflow is therefore:

.. math::

   \mathrm{MD}
   \rightarrow
   \mathrm{commensurate\ } \mathbf{Q}
   \rightarrow
   \mathrm{SED/DSF/EELS}
   \rightarrow
   \mathrm{resolution\ convolution}
   \rightarrow
   \mathrm{experimental\ comparison}.

Comparison With dynasor and pynamic
-----------------------------------

The recommended interpretation is:

- dynasor is the main theoretical and validation reference.  It formalizes the
  DSF through Van Hove and intermediate scattering functions, supports partial
  correlations, probe-specific weights, and mode projection.
- pynamic-structure-factor is a practical implementation reference for the
  direct selected-Q calculation of inelastic neutron and X-ray
  :math:`S(\mathbf{Q},\omega)`.
- pySED should combine these strengths: dynasor-style theory, pynamic-style
  selected-Q efficiency, and pySED-specific commensurate q selection, SED
  integration, EELS extended-zone visibility, and experiment comparison.

This is why the current implementation uses

.. math::

   F_{AB}(\mathbf{Q},t)
   \xrightarrow{\mathcal{F}_t}
   S_{AB}(\mathbf{Q},\omega)

as the paper-facing theory, while the code evaluates

.. math::

   \rho_w(\mathbf{Q},t)
   \xrightarrow{\mathrm{FFT}_t}
   S_w(\mathbf{Q},\omega)

as the efficient estimator.

Validation Targets
------------------

The current unit tests cover array shapes, non-negative spectra, q-point
validation, nearest-q mapping, EELS visibility, eigenvector SED, and experiment
map utilities.  The module :mod:`pySED.validation` and the script
``benchmarks/dsf_reference_validation.py`` provide a lightweight developer
entry point for checking the production direct estimator against the explicit
correlation reference estimator:

.. code-block:: powershell

   python benchmarks\dsf_reference_validation.py `
       --json validation_report.json `
       --export-case validation_case.npz

The JSON report records Python, NumPy, and pySED versions, the internal
direct-vs-correlation error, and whether optional external packages such as
dynasor or pynamic-structure-factor are importable.  The exported ``.npz`` file
contains the deterministic positions, Q points, q-dependent coherent weights,
time step, block count, and pySED reference spectra.  This makes external
cross-validation reproducible: dynasor or pynamic should be run on the same
trajectory, Q points, weights, time window, block convention, and normalization
before comparing the spectra.

The external comparison path is intentionally file-based so that dynasor or
pynamic can be run in a separate environment without becoming pySED runtime
dependencies.  A reproducible comparison has three steps:

.. code-block:: powershell

   python benchmarks\dsf_reference_validation.py `
       --export-case validation_case.npz

Run the external package on ``validation_case.npz`` using the same
``positions``, ``qpoints_cartesian``, ``coherent_weights``, ``dt``, and
``num_blocks``.  Save the external spectrum to an ``.npz`` file with a spectrum
key such as ``coherent`` and a frequency key such as ``frequencies_thz``.  The
frequency unit must match pySED's THz axis; convert angular-frequency outputs
before writing the comparison file.

.. code-block:: powershell

   python benchmarks\dsf_reference_validation.py `
       --compare-reference external_dsf.npz `
       --reference-key coherent `
       --reference-frequency-key frequencies_thz `
       --candidate-component coherent `
       --normalization least_squares `
       --comparison-json external_compare.json

``--normalization least_squares`` reports the best scalar scale factor before
computing RMSE, which is useful when two packages use different Fourier or
probe-weight prefactors.  Use ``--normalization none`` for absolute-intensity
checks after the conventions have been matched.

The internal reference check is not a substitute for external validation.  It
proves only the finite-time identity

.. math::

   \left|
   \mathcal{F}_t[\rho_w(\mathbf{Q},t)]
   \right|^2
   \equiv
   \mathcal{F}_t[
   \langle \rho_w(\mathbf{Q},t+\tau)
   \rho_w^*(\mathbf{Q},t)\rangle_t]

under the same circular-block convention.  External validation is still needed
to compare normalization, partial correlations, weighting conventions, and
trajectory handling against established packages.  The remaining scientific
validation targets are:

- compare pySED DSF with dynasor for identical trajectories, Q points, weights,
  and normalization conventions;
- compare the direct :math:`\rho(\mathbf{Q},t)\rightarrow\mathrm{FFT}`
  estimator with pynamic-structure-factor for identical trajectories;
- verify that non-commensurate q points are rejected or mapped before
  eigenvector projection;
- validate extended-zone EELS selection and interference against q-EELS
  examples and PySlice benchmarks; and
- add experimental map examples with residuals, line cuts, and peak tracking.

References
----------

.. [ScatVanHove1954] L. Van Hove, "Correlations in space and time and Born
   approximation scattering in systems of interacting particles," *Physical
   Review* **95**, 249 (1954). https://doi.org/10.1103/PhysRev.95.249

.. [ScatSquires2012] G. L. Squires, *Introduction to the Theory of Thermal Neutron
   Scattering*, Dover Publications (2012).

.. [ScatThomas2010] J. A. Thomas, J. E. Turney, R. M. Iutzi, C. H. Amon, and
   A. J. H. McGaughey, "Predicting phonon dispersion relations and lifetimes
   from the spectral energy density," *Physical Review B* **81**, 081411
   (2010). https://doi.org/10.1103/PhysRevB.81.081411

.. [ScatCromerMann1968] D. T. Cromer and J. B. Mann, "X-ray scattering factors
   computed from numerical Hartree-Fock wave functions," *Acta
   Crystallographica Section A* **24**, 321-324 (1968).
   https://doi.org/10.1107/S0567739468000550

.. [ScatMottBethe] L. Pacoste, V. M. Ignat'ev, P. M. Dominiak, and X. Zou,
   "On the structure refinement of metal complexes against 3D electron
   diffraction data using multipolar scattering factors," *IUCrJ* **11**,
   878-890 (2024). https://doi.org/10.1107/S2052252524006730

.. [ScatPeng1996] L.-M. Peng, G. Ren, S. L. Dudarev, and M. J. Whelan,
   "Robust parameterization of elastic and absorptive electron atomic
   scattering factors," *Acta Crystallographica Section A* **52**, 257-276
   (1996). https://doi.org/10.1107/S0108767395014371

.. [ScatEgerton2011] R. F. Egerton, *Electron Energy-Loss Spectroscopy in the
   Electron Microscope*, 3rd ed. (Springer, 2011).
   https://doi.org/10.1007/978-1-4419-9583-4

.. [ScatDynasor2025] E. Berger, E. Fransson, F. Eriksson, E. Lindgren,
   G. Wahnstrom, T. H. Rod, and P. Erhart, "Dynasor 2: From simulation to
   experiment through correlation functions," *Computer Physics
   Communications* **316**, 109759 (2025).
   https://doi.org/10.1016/j.cpc.2025.109759

.. [ScatPynamic2025] T. C. Sterling, "pynamic-structure-factor: a python program
   to calculate dynamic structure factors in the classical limit," manual dated
   May 8, 2025, and source repository.
   https://github.com/tyst3273/pynamic-structure-factor

.. [ScatPySlice2026] H. Walker, "PySlice: Routine Vibrational Electron Energy Loss
   Spectroscopy Prediction with Universal Interatomic Potentials" (2026).
   https://github.com/h-walk/PySlice
