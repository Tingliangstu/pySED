Theory
======

This chapter summarizes the theory used by **pySED**. The notation follows the
pySED paper, "PYSED: A tool for extracting kinetic-energy-weighted phonon
dispersion and lifetime from molecular dynamics simulations", J. Appl. Phys.
**138**, 075101 (2025) [Liang2025]_.

What pySED Computes
-------------------

The spectral energy density (SED) method projects atomic velocities from a
molecular dynamics trajectory into reciprocal space. The resulting intensity
map, :math:`\Phi'(\mathbf{q}, \omega)`, shows where kinetic energy is
concentrated as a function of wave vector :math:`\mathbf{q}` and angular
frequency :math:`\omega`.

The SED map can be used to:

- identify phonon dispersion branches directly from MD trajectories,
- compare finite-temperature MD results with lattice-dynamics calculations,
- inspect anharmonic broadening and temperature effects,
- estimate phonon lifetimes by fitting individual SED peaks.

The lifetime values obtained from SED peak fitting should normally be treated as
qualitative or semi-quantitative. For highly accurate thermal conductivity,
methods such as HNEMD, EMD, NEMD, or appropriate Boltzmann transport workflows
are usually more suitable.

Normal-Mode SED
---------------

For a periodic supercell, the velocity of basis atom :math:`b` in unit cell
:math:`l` and Cartesian direction :math:`\alpha` is written as
:math:`\dot{u}_{\alpha}(l,b,t)`. The corresponding normal-mode velocity is

.. math::

   \dot{q}(\mathbf{q},\nu,t)
   =
   \sum_{\alpha,b,l}^{3,n,N_T}
   \sqrt{\frac{m_b}{N_T}}
   \dot{u}_{\alpha}(l,b,t)
   e_{\alpha}^{*}(\mathbf{q},\nu,b)
   \exp\left(i\mathbf{q}\cdot\mathbf{r}_0(l)\right).

Here :math:`m_b` is the mass of basis atom :math:`b`,
:math:`e_{\alpha}^{*}(\mathbf{q},\nu,b)` is the complex conjugate of the
phonon eigenvector component, :math:`\nu` is the phonon branch index,
:math:`N_T` is the number of unit cells, and :math:`\mathbf{r}_0(l)` is the
equilibrium position of unit cell :math:`l`.

The time-averaged kinetic energy of a normal mode is

.. math::

   T(\mathbf{q},\nu)
   =
   \frac{1}{2\tau_0}
   \int_0^{\tau_0}
   \left|\dot{q}(\mathbf{q},\nu,t)\right|^2 dt,

where :math:`\tau_0` is the trajectory length used in the SED calculation.
Applying the Fourier transform and Parseval's theorem gives an energy
distribution in the frequency domain.

The eigenvector-based SED is then

.. math::

   \Phi(\mathbf{q},\omega)
   =
   \frac{1}{4\pi\tau_0 N_T}
   \sum_{\nu}^{3n}
   \left|
   \sum_{\alpha}^{3}
   \sum_b^n
   \int_0^{\tau_0}
   \sum_l^{N_T}
   \sqrt{m_b}\,
   \dot{u}_{\alpha}(l,b,t)
   e_{\alpha}^{*}(\mathbf{q},\nu,b)
   \exp\left(i\mathbf{q}\cdot\mathbf{r}_0(l)-i\omega t\right)dt
   \right|^2 .

This expression resolves individual branches through the eigenvectors, but it
requires a lattice-dynamics calculation and can be expensive for complex
systems.

Eigenvector-Free SED Used by pySED
----------------------------------

pySED currently implements the eigenvector-free SED expression derived by
Thomas et al. [Thomas2010]_. This expression sums over the basis atoms and Cartesian
directions without explicitly projecting onto phonon eigenvectors:

.. math::

   \Phi'(\mathbf{q},\omega)
   =
   \frac{1}{4\pi\tau_0 N_T}
   \sum_{\alpha}^{3}
   \sum_b^n
   m_b
   \left|
   \int_0^{\tau_0}
   \sum_l^{N_T}
   \dot{u}_{\alpha}(l,b,t)
   \exp\left(i\mathbf{q}\cdot\mathbf{r}_0(l)-i\omega t\right)dt
   \right|^2 .

The equivalence between the eigenvector-based and eigenvector-free total SED
follows from the eigenvector orthonormality condition

.. math::

   \sum_{\nu}^{3n}
   e_{\alpha}(\mathbf{q},\nu,b)
   e_{\beta}^{*}(\mathbf{q},\nu,b)
   =
   \delta_{\alpha\beta}.

In practice, :math:`\Phi'(\mathbf{q},\omega)` is usually sufficient for
obtaining a kinetic-energy-weighted phonon dispersion from MD trajectories. A
branch-resolved eigenvector projection can help separate nearly degenerate
modes, but that is not the current production workflow in pySED.

Partial SED
-----------

The partial SED implementation decomposes the total pySED intensity by atom
species and Cartesian direction. Let :math:`\mathcal{B}_s` be the set of basis
atoms that belong to species or type :math:`s`. The contribution from species
:math:`s` and direction :math:`\alpha` is

.. math::

   \Phi'_{s,\alpha}(\mathbf{q},\omega)
   =
   \frac{1}{4\pi\tau_0 N_T}
   \sum_{b\in\mathcal{B}_s}
   m_b
   \left|
   \int_0^{\tau_0}
   \sum_{l=1}^{N_T}
   \dot{u}_{\alpha}(l,b,t)
   \exp\left(i\mathbf{q}\cdot\mathbf{r}_0(l)-i\omega t\right)dt
   \right|^2 .

The total SED is recovered by summing over all species and directions:

.. math::

   \Phi'_{\mathrm{total}}(\mathbf{q},\omega)
   =
   \sum_s \sum_{\alpha=x,y,z}
   \Phi'_{s,\alpha}(\mathbf{q},\omega).

In the input file, use ``output_partial = 1`` during the SED-computing mode.
pySED writes one partial SED file for each atom type and direction. In plotting
mode, use ``plot_partial_SED = 3`` to plot the sum of ``x+y+z`` for type 3, or
``plot_partial_SED = 3 x`` to plot only the x-direction contribution of type 3.

Partial SED is not the same as local or spatial-bin SED. The standard SED method
requires a mapping between the supercell and the primitive cell to define
reciprocal-space wave vectors. pySED therefore does not currently support
arbitrary local SED extraction for user-defined spatial bins.

Allowed Wave Vectors
--------------------

In periodic systems, only wave vectors commensurate with the MD supercell are
allowed. pySED follows the same commensurate-supercell q-point idea used in
dynasor 2 [Dynasor2025]_, and the implementation in ``pySED/construct_BZ.py``
credits dynasor's lattice helper as the reference for this path construction.
For a one-dimensional direction, the Born-von Karman boundary condition gives

.. math::

   q_{\alpha}
   =
   \frac{2\pi}{N_{\alpha}a_{\alpha}} n_{\alpha},

where :math:`a_{\alpha}` is the lattice constant, :math:`N_{\alpha}` is the
number of unit cells, and :math:`n_{\alpha}` is an integer.

For general three-dimensional cells, pySED writes the supercell as

.. math::

   \mathbf{S} = \mathbf{P}\mathbf{p},

where :math:`\mathbf{p}` contains the primitive lattice vectors as rows,
:math:`\mathbf{S}` contains the supercell lattice vectors as rows, and
:math:`\mathbf{P}` is an integer repetition matrix.

A reduced wave vector :math:`\mathbf{q}_{\mathrm{red}}` is commensurate with
the supercell if

.. math::

   \mathbf{q}_{\mathrm{red}}\mathbf{P}^{T} \in \mathbb{Z}^{3}.

For a path between two high-symmetry points, pySED searches the line

.. math::

   \mathbf{q}(f)
   =
   \mathbf{q}_{\mathrm{start}}
   +
   f\left(\mathbf{q}_{\mathrm{end}}-\mathbf{q}_{\mathrm{start}}\right),
   \qquad 0 \le f \le 1,

and keeps only values of :math:`f` satisfying

.. math::

   \mathbf{q}(f)\mathbf{P}^{T} \in \mathbb{Z}^{3}.

The resulting reduced q-points are mapped to Cartesian reciprocal space by

.. math::

   \mathbf{q}_{\mathrm{cart}}(f)
   =
   \mathbf{q}(f)\left[2\pi(\mathbf{p}^{-1})^{T}\right].

Larger supercells provide denser commensurate q-point sampling, but also
increase trajectory size, memory use, and SED compute time.

For a diagonal repetition matrix
:math:`\mathbf{P}=\mathrm{diag}(N_x,N_y,N_z)`, a Gamma-to-boundary path along
one reduced reciprocal axis contains :math:`\lfloor N_i/2 \rfloor + 1`
commensurate q-points, or :math:`N_i/2+1` when :math:`N_i` is even. For
non-diagonal repetition matrices, the count must be obtained from the integer
condition :math:`\mathbf{q}(f)\mathbf{P}^{T}\in\mathbb{Z}^{3}` for the selected
path segment.

Lorentzian Fitting and Lifetime
-------------------------------

The linewidth of an SED peak can be fitted with a Lorentzian function:

.. math::

   \Phi(\mathbf{q},\omega), \Phi'(\mathbf{q},\omega)
   =
   \frac{I}{1+\left[(\omega-\omega_c)/\gamma\right]^2}.

Here :math:`I` is the peak magnitude, :math:`\omega_c` is the center frequency,
and :math:`\gamma` is the half-width at half-maximum (HWHM). The paper defines
the lifetime as

.. math::

   \tau(\mathbf{q},\omega)
   =
   \frac{1}{2\gamma}.

In the current code output, frequencies are handled in THz for plotting and
fitting, and the lifetime is written in ps using the relation implemented in
``pySED/FileIO.py``:

.. math::

   \tau_{\mathrm{ps}}
   =
   \frac{1}{2\pi\gamma_{\mathrm{THz}}}.

Use the single-q-point fitting workflow first to tune ``peak_height``,
``peak_prominence``, ``initial_guess_hwhm``, ``peak_max_hwhm``, and
``lorentz_fit_cutoff`` before fitting all q-points.

LO-TO Splitting
---------------

LO-TO splitting in polar materials comes from long-range Coulomb interactions
that create a macroscopic electric field near the Gamma point. pySED is a
post-processing tool: it does not add a non-analytical correction (NAC) to the
trajectory after the MD simulation.

Therefore:

- If the MD trajectory is generated with a model that does not include the
  long-range electrostatic physics, pySED will not create LO-TO splitting by
  itself.
- If the trajectory is generated with a model that includes the relevant
  long-range Coulomb effects, pySED can resolve the splitting present in that
  trajectory.
- When comparing against lattice dynamics, apply NAC in the lattice-dynamics
  reference calculation when appropriate, for example through phonopy with Born
  effective charges and dielectric constants.

Small differences between finite-temperature SED and lattice-dynamics reference
branches can come from anharmonicity, the exchange-correlation level used to
generate Born charges, or differences between the MD potential and the reference
lattice-dynamics model. For qNEP-style long-range electrostatic models and
their use in polar materials, see the qNEP article [qNEP2026]_.

References
----------

.. [Liang2025] T. Liang, W. Jiang, K. Xu, H. Bu, Z. Fan, W. Ouyang, and
   J. Xu, "PYSED: A tool for extracting kinetic-energy-weighted phonon
   dispersion and lifetime from molecular dynamics simulations," *Journal of
   Applied Physics* **138**, 075101 (2025).
   https://doi.org/10.1063/5.0278798

.. [Thomas2010] J. A. Thomas, J. E. Turney, R. M. Iutzi, C. H. Amon, and
   A. J. H. McGaughey, "Predicting phonon dispersion relations and lifetimes
   from the spectral energy density," *Physical Review B* **81**, 081411
   (2010). https://doi.org/10.1103/PhysRevB.81.081411

.. [Dynasor2025] E. Berger, E. Fransson, F. Eriksson, E. Lindgren,
   G. Wahnström, T. H. Rod, and P. Erhart, "Dynasor 2: From simulation to
   experiment through correlation functions," *Computer Physics
   Communications* **316**, 109759 (2025).
   https://doi.org/10.1016/j.cpc.2025.109759

.. [qNEP2026] Z. Fan, B. Tang, E. Berger, E. Berger, E. Fransson, K. Xu,
   Z. Yan, Z. Liu, Z. Song, H. Dong, S. Chen, L. Li, Z. Wang, Y. Zhu,
   J. Wiktor, and P. Erhart, "qNEP: A highly efficient neuroevolution
   potential with dynamic charges for large-scale atomistic simulations,"
   *Journal of Chemical Theory and Computation* **22**, 4787-4801 (2026).
   https://doi.org/10.1021/acs.jctc.6c00146
