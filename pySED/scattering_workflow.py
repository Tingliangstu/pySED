"""High-level workflows for experiment-facing scattering maps."""

from dataclasses import dataclass

import numpy as np

from pySED.compare_experiment import ScatteringMap, prepare_map_comparison
from pySED.dsf import compute_current_correlations, compute_dsf, compute_partial_dsf
from pySED.electron_kinematics import eels_qpoints_from_angle_axis
from pySED.eels import build_eels_map_from_mode_spectra, compute_mode_visibility_decomposition
from pySED.eigen_sed import compute_eigen_sed
from pySED.mode_visibility import build_one_phonon_map_from_mode_spectra
from pySED.mode_visibility import compute_one_phonon_visibility
from pySED.quantum_correction import convert_energy_axis
from pySED.q_advisor import QAdvisor


@dataclass
class EELSWorkflowResult:
    """Outputs from the MD-to-extended-zone-EELS workflow."""

    advisor: QAdvisor
    q_advice: list
    eigen_sed: object
    visibility: object
    eels_map: object
    scattering_map: ScatteringMap
    comparison: object = None
    q_plan: object = None


@dataclass
class DSFWorkflowResult:
    """Outputs from q-advised DSF calculation."""

    advisor: QAdvisor
    q_advice: list
    qpoints_reduced: np.ndarray
    qpoints_cartesian: np.ndarray
    dsf: object
    scattering_map: ScatteringMap = None
    comparison: object = None


@dataclass
class PartialDSFWorkflowResult:
    """Outputs from q-advised species-partial DSF calculation."""

    advisor: QAdvisor
    q_advice: list
    qpoints_reduced: np.ndarray
    qpoints_cartesian: np.ndarray
    partial_dsf: object
    scattering_map: ScatteringMap = None
    comparison: object = None


@dataclass
class CurrentCorrelationWorkflowResult:
    """Outputs from q-advised longitudinal/transverse current spectra."""

    advisor: QAdvisor
    q_advice: list
    qpoints_reduced: np.ndarray
    qpoints_cartesian: np.ndarray
    current: object
    scattering_map: ScatteringMap = None
    comparison: object = None


@dataclass
class OnePhononWorkflowResult:
    """Outputs from eigen-SED-to-INS/IXS one-phonon workflow."""

    advisor: QAdvisor
    q_advice: list
    eigen_sed: object
    visibility: object
    one_phonon_map: object
    scattering_map: ScatteringMap
    comparison: object = None
    q_plan: object = None


def q_path_axis(qpoints_cartesian):
    """Return cumulative path distance for ordered Cartesian Q points."""

    qpoints = np.atleast_2d(np.asarray(qpoints_cartesian, dtype=float))
    if qpoints.shape[1] != 3:
        raise ValueError("qpoints_cartesian must have shape (num_qpoints, 3)")
    if qpoints.shape[0] == 1:
        return np.zeros(1, dtype=float)
    steps = np.linalg.norm(np.diff(qpoints, axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(steps)])


def _requested_qpoints_to_advice(advisor, qpoints, qpoint_coordinates):
    coordinate_key = str(qpoint_coordinates).lower()
    if coordinate_key in ("reduced", "fractional", "crystal"):
        q_reduced = np.atleast_2d(np.asarray(qpoints, dtype=float))
        if q_reduced.shape[1] != 3:
            raise ValueError("reduced qpoints must have shape (num_qpoints, 3)")
        advice = advisor.advise_qpoints(q_reduced, preserve_image=True)
        return q_reduced, advisor.reduced_to_cartesian(q_reduced), advice

    if coordinate_key in ("cartesian", "cart", "angstrom^-1", "1/angstrom"):
        q_cart = np.atleast_2d(np.asarray(qpoints, dtype=float))
        if q_cart.shape[1] != 3:
            raise ValueError("Cartesian qpoints must have shape (num_qpoints, 3)")
        q_reduced = advisor.cartesian_to_reduced(q_cart)
        advice = advisor.advise_cartesian_qpoints(q_cart, preserve_image=True)
        return q_reduced, q_cart, advice

    raise ValueError("qpoint_coordinates must be 'reduced' or 'cartesian'")


def _enforce_q_policy(q_advice, q_policy, max_error_reduced=None, max_error_cartesian=None):
    policy = str(q_policy).lower()
    if policy not in ("strict", "nearest"):
        raise ValueError("q_policy must be 'strict' or 'nearest'")

    non_commensurate = [advice for advice in q_advice if not advice.is_commensurate]
    if policy == "strict" and non_commensurate:
        worst = max(non_commensurate, key=lambda item: item.error_cartesian)
        raise ValueError(
            "requested q points are not commensurate with the MD supercell; "
            "largest Cartesian error to nearest allowed Q is %.6g" % worst.error_cartesian
        )

    if policy == "nearest":
        if max_error_reduced is not None:
            too_large = [advice for advice in q_advice if advice.error_reduced > max_error_reduced]
            if too_large:
                raise ValueError("nearest commensurate q exceeds max_error_reduced")
        if max_error_cartesian is not None:
            too_large = [advice for advice in q_advice if advice.error_cartesian > max_error_cartesian]
            if too_large:
                raise ValueError("nearest commensurate q exceeds max_error_cartesian")

    return policy


def compute_dsf_workflow(
    positions,
    qpoints,
    primitive_cell,
    supercell_cell,
    dt,
    qpoint_coordinates="reduced",
    q_policy="strict",
    max_error_reduced=None,
    max_error_cartesian=None,
    map_component="total",
    q_axis=None,
    experimental_map=None,
    sigma_q=None,
    sigma_energy=None,
    normalization="max",
    quantum_temperature=None,
    output_energy_unit="thz",
    **dsf_kwargs,
):
    """Validate or map requested Q points before computing DSF.

    ``q_policy="strict"`` rejects non-commensurate points.  ``"nearest"`` maps
    each requested point to the nearest commensurate point while preserving the
    extended-zone reciprocal-lattice image, so ``Q = q + G`` is not folded back
    to the first Brillouin zone.
    """

    advisor = QAdvisor(primitive_cell, supercell_cell)
    q_reduced, q_cartesian, q_advice = _requested_qpoints_to_advice(
        advisor,
        qpoints,
        qpoint_coordinates,
    )
    policy = _enforce_q_policy(
        q_advice,
        q_policy,
        max_error_reduced=max_error_reduced,
        max_error_cartesian=max_error_cartesian,
    )

    if policy == "nearest":
        q_reduced = np.asarray([advice.nearest_reduced for advice in q_advice], dtype=float)
        q_cartesian = np.asarray([advice.nearest_cartesian for advice in q_advice], dtype=float)

    result = compute_dsf(positions, q_cartesian, dt, **dsf_kwargs)
    result.metadata.update(
        {
            "q_policy": policy,
            "qpoint_coordinates": qpoint_coordinates,
            "q_advisor": "finite-supercell commensurability enforced",
            "max_q_error_reduced": max((advice.error_reduced for advice in q_advice), default=0.0),
            "max_q_error_cartesian": max((advice.error_cartesian for advice in q_advice), default=0.0),
        }
    )

    scattering_map = dsf_to_scattering_map(
        result,
        component=map_component,
        q_axis=q_axis,
        source_energy_unit="thz",
        output_energy_unit=output_energy_unit,
    )
    comparison = None
    if experimental_map is not None:
        comparison = prepare_map_comparison(
            scattering_map,
            experimental_map,
            sigma_q=sigma_q,
            sigma_energy=sigma_energy,
            normalization=normalization,
            quantum_temperature=quantum_temperature,
            energy_unit=output_energy_unit,
        )

    return DSFWorkflowResult(
        advisor=advisor,
        q_advice=q_advice,
        qpoints_reduced=q_reduced,
        qpoints_cartesian=q_cartesian,
        dsf=result,
        scattering_map=scattering_map,
        comparison=comparison,
    )


def compute_partial_dsf_workflow(
    positions,
    qpoints,
    primitive_cell,
    supercell_cell,
    dt,
    atom_types,
    qpoint_coordinates="reduced",
    q_policy="strict",
    max_error_reduced=None,
    max_error_cartesian=None,
    map_component="weighted_total",
    species_pair=None,
    complex_part="real",
    q_axis=None,
    experimental_map=None,
    sigma_q=None,
    sigma_energy=None,
    normalization="max",
    quantum_temperature=None,
    output_energy_unit="thz",
    **partial_kwargs,
):
    """Validate/map requested Q points before computing species-partial DSF."""

    advisor = QAdvisor(primitive_cell, supercell_cell)
    q_reduced, q_cartesian, q_advice = _requested_qpoints_to_advice(
        advisor,
        qpoints,
        qpoint_coordinates,
    )
    policy = _enforce_q_policy(
        q_advice,
        q_policy,
        max_error_reduced=max_error_reduced,
        max_error_cartesian=max_error_cartesian,
    )

    if policy == "nearest":
        q_reduced = np.asarray([advice.nearest_reduced for advice in q_advice], dtype=float)
        q_cartesian = np.asarray([advice.nearest_cartesian for advice in q_advice], dtype=float)

    result = compute_partial_dsf(
        positions,
        q_cartesian,
        dt,
        atom_types=atom_types,
        **partial_kwargs,
    )
    result.metadata.update(
        {
            "q_policy": policy,
            "qpoint_coordinates": qpoint_coordinates,
            "q_advisor": "finite-supercell commensurability enforced",
            "max_q_error_reduced": max((advice.error_reduced for advice in q_advice), default=0.0),
            "max_q_error_cartesian": max((advice.error_cartesian for advice in q_advice), default=0.0),
        }
    )

    scattering_map = partial_dsf_to_scattering_map(
        result,
        component=map_component,
        species_pair=species_pair,
        complex_part=complex_part,
        q_axis=q_axis,
        source_energy_unit="thz",
        output_energy_unit=output_energy_unit,
    )
    comparison = None
    if experimental_map is not None:
        comparison = prepare_map_comparison(
            scattering_map,
            experimental_map,
            sigma_q=sigma_q,
            sigma_energy=sigma_energy,
            normalization=normalization,
            quantum_temperature=quantum_temperature,
            energy_unit=output_energy_unit,
        )

    return PartialDSFWorkflowResult(
        advisor=advisor,
        q_advice=q_advice,
        qpoints_reduced=q_reduced,
        qpoints_cartesian=q_cartesian,
        partial_dsf=result,
        scattering_map=scattering_map,
        comparison=comparison,
    )


def compute_current_correlation_workflow(
    positions,
    velocities,
    qpoints,
    primitive_cell,
    supercell_cell,
    dt,
    qpoint_coordinates="reduced",
    q_policy="strict",
    max_error_reduced=None,
    max_error_cartesian=None,
    map_component="longitudinal",
    q_axis=None,
    experimental_map=None,
    sigma_q=None,
    sigma_energy=None,
    normalization="max",
    quantum_temperature=None,
    output_energy_unit="thz",
    **current_kwargs,
):
    """Validate/map requested Q points before computing current spectra."""

    advisor = QAdvisor(primitive_cell, supercell_cell)
    q_reduced, q_cartesian, q_advice = _requested_qpoints_to_advice(
        advisor,
        qpoints,
        qpoint_coordinates,
    )
    policy = _enforce_q_policy(
        q_advice,
        q_policy,
        max_error_reduced=max_error_reduced,
        max_error_cartesian=max_error_cartesian,
    )

    if policy == "nearest":
        q_reduced = np.asarray([advice.nearest_reduced for advice in q_advice], dtype=float)
        q_cartesian = np.asarray([advice.nearest_cartesian for advice in q_advice], dtype=float)

    result = compute_current_correlations(
        positions,
        velocities,
        q_cartesian,
        dt,
        **current_kwargs,
    )
    result.metadata.update(
        {
            "q_policy": policy,
            "qpoint_coordinates": qpoint_coordinates,
            "q_advisor": "finite-supercell commensurability enforced",
            "max_q_error_reduced": max((advice.error_reduced for advice in q_advice), default=0.0),
            "max_q_error_cartesian": max((advice.error_cartesian for advice in q_advice), default=0.0),
        }
    )

    scattering_map = current_correlation_to_scattering_map(
        result,
        component=map_component,
        q_axis=q_axis,
        source_energy_unit="thz",
        output_energy_unit=output_energy_unit,
    )
    comparison = None
    if experimental_map is not None:
        comparison = prepare_map_comparison(
            scattering_map,
            experimental_map,
            sigma_q=sigma_q,
            sigma_energy=sigma_energy,
            normalization=normalization,
            quantum_temperature=quantum_temperature,
            energy_unit=output_energy_unit,
        )

    return CurrentCorrelationWorkflowResult(
        advisor=advisor,
        q_advice=q_advice,
        qpoints_reduced=q_reduced,
        qpoints_cartesian=q_cartesian,
        current=result,
        scattering_map=scattering_map,
        comparison=comparison,
    )


def dsf_to_scattering_map(
    dsf_result,
    component="total",
    q_axis=None,
    source_energy_unit="thz",
    output_energy_unit="thz",
):
    """Convert a DSF result component into a q-energy scattering map."""

    component_key = str(component).lower()
    if component_key not in ("total", "coherent", "incoherent"):
        raise ValueError("component must be 'total', 'coherent', or 'incoherent'")

    intensity = np.asarray(getattr(dsf_result, component_key), dtype=float)
    if intensity.ndim != 2:
        raise ValueError("DSF component must have shape (num_q, num_energy)")

    q_axis_source = "user supplied"
    if q_axis is None:
        q_axis = q_path_axis(dsf_result.qpoints_cartesian)
        q_axis_source = "Cartesian path distance"
    q_axis = np.asarray(q_axis, dtype=float)
    if q_axis.shape != (intensity.shape[0],):
        raise ValueError("q_axis length must match number of q points")

    energy_axis = np.asarray(dsf_result.frequencies_thz, dtype=float)
    if output_energy_unit != source_energy_unit:
        energy_axis = convert_energy_axis(
            energy_axis,
            source_unit=source_energy_unit,
            target_unit=output_energy_unit,
        )

    metadata = {
        "source": "pySED DSFResult",
        "component": component_key,
        "q_axis": q_axis_source,
        "energy_unit": output_energy_unit,
        "experiment": dsf_result.metadata.get("experiment") if dsf_result.metadata else None,
    }
    return ScatteringMap(q_axis, energy_axis, intensity, metadata)


def _species_pair_indices(species, species_pair):
    if species_pair is None:
        raise ValueError("species_pair must be supplied when component='partial'")
    if len(species_pair) != 2:
        raise ValueError("species_pair must contain two labels or indices")

    indices = []
    for item in species_pair:
        if isinstance(item, (int, np.integer)):
            index = int(item)
            if not 0 <= index < len(species):
                raise ValueError("species_pair index is outside the available species")
            indices.append(index)
        else:
            label = str(item)
            if label not in species:
                raise ValueError("species_pair label %s is not present in partial DSF" % label)
            indices.append(species.index(label))
    return tuple(indices)


def partial_dsf_to_scattering_map(
    partial_result,
    component="weighted_total",
    species_pair=None,
    complex_part="real",
    q_axis=None,
    source_energy_unit="thz",
    output_energy_unit="thz",
):
    """Convert a partial DSF component into a q-energy scattering map.

    ``component="weighted_total"`` returns the probe-weighted coherent total.
    ``component="partial"`` selects one species pair from the complex partial
    matrix and returns its real part, imaginary part, or magnitude.
    """

    component_key = str(component).lower()
    if component_key in ("weighted_total", "total", "coherent"):
        intensity = np.asarray(partial_result.weighted_total, dtype=float)
        selected_pair = None
        selected_part = None
    elif component_key in ("partial", "species_pair", "pair"):
        pair_indices = _species_pair_indices(partial_result.species, species_pair)
        selected = np.asarray(partial_result.partial[:, pair_indices[0], pair_indices[1], :])
        part_key = str(complex_part).lower()
        if part_key == "real":
            intensity = selected.real
        elif part_key in ("imag", "imaginary"):
            intensity = selected.imag
        elif part_key in ("abs", "magnitude"):
            intensity = np.abs(selected)
        else:
            raise ValueError("complex_part must be 'real', 'imag', or 'abs'")
        selected_pair = (
            partial_result.species[pair_indices[0]],
            partial_result.species[pair_indices[1]],
        )
        selected_part = part_key
    else:
        raise ValueError("component must be 'weighted_total' or 'partial'")

    if intensity.ndim != 2:
        raise ValueError("partial DSF map component must have shape (num_q, num_energy)")

    q_axis_source = "user supplied"
    if q_axis is None:
        q_axis = q_path_axis(partial_result.qpoints_cartesian)
        q_axis_source = "Cartesian path distance"
    q_axis = np.asarray(q_axis, dtype=float)
    if q_axis.shape != (intensity.shape[0],):
        raise ValueError("q_axis length must match number of q points")

    energy_axis = np.asarray(partial_result.frequencies_thz, dtype=float)
    if output_energy_unit != source_energy_unit:
        energy_axis = convert_energy_axis(
            energy_axis,
            source_unit=source_energy_unit,
            target_unit=output_energy_unit,
        )

    metadata = {
        "source": "pySED PartialDSFResult",
        "component": component_key,
        "species_pair": selected_pair,
        "complex_part": selected_part,
        "q_axis": q_axis_source,
        "energy_unit": output_energy_unit,
    }
    return ScatteringMap(q_axis, energy_axis, intensity, metadata)


def current_correlation_to_scattering_map(
    current_result,
    component="longitudinal",
    q_axis=None,
    source_energy_unit="thz",
    output_energy_unit="thz",
):
    """Convert a longitudinal/transverse current spectrum to a q-energy map."""

    component_key = str(component).lower()
    if component_key not in ("longitudinal", "transverse"):
        raise ValueError("component must be 'longitudinal' or 'transverse'")

    intensity = np.asarray(getattr(current_result, component_key), dtype=float)
    if intensity.ndim != 2:
        raise ValueError("current spectrum must have shape (num_q, num_energy)")

    q_axis_source = "user supplied"
    if q_axis is None:
        q_axis = q_path_axis(current_result.qpoints_cartesian)
        q_axis_source = "Cartesian path distance"
    q_axis = np.asarray(q_axis, dtype=float)
    if q_axis.shape != (intensity.shape[0],):
        raise ValueError("q_axis length must match number of q points")

    energy_axis = np.asarray(current_result.frequencies_thz, dtype=float)
    if output_energy_unit != source_energy_unit:
        energy_axis = convert_energy_axis(
            energy_axis,
            source_unit=source_energy_unit,
            target_unit=output_energy_unit,
        )

    metadata = {
        "source": "pySED CurrentCorrelationResult",
        "component": component_key,
        "q_axis": q_axis_source,
        "energy_unit": output_energy_unit,
    }
    return ScatteringMap(q_axis, energy_axis, intensity, metadata)


def eels_map_to_scattering_map(
    eels_map,
    g_index=0,
    q_axis=None,
    source_energy_unit="thz",
    output_energy_unit="thz",
):
    """Select one extended-zone G branch as a 2D q-energy map.

    The returned q axis is the cumulative Cartesian path distance unless a
    user-supplied ``q_axis`` is provided.
    """

    intensity = np.asarray(eels_map.intensity, dtype=float)
    if intensity.ndim != 3:
        raise ValueError("eels_map intensity must have shape (num_q, num_g, num_energy)")
    if not 0 <= int(g_index) < intensity.shape[1]:
        raise ValueError("g_index is outside the available extended-zone branches")

    q_cart = np.asarray(eels_map.Q_cartesian[:, int(g_index), :], dtype=float)
    q_axis_source = "user supplied"
    if q_axis is None:
        q_axis = q_path_axis(q_cart)
        q_axis_source = "Cartesian path distance"
    q_axis = np.asarray(q_axis, dtype=float)
    if q_axis.shape != (intensity.shape[0],):
        raise ValueError("q_axis length must match number of q points")

    energy_axis = np.asarray(eels_map.energy_axis, dtype=float)
    if output_energy_unit != source_energy_unit:
        energy_axis = convert_energy_axis(
            energy_axis,
            source_unit=source_energy_unit,
            target_unit=output_energy_unit,
        )

    metadata = {
        "source": "pySED EELSMap",
        "g_index": int(g_index),
        "q_axis": q_axis_source,
        "energy_unit": output_energy_unit,
    }
    return ScatteringMap(q_axis, energy_axis, intensity[:, int(g_index), :], metadata)


def one_phonon_map_to_scattering_map(
    one_phonon_map,
    g_index=0,
    q_axis=None,
    source_energy_unit="thz",
    output_energy_unit="thz",
):
    """Select one extended-zone G branch from an INS/IXS one-phonon map."""

    intensity = np.asarray(one_phonon_map.intensity, dtype=float)
    if intensity.ndim != 3:
        raise ValueError("one_phonon_map intensity must have shape (num_q, num_g, num_energy)")
    if not 0 <= int(g_index) < intensity.shape[1]:
        raise ValueError("g_index is outside the available extended-zone branches")

    q_cart = np.asarray(one_phonon_map.Q_cartesian[:, int(g_index), :], dtype=float)
    q_axis_source = "user supplied"
    if q_axis is None:
        q_axis = q_path_axis(q_cart)
        q_axis_source = "Cartesian path distance"
    q_axis = np.asarray(q_axis, dtype=float)
    if q_axis.shape != (intensity.shape[0],):
        raise ValueError("q_axis length must match number of q points")

    energy_axis = np.asarray(one_phonon_map.energy_axis, dtype=float)
    if output_energy_unit != source_energy_unit:
        energy_axis = convert_energy_axis(
            energy_axis,
            source_unit=source_energy_unit,
            target_unit=output_energy_unit,
        )

    metadata = {
        "source": "pySED OnePhononMap",
        "g_index": int(g_index),
        "q_axis": q_axis_source,
        "energy_unit": output_energy_unit,
    }
    return ScatteringMap(q_axis, energy_axis, intensity[:, int(g_index), :], metadata)


def compute_one_phonon_workflow(
    velocities,
    unitcell_vectors,
    basis_index,
    masses,
    primitive_cell,
    supercell_cell,
    eigenvectors,
    basis_positions_reduced,
    g_vectors_reduced,
    dt,
    atom_types=None,
    experiment="neutron",
    scattering_weights=None,
    frequency_power=0,
    num_blocks=1,
    sed_backend="cpu",
    target_qpoints_reduced=None,
    pairwise_g_vectors=False,
    g_index=0,
    q_axis=None,
    experimental_map=None,
    sigma_q=None,
    sigma_energy=None,
    normalization="max",
    quantum_temperature=None,
    output_energy_unit="thz",
):
    """Run commensurate-q, eigen-SED, and one-phonon INS/IXS workflow."""

    advisor = QAdvisor(primitive_cell, supercell_cell)
    qpoints = np.asarray(eigenvectors.qpoints_reduced, dtype=float)
    q_advice = advisor.advise_qpoints(qpoints)
    if any(not advice.is_commensurate for advice in q_advice):
        raise ValueError("eigenvector q points must be commensurate with the MD supercell")

    eigen_sed = compute_eigen_sed(
        velocities,
        unitcell_vectors,
        basis_index,
        masses,
        primitive_cell,
        supercell_cell,
        eigenvectors,
        dt,
        num_blocks=num_blocks,
        target_qpoints=target_qpoints_reduced,
        validate_q=True,
        backend=sed_backend,
    )
    visibility = compute_one_phonon_visibility(
        eigen_sed.qpoints_reduced,
        g_vectors_reduced,
        primitive_cell,
        basis_positions_reduced,
        masses,
        eigenvectors.eigenvectors,
        atom_types=atom_types,
        experiment=experiment,
        scattering_weights=scattering_weights,
        mode_frequencies=eigenvectors.frequencies,
        frequency_power=frequency_power,
        pairwise_g_vectors=pairwise_g_vectors,
    )
    one_phonon_map = build_one_phonon_map_from_mode_spectra(
        visibility,
        eigen_sed.sed,
        eigen_sed.frequencies_thz,
    )
    scattering_map = one_phonon_map_to_scattering_map(
        one_phonon_map,
        g_index=g_index,
        q_axis=q_axis,
        source_energy_unit="thz",
        output_energy_unit=output_energy_unit,
    )

    comparison = None
    if experimental_map is not None:
        comparison = prepare_map_comparison(
            scattering_map,
            experimental_map,
            sigma_q=sigma_q,
            sigma_energy=sigma_energy,
            normalization=normalization,
            quantum_temperature=quantum_temperature,
            energy_unit=output_energy_unit,
        )

    return OnePhononWorkflowResult(
        advisor=advisor,
        q_advice=q_advice,
        eigen_sed=eigen_sed,
        visibility=visibility,
        one_phonon_map=one_phonon_map,
        scattering_map=scattering_map,
        comparison=comparison,
    )


def compute_one_phonon_workflow_from_q_path(
    velocities,
    unitcell_vectors,
    basis_index,
    masses,
    primitive_cell,
    supercell_cell,
    eigenvectors,
    basis_positions_reduced,
    experimental_qpoints,
    dt,
    qpoint_coordinates="cartesian",
    q_policy="strict",
    max_error_reduced=None,
    max_error_cartesian=None,
    atom_types=None,
    experiment="neutron",
    scattering_weights=None,
    frequency_power=0,
    num_blocks=1,
    sed_backend="cpu",
    q_axis=None,
    experimental_map=None,
    sigma_q=None,
    sigma_energy=None,
    normalization="max",
    quantum_temperature=None,
    output_energy_unit="thz",
):
    """Run coherent INS/IXS one-phonon workflow from an experimental Q path."""

    advisor = QAdvisor(primitive_cell, supercell_cell)
    q_plan = advisor.plan_extended_zone_q_path(
        experimental_qpoints,
        coordinates=qpoint_coordinates,
        q_policy=q_policy,
        max_error_reduced=max_error_reduced,
        max_error_cartesian=max_error_cartesian,
    )
    result = compute_one_phonon_workflow(
        velocities,
        unitcell_vectors,
        basis_index,
        masses,
        primitive_cell,
        supercell_cell,
        eigenvectors,
        basis_positions_reduced,
        q_plan.g_vectors_reduced,
        dt,
        atom_types=atom_types,
        experiment=experiment,
        scattering_weights=scattering_weights,
        frequency_power=frequency_power,
        num_blocks=num_blocks,
        sed_backend=sed_backend,
        target_qpoints_reduced=q_plan.qpoints_reduced,
        pairwise_g_vectors=True,
        g_index=0,
        q_axis=q_axis,
        experimental_map=experimental_map,
        sigma_q=sigma_q,
        sigma_energy=sigma_energy,
        normalization=normalization,
        quantum_temperature=quantum_temperature,
        output_energy_unit=output_energy_unit,
    )
    result.q_plan = q_plan
    result.eigen_sed.metadata.update(
        {
            "q_plan": "experimental extended-zone Q path folded to commensurate q",
            "q_policy": q_policy,
            "max_q_error_reduced": q_plan.report.max_error_reduced,
            "max_q_error_cartesian": q_plan.report.max_error_cartesian,
        }
    )
    return result


def compute_eels_workflow(
    velocities,
    unitcell_vectors,
    basis_index,
    masses,
    primitive_cell,
    supercell_cell,
    eigenvectors,
    basis_positions_reduced,
    g_vectors_reduced,
    dt,
    num_blocks=1,
    sed_backend="cpu",
    target_qpoints_reduced=None,
    electron_form_factors=None,
    atom_types=None,
    electron_form_factor_model="unit",
    electron_form_factor_table=None,
    electron_missing="raise",
    pairwise_g_vectors=False,
    g_index=0,
    q_axis=None,
    experimental_map=None,
    sigma_q=None,
    sigma_energy=None,
    normalization="max",
    quantum_temperature=None,
    output_energy_unit="thz",
):
    """Run the commensurate-q, eigen-SED, and kinematic EELS workflow.

    This function is intentionally a workflow wrapper around the lower-level
    modules.  It keeps the paper-facing sequence explicit:

    ``q_advisor -> eigen_sed -> EELS visibility -> extended-zone map``.
    """

    advisor = QAdvisor(primitive_cell, supercell_cell)
    qpoints = np.asarray(eigenvectors.qpoints_reduced, dtype=float)
    q_advice = advisor.advise_qpoints(qpoints)
    if any(not advice.is_commensurate for advice in q_advice):
        raise ValueError("eigenvector q points must be commensurate with the MD supercell")

    eigen_sed = compute_eigen_sed(
        velocities,
        unitcell_vectors,
        basis_index,
        masses,
        primitive_cell,
        supercell_cell,
        eigenvectors,
        dt,
        num_blocks=num_blocks,
        target_qpoints=target_qpoints_reduced,
        validate_q=True,
        backend=sed_backend,
    )

    visibility = compute_mode_visibility_decomposition(
        eigen_sed.qpoints_reduced,
        g_vectors_reduced,
        primitive_cell,
        basis_positions_reduced,
        masses,
        eigenvectors.eigenvectors,
        electron_form_factors=electron_form_factors,
        atom_types=atom_types,
        electron_form_factor_model=electron_form_factor_model,
        electron_form_factor_table=electron_form_factor_table,
        electron_missing=electron_missing,
        pairwise_g_vectors=pairwise_g_vectors,
    )
    eels_map = build_eels_map_from_mode_spectra(
        visibility,
        eigen_sed.sed,
        eigen_sed.frequencies_thz,
    )
    scattering_map = eels_map_to_scattering_map(
        eels_map,
        g_index=g_index,
        q_axis=q_axis,
        source_energy_unit="thz",
        output_energy_unit=output_energy_unit,
    )

    comparison = None
    if experimental_map is not None:
        comparison = prepare_map_comparison(
            scattering_map,
            experimental_map,
            sigma_q=sigma_q,
            sigma_energy=sigma_energy,
            normalization=normalization,
            quantum_temperature=quantum_temperature,
            energy_unit=output_energy_unit,
        )

    return EELSWorkflowResult(
        advisor=advisor,
        q_advice=q_advice,
        eigen_sed=eigen_sed,
        visibility=visibility,
        eels_map=eels_map,
        scattering_map=scattering_map,
        comparison=comparison,
    )


def compute_eels_workflow_from_q_path(
    velocities,
    unitcell_vectors,
    basis_index,
    masses,
    primitive_cell,
    supercell_cell,
    eigenvectors,
    basis_positions_reduced,
    experimental_qpoints,
    dt,
    qpoint_coordinates="cartesian",
    q_policy="strict",
    max_error_reduced=None,
    max_error_cartesian=None,
    num_blocks=1,
    sed_backend="cpu",
    electron_form_factors=None,
    atom_types=None,
    electron_form_factor_model="unit",
    electron_form_factor_table=None,
    electron_missing="raise",
    q_axis=None,
    experimental_map=None,
    sigma_q=None,
    sigma_energy=None,
    normalization="max",
    quantum_temperature=None,
    output_energy_unit="thz",
):
    """Run the EELS workflow directly from an experimental extended-zone Q path.

    This wrapper performs the full experiment-facing sequence:

    ``experimental Q -> q_plan(Q=q+G) -> eigen-SED at folded q -> pairwise EELS map``.
    """

    advisor = QAdvisor(primitive_cell, supercell_cell)
    q_plan = advisor.plan_extended_zone_q_path(
        experimental_qpoints,
        coordinates=qpoint_coordinates,
        q_policy=q_policy,
        max_error_reduced=max_error_reduced,
        max_error_cartesian=max_error_cartesian,
    )

    result = compute_eels_workflow(
        velocities,
        unitcell_vectors,
        basis_index,
        masses,
        primitive_cell,
        supercell_cell,
        eigenvectors,
        basis_positions_reduced,
        q_plan.g_vectors_reduced,
        dt,
        num_blocks=num_blocks,
        sed_backend=sed_backend,
        target_qpoints_reduced=q_plan.qpoints_reduced,
        electron_form_factors=electron_form_factors,
        atom_types=atom_types,
        electron_form_factor_model=electron_form_factor_model,
        electron_form_factor_table=electron_form_factor_table,
        electron_missing=electron_missing,
        pairwise_g_vectors=True,
        g_index=0,
        q_axis=q_axis,
        experimental_map=experimental_map,
        sigma_q=sigma_q,
        sigma_energy=sigma_energy,
        normalization=normalization,
        quantum_temperature=quantum_temperature,
        output_energy_unit=output_energy_unit,
    )
    result.q_plan = q_plan
    result.eigen_sed.metadata.update(
        {
            "q_plan": "experimental extended-zone Q path folded to commensurate q",
            "q_policy": q_policy,
            "max_q_error_reduced": q_plan.report.max_error_reduced,
            "max_q_error_cartesian": q_plan.report.max_error_cartesian,
        }
    )
    return result


def compute_eels_workflow_from_angle_axis(
    velocities,
    unitcell_vectors,
    basis_index,
    masses,
    primitive_cell,
    supercell_cell,
    eigenvectors,
    basis_positions_reduced,
    angle_axis,
    beam_energy_ev,
    dt,
    direction=(1.0, 0.0),
    energy_loss_ev=0.0,
    angle_unit="mrad",
    include_longitudinal=False,
    q_policy="strict",
    max_error_reduced=None,
    max_error_cartesian=None,
    num_blocks=1,
    sed_backend="cpu",
    electron_form_factors=None,
    atom_types=None,
    electron_form_factor_model="unit",
    electron_form_factor_table=None,
    electron_missing="raise",
    q_axis=None,
    experimental_map=None,
    sigma_q=None,
    sigma_energy=None,
    normalization="max",
    quantum_temperature=None,
    output_energy_unit="thz",
):
    """Run q-EELS workflow from a calibrated detector scattering-angle axis.

    The angle axis is converted to Cartesian momentum transfer with
    :func:`pySED.electron_kinematics.eels_qpoints_from_angle_axis`, then the
    standard extended-zone q-path workflow is used.
    """

    experimental_qpoints = eels_qpoints_from_angle_axis(
        angle_axis,
        beam_energy_ev=beam_energy_ev,
        direction=direction,
        energy_loss_ev=energy_loss_ev,
        angle_unit=angle_unit,
        include_longitudinal=include_longitudinal,
    )
    result = compute_eels_workflow_from_q_path(
        velocities=velocities,
        unitcell_vectors=unitcell_vectors,
        basis_index=basis_index,
        masses=masses,
        primitive_cell=primitive_cell,
        supercell_cell=supercell_cell,
        eigenvectors=eigenvectors,
        basis_positions_reduced=basis_positions_reduced,
        experimental_qpoints=experimental_qpoints,
        dt=dt,
        qpoint_coordinates="cartesian",
        q_policy=q_policy,
        max_error_reduced=max_error_reduced,
        max_error_cartesian=max_error_cartesian,
        num_blocks=num_blocks,
        sed_backend=sed_backend,
        electron_form_factors=electron_form_factors,
        atom_types=atom_types,
        electron_form_factor_model=electron_form_factor_model,
        electron_form_factor_table=electron_form_factor_table,
        electron_missing=electron_missing,
        q_axis=q_axis,
        experimental_map=experimental_map,
        sigma_q=sigma_q,
        sigma_energy=sigma_energy,
        normalization=normalization,
        quantum_temperature=quantum_temperature,
        output_energy_unit=output_energy_unit,
    )
    result.eigen_sed.metadata.update(
        {
            "q_plan_source": "electron scattering angle axis",
            "beam_energy_ev": float(beam_energy_ev),
            "angle_unit": str(angle_unit),
            "include_longitudinal_q": bool(include_longitudinal),
        }
    )
    result.visibility.metadata.update(
        {
            "q_plan_source": "electron scattering angle axis",
            "beam_energy_ev": float(beam_energy_ev),
            "angle_unit": str(angle_unit),
            "include_longitudinal_q": bool(include_longitudinal),
        }
    )
    return result
