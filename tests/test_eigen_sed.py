import numpy as np

from pySED.eigen_sed import EigenvectorSet, compute_eigen_sed


def _single_mode_inputs():
    n_frames = 64
    dt = 1e-12
    target_bin = 5
    time = np.arange(n_frames) * dt
    freq_hz = target_bin / (n_frames * dt)

    velocities = np.zeros((n_frames, 1, 3))
    velocities[:, 0, 0] = np.cos(2 * np.pi * freq_hz * time)
    unitcell_vectors = np.zeros((1, 3))
    basis_index = np.array([1])
    masses = np.array([1.0])
    primitive = np.eye(3)
    supercell = np.eye(3)

    eig = np.zeros((1, 1, 1, 3), dtype=complex)
    eig[0, 0, 0, 0] = 1.0
    eigset = EigenvectorSet(
        qpoints_reduced=np.array([[0.0, 0.0, 0.0]]),
        frequencies=np.array([[freq_hz / 1e12]]),
        eigenvectors=eig,
    )

    return {
        "velocities": velocities,
        "unitcell_vectors": unitcell_vectors,
        "basis_index": basis_index,
        "masses": masses,
        "primitive_cell": primitive,
        "supercell_cell": supercell,
        "eigenvectors": eigset,
        "dt": dt,
        "target_bin": target_bin,
    }


def test_eigen_sed_peak_frequency_for_single_mode():
    inputs = _single_mode_inputs()
    result = compute_eigen_sed(
        inputs["velocities"],
        inputs["unitcell_vectors"],
        inputs["basis_index"],
        inputs["masses"],
        inputs["primitive_cell"],
        inputs["supercell_cell"],
        inputs["eigenvectors"],
        inputs["dt"],
    )

    peak = int(np.argmax(result.sed[:, 0, 0]))
    assert peak == inputs["target_bin"]
    assert result.metadata["backend"] == "cpu"


def test_eigen_sed_cpu_backend_matches_default():
    inputs = _single_mode_inputs()

    default = compute_eigen_sed(
        inputs["velocities"],
        inputs["unitcell_vectors"],
        inputs["basis_index"],
        inputs["masses"],
        inputs["primitive_cell"],
        inputs["supercell_cell"],
        inputs["eigenvectors"],
        inputs["dt"],
    )
    cpu = compute_eigen_sed(
        inputs["velocities"],
        inputs["unitcell_vectors"],
        inputs["basis_index"],
        inputs["masses"],
        inputs["primitive_cell"],
        inputs["supercell_cell"],
        inputs["eigenvectors"],
        inputs["dt"],
        backend="cpu",
    )

    np.testing.assert_allclose(cpu.sed, default.sed)
    np.testing.assert_allclose(cpu.frequencies_thz, default.frequencies_thz)
    assert cpu.metadata["backend"] == "cpu"
