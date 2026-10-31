"""Select the CPU or GPU implementation used by the SED driver.

Both backends are context managers: entering prepares reusable resources,
and leaving releases them even when calculation fails. The driver calls
prepare() once to check memory, compute() for each trajectory block, and
report() after all blocks. compute() fills a NumPy result array in place;
normalization and frequency folding are handled by Phonon and sed_kernel.
"""


def create_backend(name, config, reduced_qpoints, max_cores):
    """Return the backend selected by name, without starting its resources.

    Use the result in a with block. NumPy supports serial or process-based
    execution; CuPy is imported only when the GPU backend is selected.
    """
    if name == "numpy":
        from pySED.backends.cpu import CPUBackend
        return CPUBackend(config, reduced_qpoints, max_cores)
    if name == "cupy":
        from pySED.backends.gpu import GPUBackend
        return GPUBackend(config, reduced_qpoints)
    raise ValueError(f"unsupported SED backend '{name}'")
