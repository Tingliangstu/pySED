"""Select one execution backend; no plugin registry or inheritance tree.

Both backends are context managers with four operations:
prepare(vels, cell_vecs, result, accumulator, prefetch=...) reports memory;
compute(vels, cell_vecs, result) fills unnormalized NumPy SED in place;
report(total_elapsed, trajectory_wait) prints optional backend diagnostics;
close() releases owned resources, including on exceptions.
"""


def create_backend(name, config, reduced_qpoints, max_cores):
    """Load GPU support only when explicitly selected."""
    if name == "numpy":
        from pySED.backends.cpu import CPUBackend
        return CPUBackend(config, reduced_qpoints, max_cores)
    if name == "cupy":
        from pySED.backends.gpu import GPUBackend
        return GPUBackend(config, reduced_qpoints)
    raise ValueError(f"unsupported SED backend '{name}'")
