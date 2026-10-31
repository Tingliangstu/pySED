# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     pySED is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with pySED.  If not, see <http://www.gnu.org/licenses/>.
# =============================================================================

"""Select the CPU or GPU implementation used by the SED driver.

Both backends are context managers: entering prepares reusable resources,
and leaving releases them even when calculation fails. The driver calls
prepare() once to check memory, compute() for each trajectory block, and
report() after all blocks. compute() fills a NumPy result array in place;
normalization and frequency folding are handled by Phonon and sed_kernel.
"""


def create_backend(device, config, reduced_qpoints, max_cores):
    """Return the implementation for cpu or gpu without starting its resources.

    Use the result in a with block. NumPy supports serial or process-based
    execution; CuPy is imported only when the GPU backend is selected.
    """
    if device == "cpu":
        from pySED.backends.cpu import CPUBackend
        return CPUBackend(config, reduced_qpoints, max_cores)
    if device == "gpu":
        from pySED.backends.gpu import GPUBackend
        return GPUBackend(config, reduced_qpoints)
    raise ValueError(f"unsupported SED device '{device}'")
