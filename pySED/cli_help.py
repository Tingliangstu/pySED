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

"""Offline, parameter-specific command-line help for the public input options."""

from dataclasses import dataclass
from difflib import get_close_matches
import textwrap

from pySED.parameters import PARAMETERS


MANUAL_URL = "https://pysed.readthedocs.io/"
_WIDTH = 96


@dataclass(frozen=True)
class _Topic:
    syntax: str
    meaning: str
    notes: tuple[str, ...] = ()
    values: tuple[tuple[str, str], ...] = ()
    examples: tuple[str, ...] = ()
    related: tuple[str, ...] = ()
    default_note: str = ""
    required_for_compute: bool = False


# These terminal summaries follow docs/source/input_parameters/. Keep detailed
# derivations in the manual; obtain actual default values from PARAMETERS.
_TOPICS = {
    "action": _Topic(
        "action = thinking",
        "Selects the work performed by this invocation of pySED.",
        values=(
            ("thinking", "Choose the next work based on existing output files."),
            ("compute", "Recalculate SED from the trajectory and overwrite the main data files."),
            ("plot", "Draw the requested figures from existing SED data."),
            ("fit", "Detect and fit spectral peaks in existing SED data."),
        ),
        notes=(
            "When SED data are missing, thinking computes SED and creates plots in the "
            "same invocation. With SED data present, it plots or fits if the corresponding "
            "outputs are missing.",
            "thinking checks whether output files exist; it does not compare the settings "
            "used to create them with the current input. Select compute, plot, or fit "
            "explicitly after changing the relevant settings.",
            "For fitting, lorentz_fit_all_qpoint selects all q-points or only the one "
            "specified by qpoint_slice_index. Plotting and fitting require saved SED data.",
        ),
        examples=("action = compute", "action = plot", "action = fit"),
        related=("trajectory_file", "qpoint_slice_index", "lorentz_fit_all_qpoint"),
    ),
    "device": _Topic(
        "device = cpu",
        "Selects the device used to calculate SED. Plotting, peak detection, and fitting "
        "remain CPU operations.",
        values=(
            ("cpu", "Use CPU computation; max_cores controls the number of worker processes."),
            ("gpu", "Use one CUDA GPU per pySED process; a compatible CuPy installation is required."),
        ),
        notes=(
            "max_cores = 1 runs the CPU calculation serially. Larger values enable "
            "persistent worker processes. max_cores does not control GPU computation.",
            "GPU operation still needs CPU memory for trajectory data. Device selection "
            "does not change the SED definition or input trajectory formats.",
        ),
        related=("max_cores", "num_blocks", "trajectory_prefetch"),
    ),
    "num_atoms": _Topic(
        "num_atoms = 13824",
        "Number of atoms represented by the trajectory. It must agree with both the "
        "trajectory and the atom mapping in basis.in.",
        notes=("For a trajectory containing one dumped group, use the number of dumped "
               "atoms and a basis.in describing that same group.",),
        related=("basis_lattice_file", "total_num_steps", "dump_group_supercell"),
        required_for_compute=True,
    ),
    "total_num_steps": _Topic(
        "total_num_steps = 500000",
        "Number of MD steps in the production segment to analyze. This is not the number "
        "of saved frames, and it does not include a separate equilibration run.",
        notes=(
            "The requested frame count is total_num_steps / output_data_stride. The "
            "trajectory must supply that many frames, and the frame count must divide "
            "into num_blocks complete, nonempty blocks.",
            "The analyzed duration is total_num_steps * time_step / 1000 ps. Each block "
            "has that duration divided by num_blocks. Longer blocks give finer FFT "
            "frequency spacing; increasing the number of blocks coarsens that spacing.",
        ),
        examples=(
            "total_num_steps = 500000\ntime_step = 1\noutput_data_stride = 50\nnum_blocks = 5",
        ),
        related=("time_step", "output_data_stride", "num_blocks", "initial_guess_hwhm"),
        required_for_compute=True,
    ),
    "time_step": _Topic(
        "time_step = 1",
        "MD integration time step in femtoseconds (fs). Use the value from the production "
        "MD run, not the time between saved frames.",
        notes=(
            "The saved-frame interval is time_step * output_data_stride fs. For example, "
            "time_step = 1 and output_data_stride = 50 correspond to one saved frame "
            "every 50 fs.",
        ),
        related=("output_data_stride", "total_num_steps", "num_blocks"),
        required_for_compute=True,
    ),
    "output_data_stride": _Topic(
        "output_data_stride = 50",
        "Number of MD steps between saved trajectory frames. This must match the output "
        "interval used to write positions and velocities during the production run.",
        notes=(
            "The Nyquist limit in THz is 500 / (time_step * output_data_stride), with "
            "time_step in fs. Modes above this limit cannot be recovered from the saved "
            "trajectory.",
            "For time_step = 1 fs and output_data_stride = 50, the saved-frame interval "
            "is 50 fs and the Nyquist limit is 10 THz. This setting describes the input "
            "sampling; it does not request additional downsampling.",
        ),
        related=("time_step", "total_num_steps", "num_blocks", "plot_cutoff_freq"),
        required_for_compute=True,
    ),
    "num_blocks": _Topic(
        "num_blocks = 5",
        "Divide the requested trajectory into equal time blocks, calculate one spectrum "
        "per block, and average the spectra.",
        notes=(
            "The requested frame count must be at least num_blocks and divisible by it. "
            "Each block must contain enough frames for the requested analysis.",
            "More blocks provide more averaging and smaller per-block arrays, but shorten "
            "the blocks and coarsen frequency spacing. Fewer blocks improve resolution "
            "while using larger working arrays.",
            "Approximate FFT spacing in THz: 1000 * num_blocks / "
            "(total_num_steps * time_step), with time_step in fs.",
        ),
        related=("total_num_steps", "output_data_stride", "initial_guess_hwhm"),
    ),
    "trajectory_file": _Topic(
        "trajectory_file = ../gpumd_run/dump.xyz",
        "Path to one trajectory containing atomic positions, velocities, and simulation "
        "cell information. The format is detected from the file contents.",
        notes=(
            "Supported text input is extended XYZ (.xyz), or one LAMMPS custom dump with "
            "atom IDs, x y z, vx vy vz, and box information. LAMMPS positions and "
            "velocities must be in the same file.",
            "Compatible NetCDF trajectories and existing pySED NetCDF caches can be read "
            "directly. Such input is not converted and does not create another cache.",
            "Relative paths are resolved from the directory containing the pySED input "
            "file. pySED reads exactly the named file; it does not select neighboring "
            ".xyz or .nc files automatically.",
        ),
        examples=(
            "trajectory_file = ../lammps_run/trajectory.dump",
            "trajectory_file = trajectory.nc",
        ),
        related=("trajectory_read_mode", "lammps_unit", "basis_lattice_file"),
    ),
    "basis_lattice_file": _Topic(
        "basis_lattice_file = basis.in",
        "Path to the atom mapping that identifies each atom's unit cell, basis index, "
        "and mass for the SED calculation.",
        notes=(
            "The mapping must describe the atoms and ordering used by the trajectory. "
            "For a dumped atom group, provide a mapping for that same group.",
            "Relative paths are resolved from the directory containing the pySED input "
            "file. Use the structure-generation tools to prepare a consistent mapping.",
        ),
        related=("num_atoms", "prim_unitcell", "supercell_dim", "dump_group_supercell"),
    ),
    "out_files_name": _Topic(
        "out_files_name = bulk_MoS2",
        "Prefix for SED data and figures, such as bulk_MoS2.SED, bulk_MoS2.Qpts, "
        "bulk_MoS2.THz, and bulk_MoS2-SED.png.",
        notes=(
            "A relative prefix is interpreted from the directory where pySED is launched, "
            "not from the directory containing the input file. Fitting-Qpoint/ and "
            "Lifetime/ are also created in the working directory.",
            "Use the same prefix when plotting or fitting a previously calculated result.",
        ),
        related=("action", "output_partial", "plot_partial_SED"),
    ),
    "lammps_unit": _Topic(
        "lammps_unit = metal",
        "Unit style used by a LAMMPS text dump, so pySED can convert its velocities "
        "and other quantities consistently.",
        values=(("metal", "Use LAMMPS metal units."), ("real", "Use LAMMPS real units.")),
        notes=(
            "Set this to the units style of the MD run. It applies to LAMMPS text input; "
            "NetCDF input uses its own recognized unit metadata.",
        ),
        related=("trajectory_file", "time_step"),
    ),
    "trajectory_read_mode": _Topic(
        "trajectory_read_mode = cache",
        "Choose how supported text trajectories reach the SED calculation.",
        values=(
            ("cache", "Create or reuse <source-name>.pysed.nc beside the pySED input file."),
            ("direct", "Read the text trajectory directly without creating a NetCDF cache."),
        ),
        notes=(
            "Existing compatible NetCDF input is always read directly; this setting has "
            "no effect on it. Text conversion retains the original text trajectory.",
            "A cache avoids repeated text parsing on later runs. direct saves the extra "
            "cache file but parses the text again for each calculation.",
            "Both modes read time blocks. RAM use depends on block size, prefetch, and "
            "calculation buffers. Cache compression changes disk use, not block size.",
        ),
        related=("trajectory_file", "netcdf_compression_level", "netcdf_batch_size", "num_blocks"),
    ),
    "trajectory_prefetch": _Topic(
        "trajectory_prefetch = True",
        "Prepare the next trajectory block in one background thread while the current "
        "block is used for SED computation.",
        notes=(
            "At most one block is prefetched; this is a Boolean switch, not a prefetch "
            "depth. Block order is unchanged. True/False and 1/0 are accepted.",
            "Prefetch applies to direct text, NetCDF caches, and native NetCDF input. "
            "It retains one extra raw block in CPU memory; disable it when memory "
            "is limited. A speedup depends on the trajectory and filesystem.",
        ),
        related=("num_blocks", "trajectory_read_mode", "device"),
    ),
    "netcdf_batch_size": _Topic(
        "netcdf_batch_size = 64",
        "Number of text-trajectory frames parsed per batch during cache conversion or "
        "direct text reading. Must be a positive integer.",
        notes=(
            "This is neither the on-disk NetCDF/HDF5 chunk size nor the SED averaging "
            "block length. Native NetCDF input is unaffected.",
            "Smaller batches reduce temporary parser memory. Larger batches may reduce "
            "batching overhead, but do not guarantee faster reads.",
        ),
        related=("trajectory_read_mode", "num_blocks"),
    ),
    "netcdf_compression_level": _Topic(
        "netcdf_compression_level = 1",
        "Compression level for a NetCDF cache created from a text trajectory. "
        "Accepted integer values are 0 through 9.",
        values=(
            ("0", "Disable compression."),
            ("1", "Use light compression."),
            ("2-9", "Request progressively stronger compression."),
        ),
        notes=(
            "Stronger compression can reduce disk use but generally slows conversion. "
            "Changing the level causes pySED to rebuild its converted cache.",
            "This setting has no effect for direct text mode or an already named "
            "NetCDF trajectory.",
        ),
        related=("trajectory_read_mode", "trajectory_file"),
    ),
    "output_partial": _Topic(
        "output_partial = True",
        "Save SED contributions resolved by element and Cartesian direction during "
        "computation, in addition to total SED.",
        notes=(
            "Partial files are written under <out_files_name>_partial_SED/. Enable this "
            "for action = compute, then select a component with plot_partial_SED.",
            "These are element/direction contributions, not phonon-eigenvector-resolved "
            "mode spectra. True/False and 1/0 are accepted.",
        ),
        related=("plot_partial_SED", "out_files_name", "legacy_unit"),
    ),
    "legacy_unit": _Topic(
        "legacy_unit = False",
        "Choose the intensity normalization for a new SED calculation. False uses "
        "the pySED 3.0 physical scale in eV/THz; True uses the historical 2.x scale.",
        notes=(
            "Total and partial output record their scale. Plotting and fitting use "
            "the saved metadata; changing this parameter does not rescale existing data.",
            "For old 2.x output without scale metadata, set True when plotting or fitting. "
            "To change the normalization of a computed result, rerun action = compute.",
            "Frequency remains in THz and lifetime in ps. This temporary option is for "
            "historical comparisons and will be removed in a future major release.",
        ),
        related=("action", "output_partial"),
    ),
    "max_cores": _Topic(
        "max_cores = 4",
        "Maximum number of CPU worker processes for SED computation. Must be positive.",
        notes=(
            "With device = cpu, a value of 1 runs serially and larger values enable "
            "multiprocessing. More processes may increase memory use as well as speed.",
            "This setting does not control the single-GPU calculation or the number "
            "of GPUs.",
        ),
        related=("device", "num_blocks"),
    ),
    "prim_unitcell": _Topic(
        "prim_unitcell = 3.11905 0 0  -1.55952 2.70117 0  0 0 12.6227",
        "Primitive-cell lattice vectors in Angstrom. Supply nine decimal numbers, "
        "with one lattice vector per row of the resulting 3 by 3 matrix.",
        notes=(
            "Use the primitive cell consistent with basis.in and supercell_dim. "
            "The matrix must be finite and nonsingular.",
            "Use decimals for this matrix, not fraction strings such as 1/3. Fraction "
            "syntax is supported for q_path only.",
        ),
        related=("basis_lattice_file", "supercell_dim", "prim_axis", "rescale_prim"),
        required_for_compute=True,
    ),
    "prim_axis": _Topic(
        "prim_axis = 0 0.5 0.5  0.5 0 0.5  0.5 0.5 0",
        "Optional 3 by 3 transformation from the input structure cell to the "
        "primitive-cell basis expected by pySED.",
        notes=(
            "Supply nine decimal numbers. Use this only when the structure is not "
            "already in the required primitive basis, and inspect the transformed "
            "primitive cell printed during the calculation.",
        ),
        related=("prim_unitcell", "supercell_dim", "rescale_prim"),
        default_note="No extra primitive-cell transformation.",
    ),
    "supercell_dim": _Topic(
        "supercell_dim = 20 20 20",
        "Replication counts along the three cell axes. Supply three positive integers "
        "consistent with the MD supercell and its basis mapping.",
        notes=(
            "Supercell geometry determines the commensurate q-point sampling. "
            "Increasing the physical MD supercell can provide denser allowed q-points.",
        ),
        related=("prim_unitcell", "basis_lattice_file", "q_path", "rescale_prim"),
    ),
    "rescale_prim": _Topic(
        "rescale_prim = True",
        "Reconstruct the primitive cell from the actual trajectory cell, for example "
        "after NPT relaxation has changed the cell dimensions.",
        notes=(
            "For a fixed-cell run, the trajectory should be consistent with the supplied "
            "primitive cell and supercell dimensions. True/False and 1/0 are accepted.",
            "Set False when using dump_group_supercell, and supply the matching "
            "prim_unitcell explicitly.",
        ),
        related=("prim_unitcell", "supercell_dim", "dump_group_supercell"),
    ),
    "dump_group_supercell": _Topic(
        "dump_group_supercell = 80.712 0 0  40.356 69.899 0  0 0 12.356",
        "Local cell represented by basis.in when a trajectory contains one atom group "
        "dumped from a larger simulation whose stored cell describes the full system.",
        notes=(
            "Supply nine decimal values in Angstrom, with lattice vectors as rows and "
            "the same Cartesian axes as the trajectory. This changes q-point geometry, "
            "not atomic positions or velocities.",
            "Set rescale_prim = False and provide a consistent prim_unitcell and local "
            "basis.in. Do not use this override for a whole-system trajectory whose "
            "mapping already describes the full simulation cell.",
        ),
        related=("basis_lattice_file", "num_atoms", "prim_unitcell", "rescale_prim"),
        default_note="Use the simulation cell stored in the trajectory.",
    ),
    "num_qpaths": _Topic(
        "num_qpaths = 3",
        "Number of consecutive q-path segments. A path G-M-K-G has three segments "
        "and four endpoints.",
        notes=(
            "q_path must contain num_qpaths + 1 coordinate triples, and q_path_name "
            "must provide the same number of labels.",
        ),
        related=("q_path", "q_path_name"),
    ),
    "q_path_name": _Topic(
        "q_path_name = GMKG",
        "Labels for consecutive path endpoints in the SED plot. Use G for Gamma.",
        notes=(
            "The number of labels must be num_qpaths + 1 and must match the endpoint "
            "triples in q_path.",
        ),
        related=("num_qpaths", "q_path"),
    ),
    "q_path": _Topic(
        "q_path = 0 0 0  1/2 0 0  1/3 1/3 0  0 0 0",
        "Consecutive path endpoints in reduced reciprocal coordinates. Supply "
        "num_qpaths + 1 triples; these define the path, not every sampled q-point.",
        notes=(
            "Explicit fractions such as 1/3 and 7/15 are preserved exactly. Decimals "
            "within 1e-5 of a fraction with denominator at most 12 are normalized to "
            "that fraction; other decimal values are retained as written.",
            "Only q-points compatible with the MD supercell periodicity are included. "
            "A supplied endpoint is not automatically an allowed q-point. A repeated "
            "endpoint contributes a point only if it is commensurate.",
        ),
        related=("num_qpaths", "q_path_name", "supercell_dim"),
        required_for_compute=True,
    ),
    "plot_cutoff_freq": _Topic(
        "plot_cutoff_freq = 20",
        "Maximum frequency displayed in the SED plot, in THz.",
        notes=(
            "This limits the plotted range; it does not change the saved spectrum "
            "or the trajectory sampling limit.",
            "For single-q fit figures, lorentz_fit_freq_max takes precedence when set. "
            "Otherwise the figure uses plot_cutoff_freq, or the available FFT range "
            "if neither is set. This display limit does not restrict fitting.",
        ),
        related=("plot_interval", "output_data_stride", "lorentz_fit_freq_max"),
        default_note="Use the available FFT frequency range.",
    ),
    "plot_interval": _Topic(
        "plot_interval = 5",
        "Spacing between frequency-axis ticks in the SED plot, in THz.",
        notes=("Choose an interval appropriate to plot_cutoff_freq. This is a display "
               "setting and does not set the FFT frequency resolution.",),
        related=("plot_cutoff_freq", "num_blocks"),
    ),
    "plot_color": _Topic(
        "plot_color = RdBu_r",
        "Matplotlib colormap used for the SED intensity map.",
        notes=("Other colormap names include jet, Spectral, and viridis. Tune "
               "colorbar_min and colorbar_max separately to adjust intensity contrast.",),
        related=("colorbar_min", "colorbar_max", "use_contourf"),
    ),
    "colorbar_min": _Topic(
        "colorbar_min = -22",
        "Lower limit of the dimensionless natural-log SED color scale. Increase it "
        "to hide weak background, or decrease it to reveal low-intensity branches.",
        notes=(
            "Negative values are valid. The value must be smaller than colorbar_max "
            "when both limits are specified. Plot once with automatic limits before tuning.",
            "The physical scale uses ln(SED / (1 eV/THz)). Saved scale metadata "
            "determines labels for physical versus historical results.",
        ),
        related=("colorbar_max", "plot_color", "legacy_unit"),
        default_note="Use the minimum natural-log SED value.",
    ),
    "colorbar_max": _Topic(
        "colorbar_max = -10",
        "Upper limit of the dimensionless natural-log SED color scale. Increase it "
        "if the strongest branches saturate, or decrease it to emphasize weak branches.",
        notes=(
            "Negative values are valid. The value must be greater than colorbar_min "
            "when both limits are specified.",
            "The physical scale uses ln(SED / (1 eV/THz)). Saved scale metadata "
            "determines labels for physical versus historical results.",
        ),
        related=("colorbar_min", "plot_color", "legacy_unit"),
        default_note="Use the maximum natural-log SED value.",
    ),
    "use_contourf": _Topic(
        "use_contourf = True",
        "Use filled contours for a single-path SED intensity map instead of the "
        "default image-style plot.",
        notes=(
            "Multiple q-paths already use filled contours automatically. plot_slice "
            "produces a separate line spectrum and is unaffected by this choice.",
            "True/False and 1/0 are accepted.",
        ),
        related=("num_qpaths", "plot_color", "plot_slice"),
    ),
    "plot_slice": _Topic(
        "plot_slice = True",
        "Plot the spectrum at the single q-point selected by qpoint_slice_index.",
        notes=(
            "Use action = plot to inspect individual peaks before fitting. The slice "
            "shows SED against frequency with a logarithmic intensity axis.",
            "True/False and 1/0 are accepted.",
        ),
        related=("action", "qpoint_slice_index", "plot_cutoff_freq", "peak_min_significance"),
    ),
    "qpoint_slice_index": _Topic(
        "qpoint_slice_index = 2",
        "Zero-based index of the q-point selected for a single-q spectrum or fit. "
        "Use the generated .Qpts file to identify the desired point.",
        notes=(
            "Index 0 selects the first generated q-point. The index must be nonnegative "
            "and smaller than the number of q-points in the saved SED data.",
            "This selects the fitting q-point when lorentz_fit_all_qpoint is False. "
            "With action = thinking, the fit is skipped if its output already exists; "
            "use action = fit to force a refit.",
        ),
        related=("plot_slice", "lorentz_fit_all_qpoint", "re_output_total_freq_lifetime"),
    ),
    "plot_partial_SED": _Topic(
        "plot_partial_SED = O\nplot_partial_SED = O y",
        "Select an element's partial SED for plotting or fitting. An element alone "
        "sums x, y, and z; adding x, y, or z selects that Cartesian component.",
        notes=(
            "Requires partial SED files previously generated with output_partial = True. "
            "Use an element symbol present in the matching basis.in.",
            "Omit this parameter to use total SED. The value is an element and optional "
            "direction, not an atom index or a Boolean switch.",
        ),
        related=("output_partial", "basis_lattice_file", "out_files_name"),
    ),
    "if_show_figures": _Topic(
        "if_show_figures = False",
        "Choose whether figures are displayed interactively in addition to being saved.",
        notes=(
            "Leave False for unattended or headless runs. True requires a working "
            "interactive plotting display. True/False and 1/0 are accepted.",
        ),
        related=("action", "plot_slice"),
    ),
    "fitting_function": _Topic(
        "fitting_function = auto",
        "Line shape used to fit each detected spectral peak.",
        values=(
            ("auto", "Fit Lorentz and velocity-DHO candidates and compare their AICc scores."),
            ("lorentz", "Use a Lorentzian peak with fitted half-width at half maximum."),
            ("dho", "Use the velocity-spectrum damped harmonic oscillator line shape."),
        ),
        notes=(
            "auto chooses Lorentz when AICc scores differ by less than two, or when the "
            "DHO candidate is in the weak-damping limit (HWHM / center < 0.05). "
            "Both candidates use the same local fitting interval.",
            "Overdamped or incomplete-width peaks need inspection. A linewidth-derived "
            "SED lifetime is not a strict phonon lifetime in the overdamped regime.",
        ),
        related=("peak_min_significance", "initial_guess_hwhm", "peak_max_hwhm"),
    ),
    "peak_min_significance": _Topic(
        "peak_min_significance = 4.0",
        "Minimum dimensionless local-noise significance for adaptive peak detection. "
        "It compares a peak's prominence with a robust local noise estimate.",
        notes=(
            "Increase the value to reject more weak/noisy candidates; decrease it to "
            "detect weaker peaks. This is not an absolute SED intensity threshold or "
            "a formal statistical confidence level.",
            "Detection uses a Hann-smoothed log spectrum and local noise estimation. "
            "The line-shape fit uses the original, unsmoothed linear SED.",
            "Inspect a selected q-point before applying the settings to every q-point.",
        ),
        related=("plot_slice", "qpoint_slice_index", "initial_guess_hwhm", "lorentz_fit_all_qpoint"),
    ),
    "initial_guess_hwhm": _Topic(
        "initial_guess_hwhm = 0.001",
        "Initial peak half-width at half maximum (HWHM), in THz, supplied to the "
        "nonlinear Lorentz or DHO optimizer.",
        notes=(
            "The guess must be positive and no greater than peak_max_hwhm. It guides "
            "optimization; it does not set or improve frequency resolution.",
            "Approximate FFT spacing in THz is 1000 * num_blocks / "
            "(total_num_steps * time_step). For 300000 steps, 1 fs, and 5 blocks, "
            "each block lasts 60 ps and spacing is about 0.0167 THz.",
            "First tune peak_min_significance to select the intended peaks. Then inspect "
            "the fitted curves and adjust the width guess if needed. A small guess "
            "cannot make a sub-resolution linewidth physically resolved.",
        ),
        related=("peak_max_hwhm", "peak_min_significance", "total_num_steps", "num_blocks"),
    ),
    "peak_max_hwhm": _Topic(
        "peak_max_hwhm = 1.0",
        "Upper bound on fitted peak HWHM, in THz. Must be positive and at least "
        "initial_guess_hwhm.",
        notes=(
            "The large default is a loose bound. Choose a justified bound if broad "
            "backgrounds or unresolved peaks lead to unreasonable fitted widths.",
        ),
        related=("initial_guess_hwhm", "fitting_function", "peak_min_significance"),
    ),
    "lorentz_fit_freq_min": _Topic(
        "lorentz_fit_freq_min = 2.0",
        "Lower frequency bound for peak detection and fitting, in THz. The bound "
        "is inclusive and must be nonnegative.",
        notes=(
            "If both bounds are supplied, this must be smaller than lorentz_fit_freq_max. "
            "The selected interval must contain at least three sampled frequency points.",
        ),
        examples=("lorentz_fit_freq_min = 2.0\nlorentz_fit_freq_max = 5.0",),
        related=("lorentz_fit_freq_max", "peak_min_significance"),
        default_note="Start from the lowest available frequency.",
    ),
    "lorentz_fit_freq_max": _Topic(
        "lorentz_fit_freq_max = 20",
        "Upper frequency bound for peak detection and fitting, in THz. The bound "
        "is inclusive and must be positive.",
        notes=(
            "If both bounds are supplied, this must be greater than lorentz_fit_freq_min. "
            "Omitting both bounds selects the complete available frequency interval.",
            "This also sets the upper limit of single-q fit figures. When omitted, those "
            "figures use plot_cutoff_freq or the available FFT range, while fitting "
            "still uses the highest available frequency.",
        ),
        related=("lorentz_fit_freq_min", "peak_min_significance", "initial_guess_hwhm"),
        default_note="Use the highest available frequency.",
    ),
    "lorentz_fit_all_qpoint": _Topic(
        "lorentz_fit_all_qpoint = True",
        "Fit all generated q-points instead of only qpoint_slice_index, using the "
        "selected fitting_function and peak-detection settings.",
        notes=(
            "Tune representative q-points first. With action = fit, all requested fits "
            "run again. With action = thinking, all-q fitting runs only when its "
            "combined lifetime output is missing.",
            "Per-q figures go to Fitting-Qpoint/ and tables to Lifetime/. All-q fitting "
            "also rebuilds Fitting-All-Qpoints.Fre_lifetime and the frequency-lifetime "
            "figure, so re_output_total_freq_lifetime has no additional effect.",
            "Lifetime tables contain frequency (THz) and linewidth-derived lifetime (ps), "
            "not branch assignments. Inspect any peaks listed in the terminal's final "
            "quality summary before using the lifetimes.",
        ),
        related=("qpoint_slice_index", "fitting_function", "peak_min_significance", "re_output_total_freq_lifetime"),
    ),
    "modulate_factor": _Topic(
        "modulate_factor = 0",
        "Number of sampled points trimmed from each end of a detected peak's "
        "local valley-bounded fitting interval. Must be a nonnegative integer.",
        notes=(
            "Use sparingly: trimming removes measured peak tails and may change the "
            "fitted linewidth. It does not introduce a background model or a "
            "multi-peak fit.",
        ),
        related=("peak_min_significance", "initial_guess_hwhm", "fitting_function"),
    ),
    "re_output_total_freq_lifetime": _Topic(
        "re_output_total_freq_lifetime = True",
        "Rebuild the combined lifetime table and frequency-lifetime figure after "
        "refitting one selected q-point.",
        notes=(
            "Use with action = fit and lorentz_fit_all_qpoint = False. Select the "
            "q-point through qpoint_slice_index, adjust fitting controls, and set "
            "this switch to replace its contribution in the summary.",
            "All-q fitting already rebuilds the summary automatically; this switch "
            "has no effect when lorentz_fit_all_qpoint = True.",
        ),
        related=("qpoint_slice_index", "lorentz_fit_all_qpoint", "peak_min_significance"),
    ),
}


def _canonical_name(name: str) -> str:
    """Look up a public topic; case changes are accepted, legacy aliases are not."""
    names = {key.casefold(): key for key in PARAMETERS}
    folded = name.casefold()
    if folded in names:
        return names[folded]
    matches = get_close_matches(folded, names, n=1, cutoff=0.6)
    suggestion = f" Did you mean '{names[matches[0]]}'?" if matches else ""
    raise ValueError(
        f"unknown parameter '{name}'.{suggestion} "
        "Run 'pysed -h' for help and the parameter manual."
    )


def _default_text(name: str, topic: _Topic) -> str:
    """Render actual defaults while explaining unset computation inputs."""
    if topic.required_for_compute:
        return "Required for computation (action = compute or thinking)."
    if name == "plot_partial_SED":
        return "Disabled. Omit this parameter to use total SED."
    value = PARAMETERS[name].default
    if hasattr(value, "tolist"):
        value = " ".join(str(item) for item in value.tolist())
    text = str(value)
    if topic.default_note:
        text += ". " + topic.default_note
    return text


def _paragraph(text: str, indent: int = 8) -> str:
    return textwrap.fill(
        text,
        width=_WIDTH,
        initial_indent=" " * indent,
        subsequent_indent=" " * indent,
        break_long_words=False,
        break_on_hyphens=False,
    )


def format_parameter_help(name: str) -> str:
    """Return one self-contained offline help page for a public parameter.

    Topic lookup follows PARAMETERS, so removed parameters cannot survive as
    hidden help aliases. No files, network access, or user input are required.
    """
    name = _canonical_name(name)
    topic = _TOPICS[name]
    sections = [
        f"    PARAMETER: {name}",
        "    SYNTAX:\n" + textwrap.indent(topic.syntax, " " * 8),
        "    DESCRIPTION:\n" + _paragraph(topic.meaning),
    ]
    if topic.values:
        rows = ["    VALUES:"]
        for value, meaning in topic.values:
            rows.append(textwrap.fill(
                meaning,
                width=_WIDTH,
                initial_indent=f"        {value:<12}",
                subsequent_indent=" " * 20,
                break_long_words=False,
                break_on_hyphens=False,
            ))
        sections.append("\n".join(rows))
    sections.append("    DEFAULT:\n" + _paragraph(_default_text(name, topic)))
    if topic.notes:
        sections.append("    NOTES:\n" + "\n\n".join(_paragraph(note) for note in topic.notes))
    if topic.examples:
        sections.append(
            "    EXAMPLES:\n" + "\n\n".join(
                textwrap.indent(example, " " * 8) for example in topic.examples
            )
        )
    if topic.related:
        sections.append("    RELATED PARAMETERS:\n" + _paragraph(", ".join(topic.related)))
    sections.append("    MANUAL:\n        " + MANUAL_URL)
    return "\n\n".join(sections) + "\n"


__all__ = ["MANUAL_URL", "format_parameter_help"]
