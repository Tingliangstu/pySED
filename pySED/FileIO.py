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

import os
import numpy as np

from pySED import OutputPaths
from pySED.Lorentz import hwhm_to_lifetime_ps


SED_HEADER = (
    "pySED spectral energy density; unit: eV/THz; "
    "each column corresponds to one q-point "
    "(same order as the .Qpts file)"
)

def write_output(phonons, params, BZ_lattice_info):
    sed_path, qpoints_path, frequency_path, labels_path = (
        OutputPaths.sed_data_files(params.out_files_name)
    )
    output_partial = getattr(params, 'output_partial', 0)
    if output_partial:
        # Overall SED = sum over element and direction.
        sed_total = phonons.sed_avg.sum(axis=(-2, -1))
        np.savetxt(
            sed_path,
            sed_total,
            header=SED_HEADER,
        )

        dir_labels = ['x', 'y', 'z']
        for t_idx, element in enumerate(phonons.type_symbols):
            for d_idx, d_lab in enumerate(dir_labels):
                out_path = OutputPaths.partial_sed_data(
                    params.out_files_name,
                    element,
                    d_lab,
                )
                out_path.parent.mkdir(parents=True, exist_ok=True)
                np.savetxt(
                    out_path,
                    phonons.sed_avg[:, :, t_idx, d_idx],
                    header=SED_HEADER,
                )
    else:
        np.savetxt(
            sed_path,
            phonons.sed_avg,
            header=SED_HEADER,
        )

    np.savetxt(qpoints_path, BZ_lattice_info.reduced_qpoints, fmt='%.8f')
    np.savetxt(frequency_path, phonons.freq_fft, fmt='%.8f')

    with open(labels_path, 'w') as f:
        # q_distances
        f.write("Global distances along the paths:\n")
        f.write(" ".join(f"{d:.10f}" for d in BZ_lattice_info.q_distances) + "\n\n")

        # q_labels
        f.write("High-symmetry points and their distances:\n")
        for distance, label in BZ_lattice_info.q_labels:
            f.write(f"{float(distance):.10f}   {label}\n")

class load_data(object):

    def __init__(self, params):
        sed_path, qpoints_path, frequency_path, labels_path = (
            OutputPaths.sed_data_files(params.out_files_name)
        )

        if getattr(params, 'plot_partial_SED', 0):

            element = params.plot_partial_element

            if params.plot_partial_dir is None:
                # Sum x/y/z for the requested element.
                file_x = OutputPaths.partial_sed_data(params.out_files_name, element, 'x')
                file_y = OutputPaths.partial_sed_data(params.out_files_name, element, 'y')
                file_z = OutputPaths.partial_sed_data(params.out_files_name, element, 'z')
                if not all(
                    os.path.exists(path)
                    for path in (file_x, file_y, file_z)
                ):
                    raise FileNotFoundError(
                        f"partial SED files for element '{element}' are missing")
                sed_x = np.loadtxt(file_x)
                sed_y = np.loadtxt(file_y)
                sed_z = np.loadtxt(file_z)
                self.sed_avg = sed_x + sed_y + sed_z
            else:
                d = params.plot_partial_dir
                file_d = OutputPaths.partial_sed_data(params.out_files_name, element, d)
                if not os.path.exists(file_d):
                    raise FileNotFoundError(
                        f"partial SED file for element '{element}' "
                        f"and direction '{d}' is missing")
                self.sed_avg = np.loadtxt(file_d)

        else:
            self.sed_avg = np.loadtxt(sed_path)

        self.qpoints = np.loadtxt(qpoints_path)
        self.freq_fft = np.loadtxt(frequency_path)

        self.q_distances = []
        self.q_labels = {}
        with open(labels_path, 'r') as f:
            lines = f.readlines()
            # q_distances
            if lines[0].strip() == "Global distances along the paths:":
                self.q_distances = [float(x) for x in lines[1].split()]
            # q_labels
            if lines[3].strip() == "High-symmetry points and their distances:":
                for line in lines[4:]:
                    distance, label = line.split(maxsplit=1)
                    self.q_labels[float(distance)] = label.strip()


def write_phonon_lifetime(lorentz, params):
    """Write one Q-point's frequency and linewidth-derived lifetime."""

    output_file = OutputPaths.qpoint_lifetime_data(
        params.qpoint_slice_index
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_lines = [
        "Generated by pySED, Email: liangting.zj@gmail.com",
        (
            "Frequency (THz)  linewidth-derived SED lifetime (ps); "
            "tau_SED = 1/(2*pi*HWHM_THz)"
        ),
    ]

    lifetimes = hwhm_to_lifetime_ps(lorentz.popt[:, 2])
    for fit_parameters, lifetime in zip(lorentz.popt, lifetimes):
        if not np.isfinite(lifetime):  # don't output failed fits
            continue
        output_lines.append(
            f"{fit_parameters[0]:.6f} {lifetime:.8f}"
        )

    output_file.write_text(
        "\n".join(output_lines) + "\n",
        encoding="utf-8",
    )
    return str(output_file)


def read_lifetime_data(path):
    """Read a two-column pySED frequency-lifetime data file."""

    with open(path, "r", encoding="utf-8") as stream:
        data_lines = [
            line.strip() for line in stream.readlines()[2:] if line.strip()
        ]
    if not data_lines:
        return np.empty(0, dtype=float), np.empty(0, dtype=float)

    values = np.atleast_2d(np.loadtxt(data_lines, dtype=float))
    if values.shape[1] != 2:
        raise ValueError(f"invalid frequency-lifetime data in {path}")
    return values[:, 0], values[:, 1]


def deal_total_fre_lifetime(params, total_qpoints):
    """Combine all Q-point frequency-lifetime files into one data file."""

    output_lines = [
        "Generated by pySED, Email: liangting.zj@gmail.com",
        (
            "Frequency (THz)  linewidth-derived SED lifetime (ps); "
            "tau_SED = 1/(2*pi*HWHM_THz)"
        ),
    ]
    total_num_fre_lifetime = 0

    for i in range(total_qpoints):
        load_file_name = OutputPaths.qpoint_lifetime_data(i)
        try:
            frequency, lifetime = read_lifetime_data(load_file_name)
        except OSError as error:
            raise FileNotFoundError(
                f"cannot read Q-point lifetime data: {load_file_name}"
            ) from error

        if not frequency.size:
            print(
                f"\n  Warning: {load_file_name} contains no fitted "
                "lifetimes; skipping."
            )
            continue

        output_lines.extend(
            f"{freq:.6f} {tau:.8f}"
            for freq, tau in zip(frequency, lifetime)
        )
        total_num_fre_lifetime += frequency.size

    output_file = OutputPaths.combined_lifetime_data()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(
        "\n".join(output_lines) + "\n",
        encoding="utf-8",
    )
    return str(output_file), int(total_num_fre_lifetime)
