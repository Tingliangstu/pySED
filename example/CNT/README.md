# pySED Tutorial: CNT SED with GPUMD

This example computes the spectral energy density (SED) of a carbon nanotube
(CNT) using a GPUMD trajectory and pySED.

Reproduce this example before applying pySED to a new one-dimensional system.

---

## Workflow Map

`[1. Structure] -> [2. GPUMD trajectory] -> [3. pySED SED] -> [4. Plot/Fit lifetime]`

| Stage | Folder | Main files | Output |
|---|---|---|---|
| `[1. Structure]` | `structure/` | `POSCAR_CNT`, `generate_gpumd_xyz.py` | `model.xyz`, `basis.in` |
| `[2. GPUMD]` | `gpumd_run/` | `run.in`, `nep.txt` | `dump.xyz` |
| `[3. pySED]` | `SED/` | `input_SED.in` | `CNT.SED`, `CNT.Qpts`, `CNT.THz` |
| `[4. Plot/Fit]` | `SED/` | `input_SED.in` | `CNT-SED.svg`, lifetime files |

---

## [Step 1] -> Generate `model.xyz` and `basis.in`

Go to the structure folder:

```bash
cd example/CNT/structure
python generate_gpumd_xyz.py
```

The script uses `pySED.structure.generate_data.structure_maker`:

```python
from pySED.structure import generate_data

def generate_required_files(input_file, supercell):
    structure = generate_data.structure_maker(structure_file_name=input_file)
    structure.replicate_supercell(supercell=supercell)
    structure.write_xyz(filename='model.xyz', pbc="T T T")
    structure.write_lattice_basis_file()

if __name__ == "__main__":
    file_name = 'POSCAR_CNT'
    supercell = (1, 1, 160)
    generate_required_files(file_name, supercell)
```

Input structure support:

- POSCAR-style files are supported, as used here with `POSCAR_CNT`.
- `.xyz` files are also supported. For example, use
  `file_name = 'primitive_CNT.xyz'` if your starting structure is an XYZ file.

Generated files:

- `model.xyz`: GPUMD structure file.
- `basis.in`: pySED basis mapping file.

Copy or link `model.xyz` into `gpumd_run/` if your GPUMD run folder needs it.

---

## [Step 2] -> Run GPUMD and write `dump.xyz`

Go to the GPUMD folder:

```bash
cd ../gpumd_run
gpumd
```

The production part of `run.in` is:

```text
ensemble       nve
dump_exyz      8     1
run            500000
```

This writes an extended XYZ trajectory every 8 MD steps. The matching pySED
settings are:

```text
total_num_steps = 500000
time_step = 1
output_data_stride = 8
dump_xyz_file = '../gpumd_run/dump.xyz'
```

---

## [Step 3] -> Compute SED with pySED

Go to the SED folder:

```bash
cd ../SED
pysed input_SED.in
```

For the first run, use:

```text
plot_SED = 0
```

pySED writes:

- `CNT.SED`
- `CNT.Qpts`
- `CNT.THz`
- `CNT.Q_distances_and_labels`

Important input settings:

```text
num_atoms = 17920
supercell_dim = 1 1 160
q_path_name = 'GA'
q_path = 0.0 0.0 0.0  0.0 0.0 0.5
```

---

## [Step 4] -> Plot SED and fit lifetime

After the SED files exist, set:

```text
plot_SED = 1
plot_slice = 1
lorentz = 1
```

Tune the single-q-point fit first:

- `qpoint_slice_index`
- `peak_height`
- `peak_prominence`
- `initial_guess_hwhm`
- `peak_max_hwhm`

After a reasonable single-q-point fit, use:

```text
lorentz_fit_all_qpoint = 1
```

The all-q-point workflow writes `TOTAL-LORENTZ-Qpoints.Fre_lifetime`.

---

## Expected Result

![CNT SED Result](SED/CNT-SED.svg)

---

## Checklist

- [x] `num_atoms = 17920` matches `basis.in` and `dump.xyz`.
- [x] `supercell_dim = 1 1 160` matches the generated CNT supercell.
- [x] `output_data_stride = 8` matches `dump_exyz 8 1`.
- [x] The q-path `G-A` follows the tube axis.
