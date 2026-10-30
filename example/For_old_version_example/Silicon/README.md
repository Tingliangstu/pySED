# 🧩 Legacy Example: Silicon SED with pySED

It is best to reproduce the results of this case before calculating one's system.

This is an older LAMMPS-based workflow. Use it together with the
[online manual](https://pysed.readthedocs.io/en/latest/) and check the
[`input.in` parameter guide](https://pysed.readthedocs.io/en/latest/input_parameters.html)
when editing pySED inputs.

## Usage

- [x] **1. Goto the structure folder and using the Modify the generate_lammps_data.py file**

**`generate_lammps_data.py`** uses a [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR)
file in this legacy example. The current `structure_maker` can also read `.xyz`
input structures.

run **`python generate_lammps_data.py`** the lammps data and the basis.in will be generate

```python
from pySED.structure import generate_data

def generate_required_files(input_file, supercell):
    '''
    structure_maker class include functions:
    1.replicate_supercell
    2.write_xyz
    3.write_lammps_data
    4.write_lattice_basis_file
    '''
    # Generate a structure class
    structure = generate_data.structure_maker(structure_file_name = file_name)

    # Generate supercell
    structure.replicate_supercell(supercell = supercell)

    # Write xyz files for view
    structure.write_xyz()

    # Write lammps data files for run
    structure.write_lammps_data(lammps_data_types = 'metal')

    # Write basis.in files for further used
    structure.write_lattice_basis_file()

if __name__ == "__main__":

    file_name = 'POSCAR'                     ## The in file for vasp ()
    supercell = (16, 2, 2)
    generate_required_files(file_name, supercell)
```

- [x] **2. Goto the lammps_run folder and run lammps to generate trajectory.dump**

Maybe one can modify the **in.lammps**

```python

## *************** For SED run *****************
atom_style      atomic
units           metal
dimension       3
boundary        p p p

read_data       lammps_data
pair_style      tersoff
pair_coeff      *  *  SiCGe.tersoff Si(D)

timestep        0.00075				                 ## 0.75 fs

variable        T equal 500
velocity        all create ${T} 123456 mom yes rot yes dist gaussian
fix             NVT  all   nvt  temp  ${T}   ${T}  0.1

# *********************  thermo output  ************************
thermo	        2000
thermo_style    custom step temp lx ly lz press vol pxx pyy pzz pe ke etotal
run           	1000000
unfix           NVT

fix             NVE all nve
variable        dt_dump   equal 32
dump            pysed all custom ${dt_dump} trajectory.dump id type x y z vx vy vz
dump_modify     pysed sort id

run             2097152 			# 2**21
unfix           NVE


```

- [x] **3. Goto the SED folder and run pySED to get SED data**

One can modify the **input.in** and play with it.

The parameters in **input.in** are explained in the
[`input.in` parameter guide](https://pysed.readthedocs.io/en/latest/input_parameters.html).


For pySED 3.0, set `action = compute` and
`trajectory_file = ../lammps_run/trajectory.dump` in `SED/input.in`.
