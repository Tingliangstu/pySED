# pySED Tutorial: Calculating Spectral Energy Density (SED) for a CNT Example with GPUMD

This guide walks you through using **pySED** to compute the **Spectral Energy Density (SED)** of a **Carbon Nanotube (CNT)** system entirely in **GPUMD**.

It is strongly recommended that you reproduce this tutorial **exactly** before applying the method to your own system.

---

## **Workflow Overview**

1. **Generate structure input file for gpumd and basis files for pySED**
2. **Run GPUMD to produce trajectory files (`dump.exyz`)**
3. **Run pySED to calculate the SED**
 
## Usage

- [x] **1. Goto the structure folder and using the Modify the generate_gpumd_xyz.py file** 

**`generate_gpumd_xyz.py`** Now only read [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) file.

run **`python generate_gpumd_xyz.py`** the [model.xyz](https://gpumd.org/gpumd/input_files/model_xyz.html) data and the basis.in will be generate.


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
    structure = generate_data.structure_maker(structure_file_name=input_file)
    
    # Generate supercell
    structure.replicate_supercell(supercell=supercell)

    # Write xyz files for gpumd
    structure.write_xyz(filename='model.xyz', pbc="T T T")
    
    # Write basis.in files for further used
    structure.write_lattice_basis_file()

if __name__ == "__main__":

    file_name = 'POSCAR_CNT'                              ## The in file for lammps
    supercell = (1, 1, 160)
    generate_required_files(file_name, supercell)

```

- [x] **2. Goto the gpumd_run folder and run GPUMD to generate the dump.xyz** 

Maybe one can modify the **run.in**

```python

potential      nep.txt
dftd3          pbe  12  6
velocity       300

ensemble npt_mttk temp 300 300 z 0 0
time_step      1
dump_thermo    20000
dump_position  50000
run            500000

######### Nose-Hoover ###############
ensemble       nvt_nhc   300    300   100
time_step      1
dump_thermo    20000
dump_position  50000
run            500000

ensemble       nve
dump_exyz      8     1
#dump_netcdf   10     1
run            500000

```

- [x] **3. Goto the SED folder and run pySED to get SED data** 

One can modify the **input_SED.in** and play with it.

The parameters in input_sed.in are explained in detail.

One can use  **pysed -h** to find more details.



