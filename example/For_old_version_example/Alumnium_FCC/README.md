# Example to use the pySED to calculate the SED of Alumnium_FCC
 
It is best to reproduce the results of this case before calculating one's system.
 
## Usage

- [x] **1. Goto the structure folder and using the Modify the generate_lammps_data.py file** 

**`generate_lammps_data.py`** Now only read [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) file.

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

    file_name = 'POSCAR'                              ## The in file for lammps
    supercell = (30, 2, 2)
    generate_required_files(file_name, supercell)
```

- [x] **2. Goto the lammps_run folder and run lammps to generate the pos.dat and vels.dat** 

Maybe one can modify the **in.vels**

```python
## *************** For SED run *****************
atom_style      atomic
units           metal
dimension       3
boundary        p p p

read_data       lammps_data

# ---------- Define Interatomic Potential --------------------- 
pair_style      eam/alloy 
pair_coeff      *  *  Al99.eam.alloy  Al

#neighbor       2.0    nsq
#neigh_modify   every  1  delay  0  check  yes

timestep        0.00075				                 ## 0.75 fs

variable        T equal 300
velocity        all create ${T} 163456 mom yes rot yes dist gaussian
fix             NVT  all   nvt  temp  ${T}   ${T}  0.1 

# *********************  thermo output  ************************
thermo	        2000
thermo_style    custom step temp lx ly lz press vol pxx pyy pzz pe ke etotal
run           	1000000
unfix           NVT

fix             NVE all nve
variable        dt_dump   equal 32
dump            vels  all  custom  ${dt_dump}  vels.dat  id  type  vx  vy  vz
dump_modify     vels  format  line "%d  %d  %0.8g  %0.8g  %0.8g"
dump_modify     vels  sort  id
dump            pos   all  custom  ${dt_dump}  pos.dat   id  type  x  y  z
dump_modify     pos   format  line  "%d  %d  %0.8g  %0.8g  %0.8g"
dump_modify     pos   sort  id

run             2097152 			# 2**21
unfix           NVE

```

- [x] **3. Goto the SED folder and run pySED to get SED data** 

One can modify the **input_SED.in** and play with it.

The parameters in input_sed. in are explained in detail.

