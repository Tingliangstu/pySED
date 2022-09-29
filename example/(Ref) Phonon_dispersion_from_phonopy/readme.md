# Phonon Dispersion Calculator
Using Phonopy and Lammps-python interface to calculate the phonon dispersion and can be use to calculate most of different systems.

Need phonopy module and lammps-python interface installed.
 
## Files including in project

- [x] get_2or_force_constant.py (main program)
- [x] get_structure.py (get input structure from lammps) 
- [x] phonopy_link.py (call phonopy to calculate the phonon dispersion)
- [x] seekpath_link.py (call seekpath)
- [x] arrange.py

## Files need to be modified (modify only two files)

- [x] get_2or_force_constant.py (main program)
- [x] in.lammps file

## Usage
At the moment I don't want to encapsulate these files as modules, so one need to use them in the same folder

**According to one specific needs, modify the potential format in in.lammps**

- [x] **1. Modify the in.lammps file** 

```python
#read graphite  
read_data  				data.si   # Need to modify according the one's systems

#potential        function
pair_style        tersoff
pair_coeff        * *  SiCGe.tersoff  Si(D)
```

- [x] **2. Modify the get_2or_force_constant.py file** 

```python

in_file = 'in.lammps'  ## The in file for lammps

force_constants = Force_2or_lammps(in_file,
                                   supercell_matrix=np.diag([4, 4, 4]),
                                   primitive_matrix=[[1, 0.0, 0.0],
                                                     [0.0, 1, 0.0],
                                                     [0.0, 0.0, 1]])

# force_constants.my_write_force_sets()
force_constants.write_force_constants()

band_path = [[0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0]]   # Not need

band_label = [('GAMMA', 'X')]  #, ('X', 'S'), ('S', 'Y'), ('Y', 'GAMMA')]

force_constants.plot_phonon_dispersion_bands(use_seek_path=False,
                                             band_path=band_path,
                                             band_label = band_label,   # Auto determined
                                             band_resolution=101)


```

- [x] **3. python get_2or_force_constant.py** 