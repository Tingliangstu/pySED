'''
@author:
**************************  LiangTing ***************************
                      liangting.zj@gmail.com
************************ 2020/9/22 14:22:15 *********************
'''
import os, sys
import numpy as np
from pySED.structure.atoms import atom_data

def read_from_file_structure_poscar(file_name='POSCAR', number_of_dimensions = 3):
    """
    Read crystal structure from a VASP POSCAR type file
    :param file_name: POSCAR filename
    :param number_of_dimensions: number of dimensions of the crystal structure
    :return: Atoms (like phonopy) type object containing the crystal structure
    """

    # Check file exists
    if not os.path.isfile(file_name):
        print('POSCAR file does not exist!')
        exit()

    # Read from VASP POSCAR file
    print("\n********************* Reading From VASP POSCAR structure *********************")
    poscar_file = open(file_name, 'r')
    data_lines = poscar_file.read().split('\n')
    poscar_file.close()

    multiply = float(data_lines[1])      # Scale factor

    direct_cell = np.array([data_lines[i].split()
                            for i in range(2, 2 + number_of_dimensions)], dtype=float)
    direct_cell *= multiply
    scaled_positions = None
    positions = None

    try:
        # Line 7
        number_of_types = np.array(data_lines[3 + number_of_dimensions].split(), dtype=int)  # The number of atomic types
        # Line 8
        coordinates_type = data_lines[4 + number_of_dimensions][0]      # cartesian coordinates or direct (Fractional coordinates)
        # Line 8
        if coordinates_type == 'D' or coordinates_type == 'd':
            scaled_positions = np.array([data_lines[8 + k].split()[0:3]
                                         for k in range(sum(number_of_types))], dtype=float)
        else:
            positions = np.array([data_lines[8 + k].split()[0:3]
                                  for k in range(sum(number_of_types))], dtype=float)

        atomic_types = []
        for i, j in enumerate(data_lines[5].split()):
            atomic_types.append([j] * number_of_types[i])

        atomic_types = [item for sublist in atomic_types for item in sublist]

    # Old style POSCAR format
    except ValueError:
        print("\n********************* Reading Old Style POSCAR *********************")    # In vasp 4.X version
        number_of_types = np.array(data_lines[5].split(), dtype=int)
        coordinates_type = data_lines[6][0]

        # Line 7
        if coordinates_type == 'D' or coordinates_type == 'd':
            scaled_positions = np.array([data_lines[7 + k].split()[0:3]
                                         for k in range(sum(number_of_types))], dtype=float)
        else:
            positions = np.array([data_lines[7 + k].split()[0:3]
                                  for k in range(sum(number_of_types))], dtype=float)

        atomic_types = []
        for i, j in enumerate(data_lines[0].split()):
            atomic_types.append([j] * number_of_types[i])

        atomic_types = [item for sublist in atomic_types for item in sublist]

    masses = [symbol_to_mass(atomic_types[i]) for i in range(len(atomic_types))]

    structure = {'direct_cell': direct_cell,         # cell_matrix, lattice vectors in rows
                 'cell': None,
                 'scaled_positions': scaled_positions,
                 'positions': positions,
                 'masses': masses,
                 'symbols': atomic_types
                 }

    return structure

def mass_to_symbol(mass, tolerance=5e-2):

    for element in atom_data:
        if element[3] is not None and abs(mass - element[3]) < tolerance:
            return element[1]

    return 'H'           # If no match is found, use H as the wildcard

def symbol_to_mass(symbol):

    for element in atom_data:
        if element[1] is not None and symbol==element[1]:
            return element[3]

    return 1           # If no match is found, use H mass as the wildcard

def get_units(lammps_input_file):
    """
    Get the units label for LAMMPS "units" command from a list of LAMMPS input commands
    :param lammps_input_file: LAMMPS input file name (strings)
    :return units: string containing the units
    """
    lammps_commands_list = open(lammps_input_file).read().split('\n')

    for line in lammps_commands_list:
        if line.startswith('units'):
            return line.split()[1]

    return 'lj'

def get_structure_from_lammps(lammps_input_file, show_log = False):
    """
    Get the crystal structure from lammps input, if using the function one need a installed python-lammps interface on your machine,
    so it is not recommended. This function maybe only use for test.
    :param command_list: LAMMPS input file (One or more lines)
    :return: numpy array matrix with forces of atoms [Natoms x 3]
    """

    print("\n********************* Reading From LAMMPS Data **********************")
    from lammps import lammps       # Must have a installed lammps-python interface on your machine

    cmd_list = ['-log', 'read_structure_from_lammps.log']
    if not show_log:
        cmd_list += ['-echo', 'none', '-screen', 'none']

    lmp = lammps(cmdargs=cmd_list)

    if os.path.exists(lammps_input_file):
        pass
    else:
        sys.exit('\n\tFile error: The in file that lammps needs to read does not exist!\n')

    # lines = open(in_lammps,'r').readlines()
    # for line in lines: self.lmp.command(line)
    """
    Do the same as lmp.command(filename) but don't allow the script to be
    continued after quit.
    """
    lines = open(lammps_input_file, 'r').readlines()
    for line in lines:
        # print(line)
        if "quit" in line and line[0] != "#":
            sys.exit('\n\tQuit error: It is not appropriate '
                     'to exit the execution of the in.file at this time\n')  ## Quit lammps if find the keyword 'quit'
        else:
            lmp.command(line)

    lmp.command('run 0')

    natoms = lmp.get_natoms()  # Total # of atoms as int
    # Box information
    try:
        xlo = lmp.extract_global("boxxlo")
        xhi = lmp.extract_global("boxxhi")
        ylo = lmp.extract_global("boxylo")
        yhi = lmp.extract_global("boxyhi")
        zlo = lmp.extract_global("boxzlo")
        zhi = lmp.extract_global("boxzhi")
        xy = lmp.extract_global("xy")
        yz = lmp.extract_global("yz")
        xz = lmp.extract_global("xz")
    except TypeError:
        xlo = lmp.extract_global("boxxlo", 1)
        xhi = lmp.extract_global("boxxhi", 1)
        ylo = lmp.extract_global("boxylo", 1)
        yhi = lmp.extract_global("boxyhi", 1)
        zlo = lmp.extract_global("boxzlo", 1)
        zhi = lmp.extract_global("boxzhi", 1)
        xy = lmp.extract_global("xy", 1)
        yz = lmp.extract_global("yz", 1)
        xz = lmp.extract_global("xz", 1)

    except UnboundLocalError:
        boxlo, boxhi, xy, yz, xz, periodicity, box_change = lmp.extract_box()
        xlo, ylo, zlo = boxlo
        xhi, yhi, zhi = boxhi

    unitcell = np.array([[xhi - xlo, xy, xz],
                         [0, yhi - ylo, yz],
                         [0, 0, zhi - zlo]]).T

    # ************************ For lammps-python interface **********************
    '''     
    data = lmp.gather_atoms(name,type,count)  ----->>> get data information 
    count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f       
    type = 0 for integer values, 1 for double values

    coords = lmp.extract_atom(name,type)  # extract a per-atom quantity
                                          # name = "x", "type", etc
                                          # type = 0 = vector of ints
                                          #        1 = array of ints
                                          #        2 = vector of doubles
                                          #        3 = array of doubles
    '''

    # Coordinates ordered by atom ID
    coords_data = lmp.gather_atoms("x", 1, 3)  # The atomic coordinates

    masses = lmp.extract_atom("mass", 2)  # atomic  mass (vector of doubles)
    type_data = lmp.gather_atoms("type", 0, 1)  # atomic  type

    positions = np.array(coords_data[:], dtype=np.dtype('double')).reshape((natoms, 3))
    masses = np.array([masses[type_data[i]] for i in range(natoms)])  # Corresponds to the mapping in the data file

    symbols = [mass_to_symbol(masses[i]) for i in range(natoms)]

    structure = {'direct_cell': None,  # cell_matrix, lattice vectors in rows
                 'cell': unitcell,
                 'scaled_positions': None,
                 'positions': positions,
                 'masses': masses,
                 'symbols': symbols
                 }

    return structure

if __name__ == "__main__":

    file_name = 'HKUST-1_cubic.vasp'                              ## The in file for lammps

    # structure = get_structure_from_lammps(file_name)

    structure = read_from_file_structure_poscar(file_name=file_name)
    print(structure['masses'])


