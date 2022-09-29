'''
@author:
***************************  LiangTing ************************
                     liangting.zj@gmail.com
************************ 2020/8/28 14:22:15 *******************
'''

import sys
import numpy as np
from phonopy.structure.atoms import PhonopyAtoms

def mass_to_symbol(mass, tolerance=5e-2):
    from phonopy.structure.atoms import atom_data

    for element in atom_data:
        if element[3] is not None and abs(mass - element[3]) < tolerance:
            return element[1]

    return 'H'           # If no match is found, use H as the wildcard

# **************************** should be POSCAR-file *************************
def get_structure_from_lammps(lammps_input_file, show_log = False):
    """
    Get the crystal structure from lammps input
    :param command_list: LAMMPS input file (One or more lines)
    :return: numpy array matrix with forces of atoms [Natoms x 3]
    """

    from lammps import lammps

    cmd_list = ['-log', 'generate_structure.log']
    if not show_log:
        cmd_list += ['-echo', 'none', '-screen', 'none']

    lmp = lammps(cmdargs = cmd_list)

    if lammps_input_file is not None:
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
    natoms = lmp.get_natoms()   # Total # of atoms as int

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

    unitcell = np.array([[xhi-xlo, xy,  xz],
                           [0,  yhi-ylo,  yz],
                           [0,   0,  zhi-zlo]]).T

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
    coords_data = lmp.gather_atoms("x", 1, 3)     # The atomic coordinates

    masses = lmp.extract_atom("mass", 2)          # atomic  mass (vector of doubles)
    type_data = lmp.gather_atoms("type", 0, 1)    # atomic  type

    positions = np.array(coords_data[:], dtype=np.dtype('double')).reshape((natoms, 3))
    masses = np.array([masses[type_data[i]] for i in range(natoms)])   # Corresponds to the mapping in the data file

    symbols = [mass_to_symbol(masses[i]) for i in range(natoms)]

    return PhonopyAtoms(positions = positions,
                       masses = masses,
                       symbols = symbols,
                       cell = unitcell)

if __name__ == "__main__":

    in_file = 'in.lammps'                            ## The in file for lammps
    get_structure_from_lammps(in_file, False)
