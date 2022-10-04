'''
@author:
**************************  LiangTing ***************************
                      liangting.zj@gmail.com
************************ 2021/4/22 23:03:21 *********************
'''

import numpy as np
from pySED.structure import read_structure

class structure_maker(object):

    def __init__(self, structure_file_name = 'POSCAR',
                       lammps_data_flag = False,
                       lammps_infile_name = None,
                       description = 'Create LAMMPS data for SED method by using LT_Codes'):
        """
        Generate the required file for SED projects (including lammps data for running lammps and basis.in for Reciprocal space)
        :param structure_file_name: input file for primitive structure (POSCAR file is recommended).
        :param lammps_data_flag: using the lammps interface to read input file (not recommended).
        :param lammps_infile_name: if lammps_data_flag, this file is needed.
        :param description: Some description.
        """
        self.description = str(description)
        self.lammpsInfile = lammps_infile_name
        self.lammps_data_flag = lammps_data_flag
        self.input_structure = structure_file_name

        if self.lammps_data_flag:
            structure = read_structure.get_structure_from_lammps(self.lammps_infile_name)
            self.unitcell_positions = structure['positions']
            self.cell = structure['cell']

        else:
            structure = read_structure.read_from_file_structure_poscar(file_name=self.input_structure)
            self.unitcell_positions = structure['scaled_positions']
            self.cell = structure['direct_cell']
            # *********************** For VASP POSCAR *******************
            for i in range(self.unitcell_positions.shape[0]):
                self.unitcell_positions[i, :] = (
                        self.unitcell_positions[i, 0] * self.cell[0, :] +
                        self.unitcell_positions[i, 1] * self.cell[1, :] +
                        self.unitcell_positions[i, 2] * self.cell[2, :])

        self.unitcell_num_atoms = self.unitcell_positions.shape[0]

        self.masses = structure['masses']
        self.basis_atoms_symbols = structure['symbols']

        self.unitcell_index = np.ones((self.unitcell_num_atoms, 1))
        self.basis_index = np.arange(1, self.unitcell_num_atoms + 1).reshape(self.unitcell_num_atoms, 1)

    def replicate_supercell(self, supercell = (1, 1, 1), origin_position = False):

        self.supercell = supercell
        unitcells = np.copy(self.unitcell_index)                # deep copy
        basis_atom_index = np.copy(self.basis_index)            # deep copy
        basis_atom_masses = np.copy(self.masses)
        basis_atom_types = np.copy(self.basis_atoms_symbols)

        position = np.copy(self.unitcell_positions[:, :])

        # Generate supercell position
        for i in range(self.supercell[0]):
            a1 = np.copy(self.cell[0, :]) * i
            for j in range(self.supercell[1]):
                a2 = np.copy(self.cell[1, :]) * j
                for k in range(self.supercell[2]):
                    a3 = np.copy(self.cell[2, :]) * k
                    # Skip 0 indexes
                    if i == 0 and j == 0 and k == 0:
                        continue
                    new_cell = np.copy(self.unitcell_positions)
                    new_cell[:, 0] = new_cell[:, 0] + a1[0] + a2[0] + a3[0]     # x
                    new_cell[:, 1] = new_cell[:, 1] + a1[1] + a2[1] + a3[1]     # y
                    new_cell[:, 2] = new_cell[:, 2] + a1[2] + a2[2] + a3[2]     # z

                    # For unitcell index
                    unitcells = unitcells + 1
                    self.unitcell_index = np.append(self.unitcell_index, unitcells, axis=0)
                    self.basis_index = np.append(self.basis_index, basis_atom_index, axis=0)

                    self.basis_atoms_symbols = np.append(self.basis_atoms_symbols, basis_atom_types, axis=0)
                    self.masses = np.append(self.masses, basis_atom_masses, axis=0)

                    position = np.append(position, new_cell, axis=0)

        # For origin_position due to boundary origined from 0
        if origin_position:
            position[:, 0] = position[:, 0] - position[:, 0].min()
            position[:, 1] = position[:, 1] - position[:, 1].min()
            position[:, 2] = position[:, 2] - position[:, 2].min()

        self.supercell_position = position
        self.num_atoms = position.shape[0]

    def write_xyz(self, filename = 'data.xyz'):
        with open(filename, 'w') as fid:
            fid.write('{0:d}\n{1}\n'.format(self.num_atoms, self.description))
            for i, row in enumerate(self.supercell_position):
                fid.write('{0} {1:20.10f} {2:20.10f} {3:20.10f}\n'.format(self.basis_atoms_symbols[i],
                                                                          row[0], row[1], row[2]))

        fid.close()
        print('************* ' + filename + ' is written successfully' + ' ************\n')

    def write_lammps_data(self, file_name = 'lammps_data', lammps_data_types = None):

        print("Warning: This version codes can only create LAMMPS file for orthogonal lattice vectors\n")
        # Get lammps data boundary
        self.boundary = np.dot(self.cell, np.diag(self.supercell))

        # For atom types
        self.unique_types = list(np.unique(self.basis_atoms_symbols))
        self.masses_for_lammps_output = [read_structure.symbol_to_mass(i) for i in self.unique_types]

        self.num_types = len(self.unique_types)

        # Get lammps unit for generate lammps data
        if lammps_data_types is None and self.lammps_data_flag:
            self.unit = read_structure.get_units(self.input_file)
            print('Lammps data units is ' + self.unit)

        elif lammps_data_types is not None:
            self.unit = lammps_data_types
            print('************* Lammps data is for {0} units *************'.format(self.unit))

        else:
            self.unit = 'real'
            print('Lammps data units is (using the default unit) ' + self.unit)

        # write to file
        with open(file_name, 'w') as fid:
            fid.write(self.description + '\n\n')
            fid.write('{0} atoms\n\n'.format(self.num_atoms))
            fid.write('{0} atom types\n\n'.format(self.num_types))
            # Boundary
            fid.write('{0:20.10f} {1:20.10f} xlo xhi\n'.format(0, np.linalg.norm(self.boundary[0])))
            fid.write('{0:20.10f} {1:20.10f} ylo yhi\n'.format(0, np.linalg.norm(self.boundary[1])))
            fid.write('{0:20.10f} {1:20.10f} zlo zhi\n\n'.format(0, np.linalg.norm(self.boundary[2])))
            # Masses
            fid.write('Masses\n\n')
            for i in range(self.num_types):
                fid.write('{0} {1:20.10f}  # {2}\n'.format(i + 1, self.masses_for_lammps_output[i], self.unique_types[i]))

            fid.write('\nAtoms\n\n')

            if self.unit == 'metal':
                for i, row in enumerate(self.supercell_position):
                    fid.write('{0} {1:d} {2:20.10f} {3:20.10f} {4:20.10f}  # {5}\n'.format(i + 1,
                                                                     self.unique_types.index(self.basis_atoms_symbols[i])+1,
                                                                     row[0], row[1], row[2],
                                                                     self.basis_atoms_symbols[i]))

            elif self.unit == 'real':
                for i, row in enumerate(self.supercell_position):
                    fid.write('{0} {1} {2:d} {3} {4:20.10f} {5:20.10f} {6:20.10f}  # {7}\n'.format(i + 1, 1, # ids and molecular's ids
                                                                     self.unique_types.index(self.basis_atoms_symbols[i]) + 1, # atom types
                                                                     0.0,   # charge
                                                                     row[0], row[1], row[2],
                                                                     self.basis_atoms_symbols[i]))

            else:
                print('This version codes can\'t deal with more lammps unit')
                exit()

        fid.close()

        print('\n************* ' + file_name + ' is written successfully' + ' ************\n')

    ## Most important file
    def write_lattice_basis_file(self, file_name='basis.in'):
        basis_file = self.description + '\n'
        basis_file += "atoms_ids unitcell_index basis_index mass_types\n"
        for i in range(self.num_atoms):
            basis_file += '{0} {1:d} {2:d} {3:10.6f}\n'.format(i+1, int(self.unitcell_index[i][0]),
                                                               int(self.basis_index[i][0]), self.masses[i])
        # write the file
        f = open(file_name, 'w')
        f.write(basis_file)
        f.close()

        print('************* ' + file_name + ' is written successfully' + ' ************\n')

if __name__ == "__main__":

    file_name = 'HKUST-1_cubic.vasp'                              ## The in file for lammps

    si = structure_maker(structure_file_name = file_name)
    si.replicate_supercell(supercell = (1, 1, 60))
    si.write_xyz()
    si.write_lammps_data(lammps_data_types = 'metal')
    si.write_lattice_basis_file()







