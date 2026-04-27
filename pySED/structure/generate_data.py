'''
@author:
**************************  LiangTing ***************************
                      liangting.zj@gmail.com
************************ 2021/4/22 23:03:21 *********************
'''

import numpy as np
import sys
from pySED.structure import read_structure

class structure_maker(object):

    def __init__(self, structure_file_name='POSCAR',
                       lammps_data_flag=False,
                       lammps_infile_name=None,
                       description='Create LAMMPS or GPUMD data for SED method by using pySED Codes'):
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
            if self.lammpsInfile is None:
                print('*********** [ERROR] LAMMPS infile required! **********')
                sys.exit(1)
            structure = read_structure.get_structure_from_lammps(self.lammpsInfile, show_log=False)

        elif self.input_structure.endswith('.xyz'):
            structure = read_structure.read_from_file_structure_xyz(file_name=self.input_structure)

        else:
            structure = read_structure.read_from_file_structure_poscar(file_name=self.input_structure)

        self.unitcell_positions = structure['positions']
        self.cell = structure['cell']

        self.unitcell_num_atoms = self.unitcell_positions.shape[0]

        self.masses = structure['masses']
        self.basis_atoms_symbols = structure['symbols']

        self.unitcell_index = np.ones((self.unitcell_num_atoms, 1))
        self.basis_index = np.arange(1, self.unitcell_num_atoms + 1).reshape(self.unitcell_num_atoms, 1)

    def replicate_supercell(self, supercell=(1, 1, 1), origin_position=False):

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

    def write_xyz(self, filename='model.xyz', pbc="T T T"):
        lattice = self.calculate_supercell_matrix(self.cell, self.supercell)
        with open(filename, 'w') as fid:
            fid.write('{0:d}\n'.format(self.num_atoms))

            lattice_str = 'Lattice="{0:10.10f} {1:10.10f} {2:10.10f} {3:10.10f} {4:10.10f} {5:10.10f} {6:10.10f} ' \
                          '{7:10.10f} {8:10.10f}" Properties=species:S:1:pos:R:3 config_type="SED by pySED codes"'.format(
                           lattice[0, 0], lattice[0, 1], lattice[0, 2],
                           lattice[1, 0], lattice[1, 1], lattice[1, 2],
                           lattice[2, 0], lattice[2, 1], lattice[2, 2])

            fid.write(lattice_str + f' pbc="{pbc}"' + '\n')
            for i, row in enumerate(self.supercell_position):
                fid.write('{0:2} {1:20.10f} {2:20.10f} {3:20.10f}\n'.format(self.basis_atoms_symbols[i], row[0], row[1], row[2]))

        print(f'\n************* {filename} is written successfully ************\n')

    def write_lammps_data(self, file_name='lammps_data', lammps_data_types=None, atom_order=None):

        # Get lammps data boundary
        restricted_cell, rotation = self.calculate_restricted_cell(self.calculate_supercell_matrix(self.cell, self.supercell),
                                                                   return_rotation=True)
        xhi, yhi, zhi, xy, xz, yz = self.lammps_parameters_from_restricted_cell(restricted_cell)
        xlo, ylo, zlo = 0.0, 0.0, 0.0
        lammps_positions = self.supercell_position @ rotation

        # For atom types
        if atom_order is not None:
            # Convert basis_atoms_symbols to a plain set of strings
            missing = set(map(str, self.basis_atoms_symbols)) - set(map(str, atom_order))
            if missing:
                missing_str = ", ".join(missing)
                raise ValueError(
                    f"\n**************** The 'atom_order' list is missing the following atom types: {missing_str}. ****************\n"
                    f"************* Please check your 'atom_order' and make sure all species are included. ************\n"
                )
            self.unique_types = list(atom_order)

        else:
            self.unique_types = list(np.unique(self.basis_atoms_symbols))

        self.masses_for_lammps_output = [read_structure.symbol_to_mass(i) for i in self.unique_types]

        self.num_types = len(self.unique_types)

        # Get lammps unit for generate lammps data
        if lammps_data_types is None and self.lammps_data_flag:
            self.unit = read_structure.get_units(self.input_file)
            print('Lammps data units is ' + self.unit)

        elif lammps_data_types is not None:
            self.unit = lammps_data_types
            print('\n************* Lammps data is for {0} units *************'.format(self.unit))

        else:
            self.unit = 'real'
            print('\n*********** Lammps data units is {0} (using the default unit) *************'.format(self.unit))

        # write to file
        with open(file_name, 'w') as fid:
            fid.write(self.description + '\n\n')
            fid.write('{0} atoms\n\n'.format(self.num_atoms))
            fid.write('{0} atom types\n\n'.format(self.num_types))
            # Boundary
            fid.write(f'{xlo:20.10f} {xhi:20.10f} xlo xhi\n')
            fid.write(f'{ylo:20.10f} {yhi:20.10f} ylo yhi\n')
            fid.write(f'{zlo:20.10f} {zhi:20.10f} zlo zhi\n')
            fid.write(f'{xy:20.10f} {xz:20.10f} {yz:20.10f} xy xz yz\n\n')
            # Masses
            fid.write('Masses\n\n')
            for i in range(self.num_types):
                fid.write('{0} {1:20.10f}  # {2}\n'.format(i + 1, self.masses_for_lammps_output[i], self.unique_types[i]))

            fid.write('\nAtoms\n\n')

            if self.unit == 'metal':
                for i, row in enumerate(lammps_positions):
                    fid.write('{0} {1:d} {2:20.10f} {3:20.10f} {4:20.10f}  # {5}\n'.format(i + 1,
                                                                     self.unique_types.index(self.basis_atoms_symbols[i])+1,
                                                                     row[0], row[1], row[2],
                                                                     self.basis_atoms_symbols[i]))

            elif self.unit == 'real':
                # Counters for atom and molecule IDs
                atom_id = 0
                molecule_id = 0
                for i, row in enumerate(lammps_positions):
                    atom_id += 1
                    if i > 0 and i % self.unitcell_num_atoms == 0: # assign by supercell of z direction
                        molecule_id += 1
                    fid.write('{0} {1} {2:d} {3} {4:20.10f} {5:20.10f} {6:20.10f}  # {7}\n'.format(i + 1, molecule_id + 1, # ids and molecular's ids
                                                                     self.unique_types.index(self.basis_atoms_symbols[i]) + 1,  # atom types
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

    @staticmethod
    def calculate_supercell_matrix(cell, supercell):
        cell = np.asarray(cell, dtype=float)
        supercell = np.asarray(supercell, dtype=float)

        if cell.shape != (3, 3):
            raise ValueError('cell must be a 3x3 matrix with lattice vectors as rows.')
        if supercell.shape != (3,):
            raise ValueError('supercell must contain three replication factors.')
        if np.any(supercell <= 0):
            raise ValueError('supercell replication factors must be positive.')

        return np.diag(supercell) @ cell

    @staticmethod
    def calculate_restricted_cell(cell, return_rotation=False):
        """
        Convert any right-handed triclinic cell to the lower-triangular
        restricted representation used by LAMMPS. The returned rotation maps
        row-vector Cartesian coordinates from the input cell frame to the
        restricted frame: new_positions = old_positions @ rotation.
        """
        cell = np.asarray(cell, dtype=float)
        if cell.shape != (3, 3):
            raise ValueError('cell must be a 3x3 matrix with lattice vectors as rows.')

        avec, bvec, cvec = cell
        ax = np.linalg.norm(avec)
        if np.isclose(ax, 0.0):
            raise ValueError('The first lattice vector has zero length.')

        ex = avec / ax
        bx = np.dot(bvec, ex)
        by_vec = bvec - bx * ex
        by = np.linalg.norm(by_vec)
        if np.isclose(by, 0.0):
            raise ValueError('The first two lattice vectors are collinear.')

        ey = by_vec / by
        ez = np.cross(ex, ey)
        ez_norm = np.linalg.norm(ez)
        if np.isclose(ez_norm, 0.0):
            raise ValueError('The lattice vectors do not form a valid 3D cell.')
        ez = ez / ez_norm

        cx = np.dot(cvec, ex)
        cy = np.dot(cvec, ey)
        cz = np.dot(cvec, ez)
        if cz <= 0.0 and np.isclose(cz, 0.0):
            raise ValueError('The third lattice vector is coplanar with the first two.')
        if cz < 0.0:
            raise ValueError('The cell must be right-handed for LAMMPS restricted triclinic output.')

        restricted_cell = np.array([[ax, 0.0, 0.0],
                                    [bx, by, 0.0],
                                    [cx, cy, cz]])
        restricted_cell[np.isclose(restricted_cell, 0.0, atol=1e-12)] = 0.0

        rotation = np.vstack((ex, ey, ez)).T
        if return_rotation:
            return restricted_cell, rotation
        return restricted_cell

    @staticmethod
    def lammps_parameters_from_restricted_cell(restricted_cell):
        restricted_cell = np.asarray(restricted_cell, dtype=float)
        if restricted_cell.shape != (3, 3):
            raise ValueError('restricted_cell must be a 3x3 matrix.')

        xhi = restricted_cell[0, 0]
        yhi = restricted_cell[1, 1]
        zhi = restricted_cell[2, 2]
        xy = restricted_cell[1, 0]
        xz = restricted_cell[2, 0]
        yz = restricted_cell[2, 1]

        return xhi, yhi, zhi, xy, xz, yz

    @staticmethod
    def calculate_lattice_parameters(cell, supercell):
        supercell_matrix = structure_maker.calculate_supercell_matrix(cell, supercell)
        restricted_cell = structure_maker.calculate_restricted_cell(supercell_matrix)
        xhi, yhi, zhi, xy, xz, yz = structure_maker.lammps_parameters_from_restricted_cell(restricted_cell)

        return xhi, yhi, zhi, xy, xz, yz

if __name__ == "__main__":

    file_name = 'POSCAR_MoS2'                              ## The in file for lammps
    si = structure_maker(structure_file_name=file_name)
    si.replicate_supercell(supercell=(8, 8, 20))
    si.write_xyz()
    si.write_lammps_data()
    si.write_lattice_basis_file()







