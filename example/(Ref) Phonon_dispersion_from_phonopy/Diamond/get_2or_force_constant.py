'''
@author:
***************************  LiangTing ************************
                      liangting.zj@gmail.com
************************ 2020/8/27 14:22:15 *******************
'''
import sys
import numpy as np
import get_structure
from arrange import get_correct_arrangement
from phonopy_link import get_phonon, obtain_phonon_dispersion_bands

from phonopy.file_IO import write_FORCE_CONSTANTS, write_force_constants_to_hdf5, write_FORCE_SETS

# define the force unit conversion factors to LAMMPS metal style (eV/Angstrom)
unit_factors = {'real': 4.336410389526464e-2,
                'metal': 1.0,
                'si': 624150636.3094}

############################################################################
#########   LAMMPS-Class Second Order Force Constant Calculator   ##########
############################################################################

class Force_2or_lammps(object):
    def __init__(self,
                 lammps_input_file,
                 supercell_matrix = np.identity(3, dtype = int),   # must be integers in lammps
                 primitive_matrix = np.identity(3, dtype = float),
                 displacement_distance = 0.01,
                 show_log = False,
                 show_progress = True,
                 symmetrize=True):
        """
        Main Force_2or_LAMMPS class.
        ***** Force_2or_LAMMPS uses the Lammps-Python interface to output the force,
        and then uses the phonopy-interface to calculate the second-order force constant to obtain phonon information. *******
        :param lammps_input (str): LAMMPS input file name or list of commands
        :param supercell_matrix (float):  3x3 matrix supercell
        :param primitive cell (float):  3x3 matrix primitive cell
        :param displacement_distance: displacement distance in Angstroms (The default is 0.01 angstrom)
        :param show_log (boolean): Set true to display lammps log info
        :param show_progress (boolean): Set true to display progress of calculation
        """

        # Check if input is file or list of commands
        if type(lammps_input_file) is str:
            # read from file name
            self._lammps_input_file = lammps_input_file
            self._lammps_commands_list = open(lammps_input_file).read().split('\n')

        else:
            sys.exit('\n\tInput File error: Only LAMMPS files are supported\n')

        self._structure = get_structure.get_structure_from_lammps(self._lammps_input_file)

        self._supercell_matrix = supercell_matrix
        self._primitive_matrix = primitive_matrix
        self._displacement_distance = displacement_distance
        self._show_log = show_log
        self._show_progress = show_progress
        self._symmetrize = symmetrize

        self._force_constants = None
        self._data_set = None

        self.units = self._get_units(self._lammps_commands_list)

        if self.units not in unit_factors.keys():
            print('Units style not supported, use: {}'.format(unit_factors.keys()))
            exit()


    def _get_units(self, commands_list):
        """
        Get the units label for LAMMPS "units" command from a list of LAMMPS input commands
        :param commands_list: list of LAMMPS input commands (strings)
        :return units: string containing the units
        """
        for line in commands_list:
                if line.startswith('units'):
                    return line.split()[1]
        return 'lj'


    def get_force_constants(self, include_data_set = False):
        """
        calculate the force constants with phonopy using lammps to calculate forces
        :param include_data_set (Bool): Whether to return data_set
        :return: ForceConstants type object containing force constants
        """

        if self._force_constants is None:
            phonon = get_phonon(self._structure,
                                setup_forces = False,
                                super_cell_phonon = self._supercell_matrix,
                                primitive_matrix = self._primitive_matrix,
                                symmetrize = self._symmetrize)

            '''
            is_plusminus : 'auto', True, or False, optional
            For each atom, displacement of one direction (False), both
            direction, i.e., one directiona and its opposite direction
            (True), and both direction if symmetry requires ('auto').
            Default is 'auto'.
            '''

            phonon.generate_displacements(distance = self._displacement_distance, is_plusminus = 'auto')
            cells_with_disp = phonon.get_supercells_with_displacements()

            #cells_with_disp_1 = [cell.get_positions() for cell in cells_with_disp]

            '''
             """Displacement dataset for fc2
            dict
            Displacements in supercells. There are two types of formats.
            Type 1. Two atomic displacement in each supercell:
                {'natom': number of atoms in supercell,
                 'first_atoms': [
                   {'number': atom index of first displaced atom,
                    'displacement': displacement in Cartesian coordinates,
                    'forces': forces on atoms in supercell} ... ]}
            '''

            data_set = phonon.get_displacement_dataset()
            
            # Check forces for non displaced supercell (new in phonolammps)
            forces_supercell = self.get_forces(phonon.get_supercell())
            if np.max(forces_supercell) > 1e-1:
                warnings.warn('Large atomic forces found for non displaced structure: '
                              '{}. Make sure your unit cell is properly optimized'.format(np.max(forces_supercell)))
            
            # Get forces from lammps

            print('************ Starting to generate atomic displacements ************')

            for i, cell in enumerate(cells_with_disp):

                if self._show_progress:
                    print('              displacement atoms {} / {}'.format(i + 1, len(cells_with_disp)))

                forces = self.get_forces(cell)
                data_set['first_atoms'][i]['forces'] = forces

            print('************************ Done !!! ***********************')

            phonon.set_displacement_dataset(data_set)
            phonon.produce_force_constants()

            self._force_constants = phonon.get_force_constants()
            self._data_set = data_set

        if include_data_set:
            return self._force_constants, self._data_set
        else:
            return self._force_constants

    def get_forces(self, cell_with_disp = None):
        """
        Calculate the forces of a supercell using lammps
        :param cell_with_disp: supercell from which determine the forces
        :return: numpy array matrix with forces of atoms [Natoms x 3]
        """
        from lammps import lammps

        supercell_sizes = np.diag(self._supercell_matrix)

        #print(supercell_sizes)

        cmd_list = ['-log', 'generate_force.log']
        if not self._show_log:
            cmd_list += ['-echo', 'none', '-screen', 'none']

        lmp = lammps(cmdargs = cmd_list)
        lmp.commands_list(self._lammps_commands_list)
        lmp.command('replicate {} {} {}'.format(*supercell_sizes))
        lmp.command('run 0')

        natoms = lmp.get_natoms()

        coords_data = lmp.gather_atoms("x", 1, 3)

        reference = np.array(coords_data[:], dtype=np.dtype('double')).reshape((natoms, 3))
        template = get_correct_arrangement(reference, self._structure, self._supercell_matrix)
        indexing = np.argsort(template)

        coordinates = cell_with_disp.get_positions()

        # print(coordinates)
        for i in range(natoms):
            lmp.command('set atom {} x {} y {} z {}'.format(i + 1,
                                                            coordinates[template[i], 0],
                                                            coordinates[template[i], 1],
                                                            coordinates[template[i], 2]))
            lmp.command('print "   "')

        lmp.command('print "displacement atoms successfully"')

        lmp.command('run 0')

        forces = lmp.gather_atoms("f", 1, 3)

        forces = np.array([forces[i] for i in range(natoms * 3)]).reshape((natoms, 3))[indexing, :]
        forces = forces * unit_factors[self.units]

        lmp.close()

        return forces

    def write_force_constants(self, filename = 'FORCE_CONSTANTS', hdf5 = False):
        """
        Write the force constants in a file in phonopy plain text format
        :param filename: Force constants filename
        """

        force_constants = self.get_force_constants()
        if hdf5:
            write_force_constants_to_hdf5(force_constants,
                                          filename = 'fc2.hdf5', 
                                          p2s_map = None,
                                          physical_unit = None,
                                          compression = "gzip")
        else:
            write_FORCE_CONSTANTS(force_constants, filename = filename)

        print('\n******** {} file is written successfully **********'.format(filename))
        print('**** One can use it to generate Phonon dispersion with PHONOPY command line **** \n')

    def my_write_force_sets(self, filename='FORCE_SETS'):
        """
        Write the force sets in a file in phonopy plain text format
        :param filename: Force sets filename
        """

        _, data_set = self.get_force_constants(include_data_set=True)

        write_FORCE_SETS(data_set, filename = filename)

        print('\n*********** FORCE_SETS file is written successfully **********')

    def plot_phonon_dispersion_bands(self, use_seek_path = False,
                                     band_path = None,
                                     band_label = None,
                                     band_resolution = 51,
                                     save_band_structure_file = True,
                                     save_gv_file = True,
                                     save_ave_gv_file = True,
                                     write_band_yaml = True,
                                     show_dispersion = False):
        """
        Plot phonon band structure using seekpath automatic k-path
        Warning: The labels may be wrong if the structure is not standarized
        """
        
        import matplotlib.pyplot as plt

        def replace_list(text_string):
            substitutions = {'GAMMA': u'$\Gamma$'}

            for item in substitutions.items():
                text_string = text_string.replace(item[0], item[1])
   
            return text_string

        force_constants = self.get_force_constants()

        from seekpath_link import get_path_using_seek_path

        if use_seek_path:
            bands_and_labels = get_path_using_seek_path(phonopy_structure = self._structure)

        elif band_path is not None:
            bands_ranges = []
            for i in range(len(band_path)):
                if i < int(len(band_path))-1:
                    bands_ranges.append(np.array([band_path[i], band_path[i+1]]))
                else:
                    break

            bands_and_labels = {'ranges': bands_ranges,
                                'labels': band_label
                                }

        else:
            print("Please provide high symmetry points paths")

        print("****** Phonon dispersion relation data is being generated, please wait ******")

        _bands = obtain_phonon_dispersion_bands(self._structure,
                                                bands_and_labels['ranges'],
                                                force_constants,
                                                self._supercell_matrix,
                                                primitive_matrix = self._primitive_matrix,
                                                band_resolution = band_resolution,
                                                write_band_yaml = write_band_yaml)

        band_distances = []
        frequency = []
        group_velocities = []
        eigenvectors = []

        for i, freq in enumerate(_bands[1]):
            plt.plot(_bands[1][i], _bands[2][i], color='r')  # bands[1] = q_distance, _band[2] = frequency, _bands[3] = eigenvectors, _band[4] = group_velocities
            
            band_distances.append(_bands[1][i])
            frequency.append(_bands[2][i])
            eigenvectors.append(_bands[3][i])
            group_velocities.append(_bands[4][i])

        if save_band_structure_file:
            band_distances = np.reshape(band_distances, (np.size(band_distances, 0) * np.size(band_distances, 1), -1))
            frequency_1 = np.reshape(frequency, (np.size(frequency, 0) * np.size(frequency, 1), -1))

            np.savetxt('band_structure.txt', np.column_stack((band_distances, frequency_1)), delimiter='   ')
            print("********* Save band-structure file successfully *********")

        # ************* if save_ave_gv_file = True, save_gv_file must be True **************
        if save_ave_gv_file:
            save_gv_file = True

        if save_gv_file:
            frequency_2 = np.reshape(frequency, (-1, 1))
            group_velocities = np.reshape(group_velocities, (-1, 3))

            total_group_velocities= []
            for i in range(np.size(group_velocities, 0)):
                group_velocity = np.sqrt((((group_velocities[i][0]) * 100) ** 2 + ((group_velocities[i][1]) * 100) ** 2 + ((group_velocities[i][2]) * 100) ** 2))
                total_group_velocities.append(group_velocity)

            # Group velocity is in the unit of m/s
            np.savetxt('group_velocities.txt', np.column_stack((frequency_2, total_group_velocities)), delimiter='   ')
            print("********* Save group-velocities file successfully *********")

            if save_ave_gv_file:
                # ************** Use the API of more_itertools **************
                from more_itertools import chunked

                sort_freq = np.sort(frequency_2, axis=0)
                Gv_sort_index = np.argsort(frequency_2, axis=0)
                sort_group_velocities = np.array(total_group_velocities)[Gv_sort_index]

                ave_interval = np.size(frequency, 0) * np.size(frequency, 2)

                # Take an average every 'ave_interval' frequencies
                ave_group_velocities = [sum(x) / len(x) for x in chunked(sort_group_velocities, ave_interval)]
                ave_frequency = [sum(x) / len(x) for x in chunked(sort_freq, ave_interval)]

                # Group velocity is in the unit of m/s
                np.savetxt('group_ave_velocities.txt', np.column_stack((ave_frequency, ave_group_velocities)),
                           delimiter='   ')

                print("********* Save average group-velocities(frequency-dependent) file successfully *********")

        #plt.axes().get_xaxis().set_ticks([])
        plt.ylabel('Frequency (THz)', fontsize='x-large')
        plt.xlabel('Wave vector', fontsize='x-large')
        plt.xlim([0, _bands[1][-1][-1]])
        
        plt.axhline(y=0, color='k', ls='dashed')
        plt.suptitle('Phonon dispersion', fontsize='x-large')

        if 'labels' in bands_and_labels and bands_and_labels['labels'] is not None:

            # plt.rcParams.update({'mathtext.default': 'regular'})
            labels = bands_and_labels['labels']

            labels_e = []
            x_labels = []

            for i, freq in enumerate(_bands[1]):
                if labels[i][0] == labels[i - 1][1]:
                    labels_e.append(replace_list(labels[i][0]))
                else:
                    labels_e.append(
                        replace_list(labels[i - 1][1]) + '/' + replace_list(labels[i][0]))
                x_labels.append(_bands[1][i][0])
                
            x_labels.append(_bands[1][-1][-1])
            labels_e.append(replace_list(labels[-1][1]))
            labels_e[0] = replace_list(labels[0][0])
            
            ## ********** Output qx_axis_data for plot figures ***********
            np.savetxt('qx_axis_data.txt', np.column_stack(x_labels), delimiter='\n')

            plt.xticks(x_labels, labels_e, rotation='horizontal', fontsize='x-large')

        plt.savefig('Phonon_dispersion.png', format='png', dpi=1000, bbox_inches='tight')

        if show_dispersion:
            plt.show()


if __name__ == "__main__":
    in_file = 'in.lammps'  ## The in file for lammps

    force_constants = Force_2or_lammps(in_file,
                                       supercell_matrix=np.diag([4, 4, 4]),
                                       primitive_matrix=[[0.0, 0.5, 0.5],
                                                         [0.5, 0.0, 0.5],
                                                         [0.5, 0.5, 0.0]])

    # force_constants.my_write_force_sets()
    force_constants.write_force_constants()
    band_path = [[0.0, 0.0, 0.0],
                 [0.5, 0.0, 0.5]]   # Not need

    # band_label = [('GAMMA', 'X')]  #, ('X', 'S'), ('S', 'Y'), ('Y', 'GAMMA')]

    force_constants.plot_phonon_dispersion_bands(use_seek_path=True,
                                                 band_path=band_path,
                                                 # band_label = band_label,   # Auto determined
                                                 band_resolution=101)

    print('\n******************* COMPUTE ALL DONE !!!********************\n')
