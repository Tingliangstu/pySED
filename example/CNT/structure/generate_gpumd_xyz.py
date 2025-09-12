'''
Call pySED module to generate model.xyz for GPUMD run and the required basis.in file for pySED
@author:
**************************  LiangTing ***************************
                      liangting.zj@gmail.com
************************ 2022/4/22 23:03:21 *********************
'''

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
    
    # Write lammps data files for run (lammps_data_types can be 'metal' or 'real')
    #structure.write_lammps_data(lammps_data_types='metal')
    
    # Write basis.in files for further used
    structure.write_lattice_basis_file()

if __name__ == "__main__":

    file_name = 'POSCAR_CNT'                              ## The in file for lammps
    supercell = (1, 1, 160)
    generate_required_files(file_name, supercell)
