import numpy as np

def get_path_using_seek_path(phonopy_structure):
    """
    Obtain the path in reciprocal space to plot the phonon band structure
    :return: dictionary with list of q-points and labels of high symmetry points
    """
    try:
        import seekpath

        cell = phonopy_structure.get_cell()
        positions = phonopy_structure.get_scaled_positions()
        numbers = np.unique(phonopy_structure.get_chemical_symbols(), return_inverse=True)[1]

        path_data = seekpath.get_path((cell, positions, numbers))

        labels = path_data['point_coords']
        band_ranges = []

        for set in path_data['path']:
            band_ranges.append([labels[set[0]], labels[set[1]]])

        return {'ranges': band_ranges,
                'labels': path_data['path']}

    except ImportError:

        print('Seekpath not installed. Autopath is deactivated')
        band_ranges = ([[[0.0, 0.0, 0.0], [0.5, 0.0, 0.5]]])

        return {'ranges': band_ranges,
                'labels': [['GAMMA', '1/2 0 1/2']]}