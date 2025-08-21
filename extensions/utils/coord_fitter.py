'''
This module collects multiple functions that are used to fit a set of coordinates
'''

import numpy as np

def svd_fit(coord_list):
    '''
    This function uses the SVD method to fit a set of coordinates
    Return:
    - the fitted orientation (a normalize vector),
    - the maximum projection of the coordinates onto the orientation
    - the minimum projection of the coordinates onto the orientation
    Params:
    - coord_list: list of coordinates like [[x1, y1, z1], [x2, y2, z2], ...]
    '''
    np_coord_list = np.array(coord_list) - np.array(coord_list).mean(0)
    U, s, Vh = np.linalg.svd(np_coord_list)
    nm_vec = Vh[0] / np.linalg.norm(Vh[0])
    max_proj = np.dot(np_coord_list, nm_vec).max()
    min_proj = np.dot(np_coord_list, nm_vec).min()
    return nm_vec, max_proj, min_proj