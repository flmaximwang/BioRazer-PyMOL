from pymol import cmd
import numpy as np

def orient_vec_to(selection, ori_vec: np.ndarray, axis: str):
    '''
    Orient a selection with a transformation matrix,
    with which the orientation vector is aligned to the given axis
    '''
    axis_dict = {
        'x': np.array([1, 0, 0]),
        'y': np.array([0, 1, 0]),
        'z': np.array([0, 0, 1]),
    }
    if axis not in axis_dict:
        raise ValueError(f'Unknown axis: {axis}, available axes: {list(axis_dict.keys())}')
    ori_vec = np.array(ori_vec)
    ori_vec /= np.linalg.norm(ori_vec)
    axis_vec = axis_dict[axis]
    
    # Calculate the rotation matrix
    v = np.cross(ori_vec, axis_vec) # Calculate the rotation axis
    s = np.linalg.norm(v)
    c = np.dot(ori_vec, axis_vec)
    v_matrix = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0],
    ])
    R = np.eye(3) + v_matrix + np.dot(v_matrix, v_matrix) * (1 - c) / (s ** 2)
    
    