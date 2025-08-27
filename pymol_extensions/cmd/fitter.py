from pymol import cmd
from .visualizers import *
from ..utils.coord_fitter import *

ORI_FITTER = svd_fit

def set_orientation_fitter(fitter):
    '''
    Set the orientation fitter
    fitter: the name fot fitters, ['svd', ...]
    '''
    switcher = {
        'svd': svd_fit,
    }
    if fitter not in switcher:
        raise ValueError(f'Unknown fitter: {fitter}, available fitters: {list(switcher.keys())}')
    global ORI_FITTER
    ORI_FITTER = switcher[fitter]

def fit_orientation(selection, symmetric=True, quiet=False, reverse=False, **kwargs):
    '''
    Fit the orientation of a set of coordinates and visualize it
    selection: selection string
    kwargs: additional arguments for visualize_orientation
    '''
    kwargs.pop('_self')
    coord_list = cmd.get_coords(selection)
    nm_vec, max_proj, min_proj = ORI_FITTER(coord_list)
    if reverse:
        nm_vec = -nm_vec
    visualize_orientation(
        nm_vec, 
        center=np.mean(coord_list, axis=0),
        scale = max_proj - min_proj,
        symmetric=symmetric,
        **kwargs
    )
    if not quiet:
        print('Fitted orientation:', nm_vec)
    return nm_vec

cmd.extend('fit_orientation', fit_orientation)