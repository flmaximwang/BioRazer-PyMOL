from pymol import cmd
import pandas as pd
from . import selectors
from . import visualizers
from ..utils import aa_constants

def visualize_interaction_matrix(interaction_matrix_path):
    '''
    Visualize the interaction matrix.
    The matrix is like:
    | atom1                          | atom2                          | distance (Å)                                       | Comments |
    | ------------------------------ | ------------------------------ | -------------------------------------------------- | -------- |
    | /8e5f_formatted//A/CYS`88/CB   | /8e5f_formatted//B/THR`54/OG1  | 3.3093292713165283                                 |          |
    | /8e5f_formatted//A/HIS`89/ND1  | /8e5f_formatted//B/LEU`53/O    | 3.0487306118011475                                 |          |
    
    Only columns "atom1", "atom2" and "distance (Å)" are required.
    '''
    
    interaction_matrix = pd.read_csv(interaction_matrix_path)
    for i in range(len(interaction_matrix)):
        atom1 = interaction_matrix.loc[i, 'atom1']
        atom2 = interaction_matrix.loc[i, 'atom2']
        atom1_abbr = selectors.macro_to_resn_resi(atom1)
        atom2_abbr = selectors.macro_to_resn_resi(atom2)
        cmd.distance(f'{atom1_abbr}_{atom2_abbr}', atom1, atom2)
        atom_1_main_flag = selectors.macro_splitter(atom1)["name"] in ["N", "O", "C"]
        atom_2_main_flag = selectors.macro_splitter(atom2)["name"] in ["N", "O", "C"]
        if atom_1_main_flag and atom_2_main_flag:
            visualizers.show_res_sticks(f"{atom1} or {atom2}")
        elif atom_1_main_flag:
            visualizers.show_res_sticks(atom1)
            visualizers.show_side_chain_sticks(atom2)
        elif atom_2_main_flag:
            visualizers.show_side_chain_sticks(atom1)
            visualizers.show_res_sticks(atom2)
        else:
            visualizers.show_side_chain_sticks(f"{atom1} or {atom2}")

all_funcs = [visualize_interaction_matrix]
__all__ = [func.__name__ for func in all_funcs]

for func in all_funcs:
    cmd.extend(func.__name__, func) 