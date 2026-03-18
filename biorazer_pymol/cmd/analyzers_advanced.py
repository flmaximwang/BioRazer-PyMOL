from pymol import cmd
import os, re
from pymol import stored
import numpy as np
import pandas as pd
from . import selectors
from . import visualizers
from ..utils import aa_constants
import biotite.structure as struc

def find_interactions_between(selection1, selection2, cutoff=5, show_interactions=True, name="I", out_dir=".", log=False):
    '''
    find_interactions_between selection1, selection2, [cutoff=5, [show_interactions=True, [name="I", [out_dir="de]]]]
    Automatically find interactions between two selections. This functions doesn't work if you use default PyMOL, because pandas, ammolite and biotite are not included in the default PyMOL.

    '''
    

    sele1_interaction_name = "stored.sele1_interaction"
    sele2_interaction_name = "stored.sele2_interaction"
    exec(f"{sele1_interaction_name} = set()")
    exec(f"{sele2_interaction_name} = set()")
    cmd.iterate(f"byres ({selection1} and not elem H) within {cutoff} of ({selection2} and not elem H)", f"{sele1_interaction_name}.add((model, segi, chain, resn, resi))")
    cmd.iterate(f"byres ({selection2} and not elem H) within {cutoff} of ({selection1} and not elem H)", f"{sele2_interaction_name}.add((model, segi, chain, resn, resi))")
    
    # Because dicts are not hashable, we need to convert them to tuples, and then convert them back to dicts
    exec(f"{sele1_interaction_name} = list(map(lambda x: {{'model': x[0], 'segi': x[1], 'chain': x[2], 'resn': x[3], 'resi': x[4]}}, {sele1_interaction_name}))")
    exec(f"{sele2_interaction_name} = list(map(lambda x: {{'model': x[0], 'segi': x[1], 'chain': x[2], 'resn': x[3], 'resi': x[4]}}, {sele2_interaction_name}))")
    sele1_interactions = list(eval(sele1_interaction_name))
    sele2_interactions = list(eval(sele2_interaction_name))
    sele1_interactions.sort(key=lambda x: int(x['resi']))
    sele2_interactions.sort(key=lambda x: int(x['resi']))
    
    res_suffix = f"{name}_{selection1}|{selection2}".replace(" ", "_")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    contact_matrix = generate_contact_matrix(sele1_interactions, sele2_interactions, result_filepath=os.path.join(out_dir, f"contact_matrix {res_suffix}"))
    contacting_atom_list = generate_contacting_atom_list(contact_matrix, result_file_path=os.path.join(out_dir, f"contacting_atom_list {res_suffix}"))
    
    if show_interactions:
        main_chain_resi_set = set()
        pros = set()
        glys = set()
        for index, row in contacting_atom_list.iterrows():
            atom_pattern = re.compile(r"/([^/]*)/([^/]*)/([^/]*)/([^/]*)`([^/]*)/([^/]*)")
            atoms = []
            for i in ['atom1', 'atom2']:
                atoms.append(row[i])
            atom_labels = []
            for atom in atoms:
                atom_match = atom_pattern.match(atom)
                if atom_match.group(4) in aa_constants.THREE_TO_ONE_MAP:
                    atom_labels.append(aa_constants.THREE_TO_ONE_MAP[atom_match.group(4)] + atom_match.group(5))
                else:
                    atom_labels.append(atom_match.group(4) + atom_match.group(5))
                cmd.show("sticks", f"byres {atom}")
                visualizers.hide_main_chain_sticks(f"byres {atom}")
                with open("output.log", "a") as log_output:
                    log_output.write(atom_match.group(4) + "\n")
                if atom_match.group(6) in ["C", "O", "N"]:
                    main_chain_resi_set.add(selectors.generate_macro(atom_match.group(1), atom_match.group(2), atom_match.group(3), atom_match.group(4), atom_match.group(5), "*"))
                if atom_match.group(4) == "PRO":
                    pros.add(selectors.generate_macro(atom_match.group(1), atom_match.group(2), atom_match.group(3), atom_match.group(4), atom_match.group(5), "*"))
                if atom_match.group(4) == "GLY":
                    glys.add(selectors.generate_macro(atom_match.group(1), atom_match.group(2), atom_match.group(3), atom_match.group(4), atom_match.group(5), "*"))
            cmd.distance(f"{name}_{atom_labels[0]}-{atom_labels[1]}", atoms[0], atoms[1])
        for resi in main_chain_resi_set:
            visualizers.show_main_chain_sticks(resi)
        for pro in pros:
            cmd.show("sticks", f"(byres {pro}) and name N")
        for gly in glys:
            cmd.show("sticks", f"(byres {gly}) and name O")
        if log:
            for pro in pros:
                with open("output.log", "a") as log_output:
                    log_output.write(f"Showing main chain sticks for {pro}\n")

def generate_distance_matrix(list1, list2):
    '''
    Receive 2 lists of atoms and returns a numpy distance matrix.
    Several formats of lists are acceptable:
    - selection
    - atom array
    
    Hydrogen atoms are ignored.
    '''
    
    if not type(list1) == type(list2):
        raise ValueError("The two lists must be of the same type.")
    atomarrays = []
    if isinstance(list1, str):
        for selection in [list1, list2]:
            atomarrays.append(ammolite.convert_to_atom_array(cmd.get_model(selection)))
    elif isinstance(list1, struc.AtomArray):
        for atomarray in [list1, list2]:
            atomarrays.append(atomarray)
    # print(atomarrays[0])
    distance_matrix = np.zeros((len(atomarrays[0]), len(atomarrays[1])))
    for i in range(len(atomarrays[0])):
        for j in range(len(atomarrays[1])):
            distance_matrix[i, j] = struc.distance(atomarrays[0][i], atomarrays[1][j])
    return distance_matrix

def generate_contact_matrix(list1, list2, result_filepath = "contact_matrix"):
    '''
    Receive 2 lists of residues. Elements in each list is like {'chain': chain, 'resi': resi, 'resn': resn}.
    
    Return the contact matrix (pandas.Dataframe) between two lists of residues.
    A contacting matrix is like, while each element in the table is a distance in Å:
    |                                  | /model/segi/chain/resi`resn/name | /model/segi/chain/resi`resn/name |
    |----------------------------------|----------------------------------|----------------------------------|
    | /model/segi/chain/resi`resn/name | 0.0                              | 0.0                              |
    '''
    
    result = pd.DataFrame()
    for i in list1:
        for j in list2:
            shortest_length, atom1, atom2 = get_nearest_distance_between(
                selectors.generate_macro(i['model'], i['segi'], i['chain'], i['resn'], i['resi'], "*"),
                selectors.generate_macro(j['model'], j['segi'], j['chain'], j['resn'], j['resi'], "*")
            )
            result.loc[selectors.atom_to_macro(atom1, i['model'], i['segi']), selectors.atom_to_macro(atom2, j['model'], j['segi'])] = shortest_length
    result.to_csv(f"{result_filepath}.csv")
    return result

def generate_contacting_atom_list(contact_matrix, result_file_path = "contacting_atom_list"):
    '''
    Receive a contacting matrix
    A contacting matrix is like, while each element in the table is a distance in Å:
    |                     | ///chain/resi`resn/ | ///chain/resi`resn/ | ///chain/resi`resn/ | ///chain/resi`resn/ |
    |---------------------|---------------------|---------------------|---------------------|---------------------|
    | ///chain/resi`resn/ | 0.0                 | 0.0                 | 0.0                 | 0.0                 |
    
    Return a table of contacting atoms.
    | atom1                   | atom2                   | distance (Å)        |
    |-------------------------|-------------------------|---------------------|
    | ///chain/resi`resn/name | ///chain/resi`resn/name | 0.0                 |
    '''
    
    result = pd.DataFrame({
        "atom1": [],
        "atom2": [],
        "distance (Å)": []
    })
    for i in contact_matrix.index:
        for j in contact_matrix.columns:
            if contact_matrix.loc[i, j] <= 3.5:
                result = pd.concat([result,
                    pd.DataFrame({
                        "atom1": [i],
                        "atom2": [j],
                        "distance (Å)": [contact_matrix.loc[i, j]]
                    })
                ])
    result.to_csv(f"{result_file_path}.csv", index=False)
    return result

def get_nearest_distance_between(selection1, selection2, message=False, ignore_hydrogens=True):
    '''
    Receive 2 selections
    Return the nearest distance between two selections, as well as the corresponding 2 atoms
    '''
    if ignore_hydrogens:
        model1 = cmd.get_model(selection1 + " and not elem H")
        model2 = cmd.get_model(selection2 + " and not elem H")
    else:
        model1 = cmd.get_model(selection1)
        model2 = cmd.get_model(selection2)
    model1 = ammolite.convert_to_atom_array(model1)
    model2 = ammolite.convert_to_atom_array(model2)
    shortest_length = 10000
    atom1 = None
    atom2 = None
    for i in model1:
        for j in model2:
            temp_length = struc.distance(i, j)
            if temp_length < shortest_length:
                shortest_length = temp_length
                atom1 = i
                atom2 = j
    if message:
        print(f"Shortest length between {selection1} and {selection2} is", shortest_length, "Å")
    return shortest_length, atom1, atom2

all_funcs = [
    find_interactions_between,
    generate_contact_matrix,
    generate_contacting_atom_list,
    get_nearest_distance_between,
    generate_distance_matrix,
]
__all__ = [func.__name__ for func in all_funcs]

for func in all_funcs:
    cmd.extend(func.__name__, func)