from pymol import cmd
import json
from pymol import stored
import biotite.structure as struc
import ammolite
import numpy as np
import pandas as pd
from . import analyzers_advanced
from . import selectors
from ..utils.name_utils import get_unused_name
from ..utils import macro
from ..utils import aa_constants


def hide_main_chain_sticks(selection):
    cmd.hide("sticks", f"name C+N+O and {selection}")


def show_main_chain_sticks(selection):
    cmd.show("sticks", f"name C+N+O and {selection}")


def show_side_chain_sticks(selection):
    cmd.show("sticks", f"not name C+N+O and byres {selection}")


def show_res_sticks(selection):
    cmd.show("sticks", f"byres {selection}")


def visualize_orientation(
    direction, center=[0, 0, 0], scale=1.0, symmetric=False, color="green", color2="red"
):
    """
    Draw an arrow. Helper function for "helix_orientation" etc.
    """
    from pymol import cgo

    color_list = cmd.get_color_tuple(color)
    color2_list = cmd.get_color_tuple(color2)
    if symmetric:
        scale *= 0.5
    end = np.array(center) + np.array(direction) * scale
    radius = 0.3
    obj = [cgo.SAUSAGE]
    obj.extend(center)
    obj.extend(end)
    obj.extend(
        [
            radius,
            0.8,
            0.8,
            0.8,
        ]
    )
    obj.extend(color_list)
    if symmetric:
        start = np.array(center) - np.array(direction) * scale
        obj.append(cgo.SAUSAGE)
        obj.extend(center)
        obj.extend(start)
        obj.extend(
            [
                radius,
                0.8,
                0.8,
                0.8,
            ]
        )
        obj.extend(color2_list)
    coneend = np.array(end) + np.array(direction) * 4.0 * radius
    if cmd.get_version()[1] >= 1.2:
        obj.append(cgo.CONE)
        obj.extend(end)
        obj.extend(coneend)
        obj.extend(
            [
                radius * 1.75,
                0.0,
            ]
        )
        obj.extend(color_list * 2)
        obj.extend(
            [
                1.0,
                1.0,  # Caps
            ]
        )
    cmd.load_cgo(obj, get_unused_name("oriVec"), zoom=0)


def display_cytc_multimer(
    selection="all",
    monomer_num=1,
    color1="green",
    color2="smudge",
    color3="salmon",
    color4="brown",
):
    # Color the monomer of the cytochrome c multimer
    if not isinstance(monomer_num, int):
        try:
            monomer_num = int(monomer_num)
        except:
            raise ValueError(
                f"monomer_num should be an integer or something that can be converted to an integer. While the input is: {monomer_num}"
            )
    for i in range(monomer_num):
        chain_id = chr(ord("A") + i)
        if i % 2 == 0:
            c_color = color1
            hec_color = color3
        else:
            c_color = color2
            hec_color = color4
        cmd.color(c_color, f"({selection}) and chain {chain_id}")
        cmd.color(hec_color, f"({selection}) and chain {chain_id} and resn HEC")
    color_NOS(selection)

    # 显示所有 Fe 附近 his 的 侧链
    cmd.show("sticks", "(byres(resn HIS and (name Fe around 3.5))) and not name C+O+N")
    cmd.color(
        "lightblue", "(byres(resn HIS and (name Fe around 3.5))) and not name C+O+N"
    )
    color_NOS(selection)

    # 显示所有 Heme 与 cys 的共价键
    cmd.show("sticks", "(byres(resn CYS and (resn HEC around 3.5))) and not name C+O+N")
    cmd.color(
        "gray70", "(byres(resn CYS and (resn HEC around 3.5))) and not name C+O+N"
    )
    color_NOS(selection)


def display_disulfide_bonds(selection="all", color="gray70", quiet=False, output=None):
    """
    Display disulfide bonds in the selection with the specified color.
    """
    sulfur_model = cmd.get_model(f"resn cys and name SG and {selection}")
    if len(sulfur_model.get_coord_list()) < 2:
        print("Not enough Cys residues to form disulfide bonds.")
        return
    cys_S_array = ammolite.convert_to_atom_array(sulfur_model)
    S_distance_matrix = analyzers_advanced.generate_distance_matrix(
        cys_S_array, cys_S_array
    )
    # print(list(map(lambda x: x.res_id, cys_S_array)))
    # print(S_distance_matrix)
    disulfide_bond_list = []
    # 只对上三角 (不含对角线) 迭代
    for i in range(0, S_distance_matrix.shape[0]):
        for j in range(i + 1, S_distance_matrix.shape[1]):
            if S_distance_matrix[i, j] < 2.5:
                disulfide_bond_list.append((cys_S_array[i], cys_S_array[j]))
    for bond in disulfide_bond_list:
        atom1, atom2 = bond
        cmd.show(
            "sticks",
            f"byres {selectors.atom_to_macro(atom1)} or byres {selectors.atom_to_macro(atom2)}",
        )
        cmd.color(
            color,
            f"((byres {selectors.atom_to_macro(atom1)}) and name SG) or ((byres {selectors.atom_to_macro(atom2)}) and name SG)",
        )

    if not quiet:
        for i, bond in enumerate(disulfide_bond_list):
            atom1, atom2 = bond
            print(
                f"Disulfide {i}: {selectors.atom_to_macro(atom1)} || {selectors.atom_to_macro(atom2)}"
            )
    if output is not None:
        with open(output, "a") as f:
            f.write(
                f"Disulfide {i}: {selectors.atom_to_macro(atom1)} || {selectors.atom_to_macro(atom2)}\n"
            )


def display_model_terminals(
    model, chains=[], N_term_color="blue", C_term_color="red", object_name="terminals"
):
    """
    分别在 model 的 N 端 CA 和 C 端 CA 创建一个球体，用于标记 N/C 端。
    """
    if len(chains) == 0:
        chains = cmd.get_chains(model)
    for chain in chains:
        # 获取 chain 的 N/C 端
        print(chain)
        res_indices = set(
            map(
                lambda x: int(x.resi),
                cmd.get_model(f"c. {chain} and n. CA and {model}").atom,
            )
        )
        N_term_resi = min(res_indices)
        C_term_resi = max(res_indices)
        N_term_CA = cmd.get_model(
            f"c. {chain} and n. CA and resi {N_term_resi} and {model}"
        ).atom[0]
        C_term_CA = cmd.get_model(
            f"c. {chain} and n. CA and resi {C_term_resi} and {model}"
        ).atom[0]

        # 创建 N/C 端的球体
        cmd.pseudoatom(
            object_name,
            pos=N_term_CA.coord,
            chain=chain,
            resi=N_term_resi,
            color=N_term_color,
            name="N",
            label="N terminus",
        )
        cmd.pseudoatom(
            object_name,
            pos=C_term_CA.coord,
            chain=chain,
            resi=C_term_resi,
            color=C_term_color,
            name="C",
            label="C terminus",
        )
    cmd.show("spheres", object_name)


def display_terminals(
    selection="all", chains=[], N_term_color="blue", C_term_color="red"
):
    """
    分别在 selection 的 N 端 CA 和 C 端 CA 创建一个球体，用于标记 N/C 端。
    """
    models = cmd.get_object_list(selection)
    for model in models:
        display_model_terminals(
            str(model),
            chains=chains,
            N_term_color=N_term_color,
            C_term_color=C_term_color,
            object_name=f"{model}_terminals",
        )


def load_contacts(contact_filepath, specificity_level=0):
    """
    The function `load_contacts` reads a CSV file containing contact information, filters the contacts
    based on a specificity level, and generates distance measurements in a molecular visualization
    software.

    :param contact_filepath: The `contact_filepath` parameter in the `load_contacts` function is the
    file path to the contacts data that you want to load and process. This function reads the contacts
    data from the specified file and performs certain operations based on the `specificity_level`
    parameter
    :param specificity_level: The `specificity_level` parameter in the `load_contacts` function is used
    to filter the contacts based on a specificity index. If a `specificity_level` is provided, only
    contacts with a specificity index greater than or equal to the specified level will be included in
    the analysis, defaults to 0 (optional)
    """

    my_contacts = pd.read_csv(contact_filepath)
    if isinstance(specificity_level, str):
        specificity_level = float(specificity_level)
    if specificity_level > 0:
        my_contacts = my_contacts.loc[
            my_contacts["Specificity Index"] >= specificity_level, :
        ]
    for i in my_contacts.index:
        atom1_dict = macro.read_macro(my_contacts.loc[i, "Atom1"])
        atom2_dict = macro.read_macro(my_contacts.loc[i, "Atom2"])
        atom1_identifier = f"{aa_constants.THREE_TO_ONE_MAP[atom1_dict['resn']]}{atom1_dict['resi']}_{atom1_dict['name']}"
        atom2_identifier = f"{aa_constants.THREE_TO_ONE_MAP[atom2_dict['resn']]}{atom2_dict['resi']}_{atom2_dict['name']}"
        measurement_name = f"{atom1_identifier}-{atom2_identifier}"
        cmd.distance(
            measurement_name, my_contacts.loc[i, "Atom1"], my_contacts.loc[i, "Atom2"]
        )
        cmd.show(
            "sticks",
            f"byres {my_contacts.loc[i, 'Atom1']} or byres {my_contacts.loc[i, 'Atom2']}",
        )


def load_interface_annotations(
    annotation_json,
    selection="all",
    object_name="interface_annotations",
    method="create",
):
    with open(annotation_json, "r") as f:
        interface_annotations = json.load(f)
    res = {}
    for chain_id, res_id in interface_annotations["chain_id, res_id"]:
        if chain_id not in res:
            res[chain_id] = []
        res[chain_id].append(str(res_id))
    selection_list = []
    for chain_id, res_ids in res.items():
        selection_list.append(f"(chain {chain_id} and resi {'+'.join(res_ids)})")
    selector = selection + " and " + f'({" or ".join(selection_list)})'
    print(selector)
    if method == "create":
        cmd.create(object_name, selector)
        cmd.hide("everything", object_name)
        cmd.show("sticks", object_name + " and not name C+O+N")
    else:
        cmd.select(object_name, selector)
        cmd.show("sticks", object_name + " and not name C+O+N")


all_funcs = [
    hide_main_chain_sticks,
    show_main_chain_sticks,
    show_side_chain_sticks,
    show_res_sticks,
    visualize_orientation,
    display_cytc_multimer,
    display_disulfide_bonds,
    load_contacts,
    display_model_terminals,
    display_terminals,
    load_interface_annotations,
]
__all__ = [func.__name__ for func in all_funcs]

for func in all_funcs:
    cmd.extend(func.__name__, func)
