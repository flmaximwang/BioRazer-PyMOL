from pymol import cmd
from ..utils import aa_constants

def select_hydrophobic_res(gly=True):
    """
    DESCRIPTION

    Select hydrophobic residues in the current object.

    USAGE

    select_hydrophobic_residues

    EXAMPLES

    select_hydrophobic_residues

    SEE ALSO

    select_hydrophobic_atoms
    """
    # Select hydrophobic residues
    if gly:
        selection = ' or '.join(['resn %s' % aa for aa in aa_constants.HYDROPHOBIC_AAS])
    else:
        selection = ' or '.join(['resn %s' % aa for aa in aa_constants.HYDROPHOBIC_AAS if aa != 'GLY'])
    cmd.select('hydrophobic_residues', selection)

    return selection

def atom_to_macro(atom, model="*", segi="*"):
    return f"/{model}/{segi}/{atom.chain_id}/{atom.res_name}`{atom.res_id}/{atom.atom_name}"

def generate_macro(model, segi, chain, resn, resi, name):
    return f"/{model}/{segi}/{chain}/{resn}`{resi}/{name}"

def macro_splitter(macro: str):
    """
    Split a macro into its components.

    Parameters
    ----------
    macro : str
        The macro to be split, like /model/segi/chain/resn`resi/name.

    Returns
    -------
    dict
        A dict containing the model, segi, chain, resn, resi and name.
    """
    empty, model, segi, chain, resn_resi, name = macro.split("/")
    print(macro.split("/"))
    resn, resi = resn_resi.split("`")
    res = {"model": model, "segi": segi, "chain": chain, "resn": resn, "resi": resi, "name": name}
    return res

def macro_to_resn_resi(macro: str):
    """
    Get the resn and resi from a macro.

    Parameters
    ----------
    macro : str
        The macro to be split, like /model/segi/chain/resn`resi/name.

    Returns
    -------
    tuple
        A tuple containing the resn and resi.
    """
    if macro_splitter(macro)["resn"] in aa_constants.THREE_TO_ONE_MAP:
        return aa_constants.THREE_TO_ONE_MAP[macro_splitter(macro)["resn"]] + macro_splitter(macro)["resi"]
    else:
        return macro_splitter(macro)["resn"] + macro_splitter(macro)["resi"]

# 批量导出该模块中定义的所有函数
all_functions = [
    select_hydrophobic_res,
    atom_to_macro,
    generate_macro
]
__all__ = [func.__name__ for func in all_functions]

for func in all_functions:
    cmd.extend(func.__name__, func)