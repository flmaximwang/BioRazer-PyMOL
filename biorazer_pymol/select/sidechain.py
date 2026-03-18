from pymol import cmd


def select_sidechain(selection, selection_name="sele"):
    """
    Select sidechain atoms from a given selection.

    Parameters:
    selection (str): The selection string to filter sidechain atoms from.
    selection_name (str): The name for the new selection. Default is "sele".

    Returns:
    None
    """

    # Define the selection string for sidechain atoms
    sidechain_sel_str = f"({selection}) and not name N+C+CA+O"
    # Create the new selection in PyMOL
    cmd.select(selection_name, sidechain_sel_str)


cmd.extend("select_sidechain", select_sidechain)


def select_sidechain_plus_ca(selection, selection_name="sele"):
    """
    Select sidechain atoms excluding CA from a given selection.

    Parameters:
    selection (str): The selection string to filter sidechain atoms from.
    selection_name (str): The name for the new selection. Default is "sele".

    Returns:
    None
    """

    # Define the selection string for sidechain atoms excluding CA
    sidechain_no_ca_sel_str = f"({selection}) and not name N+C+O"
    # Create the new selection in PyMOL
    cmd.select(selection_name, sidechain_no_ca_sel_str)


cmd.extend("select_sidechain_plus_ca", select_sidechain_plus_ca)


def select_backbone(selection, selection_name="sele"):
    """
    Select backbone atoms from a given selection.

    Parameters:
    selection (str): The selection string to filter backbone atoms from.
    selection_name (str): The name for the new selection. Default is "sele".

    Returns:
    None
    """

    # Define the selection string for backbone atoms
    backbone_sel_str = f"({selection}) and name N+C+CA+O"
    # Create the new selection in PyMOL
    cmd.select(selection_name, backbone_sel_str)


cmd.extend("select_backbone", select_backbone)


def select_backbone_no_ca(selection, selection_name="sele"):
    """
    Select backbone atoms excluding CA from a given selection.

    Parameters:
    selection (str): The selection string to filter backbone atoms from.
    selection_name (str): The name for the new selection. Default is "sele".

    Returns:
    None
    """

    # Define the selection string for backbone atoms excluding CA
    backbone_no_ca_sel_str = f"({selection}) and name N+C+O"
    # Create the new selection in PyMOL
    cmd.select(selection_name, backbone_no_ca_sel_str)


cmd.extend("select_backbone_no_ca", select_backbone_no_ca)
