from pymol import cmd


def select_hydrophobic(selection, selection_name="sele"):
    """
    Select hydrophobic atoms from a given selection.

    Parameters:
    selection (str): The selection string to filter hydrophobic atoms from.
    selection_name (str): The name for the new selection. Default is "sele".

    Returns:
    None
    """

    # Define the selection string for hydrophobic atoms
    hydrophobic_sel_str = f"({selection}) and (resn ALA+VAL+ILE+LEU+MET+PHE+TRP+PRO)"
    # Create the new selection in PyMOL
    cmd.select(selection_name, hydrophobic_sel_str)


def select_hydrophilic(selection, selection_name="sele"):
    """
    Select hydrophilic atoms from a given selection.

    Parameters:
    selection (str): The selection string to filter hydrophilic atoms from.
    selection_name (str): The name for the new selection. Default is "sele".

    Returns:
    None
    """

    # Define the selection string for hydrophilic atoms
    hydrophilic_sel_str = (
        f"({selection}) and (resn ARG+ASN+ASP+GLN+GLU+HIS+LYS+SER+THR+TYR+CYS+GLY)"
    )
    # Create the new selection in PyMOL
    cmd.select(selection_name, hydrophilic_sel_str)
