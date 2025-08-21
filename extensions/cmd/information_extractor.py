from pymol import cmd

def extract_residue_indexes(selection):
    """
    Extracts the indexes of residues from a given selection in PyMOL.

    Args:
        selection (str): The selection string to extract residue indexes from.

    Returns:
        list: A list of residue indexes.
    """
    # Get the selection object
    sel = cmd.get_model(selection)
    
    # Extract residue indexes
    residue_indexes = [residue for residue in sel.get_residues()]
    
    return residue_indexes