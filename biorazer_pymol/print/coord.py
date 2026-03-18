from pymol import cmd


def print_coord(selection: str):
    """
    Print the coordinates of atoms in the given selection.

    Args:
        selection (str): The selection string to specify which atoms to print.
    """
    model = cmd.get_model(selection)
    for atom in model.atom:
        print(
            f"{atom.name}: {atom.coord[0]:.3f}, {atom.coord[1]:.3f}, {atom.coord[2]:.3f}"
        )


cmd.extend("print_coord", print_coord)
