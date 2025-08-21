from pymol import cmd, util


def color_by_residue_property():
    """
    Color the residues in the current selection by their properties.
    """
    util.cbac("resn ALA+ILE+LEU+VAL+MET+CYS+PRO")
    util.cbaw("resn TRP+TYR+PHE")
    util.cbag("resn GLU+ASP")
    util.cbas("resn HIS+ARG+LYS")
    util.cbac("resn SER+THR+GLY+ASN+GLN")


cmd.extend("color_by_residue_property", color_by_residue_property)
