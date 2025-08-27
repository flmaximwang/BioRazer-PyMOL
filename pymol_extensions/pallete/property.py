from pymol import cmd, util


def color_by_element(sele="all"):
    cmd.color("red", f"({sele}) and elem O")
    cmd.color("blue", f"({sele}) and elem N")
    cmd.color("yellow", f"({sele}) and elem S")


cmd.extend("color_element", color_by_element)


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


def color_by_plddt(selection="all"):
    """
    The function `color_AF_pLDDT` assigns colors to different selection ranges based on a specified
    threshold for a molecular visualization program.

    :param selection: The `color_AF_pLDDT` function you provided seems to be a PyMOL script for coloring
    protein structures based on the pLDDT score ranges. The function defines color ranges for blue,
    light blue, yellow, and orange, and then assigns colors to different pLDDT, defaults to all
    (optional)
    """

    print(">>> Blue = very high (pLDDT > 90)")
    print(">>> Light blue = confident (90 > pLDDT > 70)")
    print(">>> Yellow = low (70 > pLDDT > 50)")
    print(">>> Orange = very low (pLDDT < 50)")

    blue_rgb = [0, 76, 202]
    blue = []
    for c in blue_rgb:
        blue.append(c / 255.0)

    lightblue_rgb = [73, 196, 238]
    lightblue = []
    for c in lightblue_rgb:
        lightblue.append(c / 255.0)

    yellow_rgb = [255, 213, 57]
    yellow = []
    for c in yellow_rgb:
        yellow.append(c / 255.0)

    orange_rgb = [255, 113, 67]
    orange = []
    for c in orange_rgb:
        orange.append(c / 255.0)

    # select and colour blue
    blue_upper = 100.0
    blue_lower = 90.0

    blue_sel_str = (
        selection + " & ! b < " + str(blue_lower) + " & ! b > " + str(blue_upper)
    )
    cmd.select("very_high", blue_sel_str)
    cmd.set_color("blue_plddt", blue)
    cmd.color("blue_plddt", "very_high")

    # select and colour lightblue
    lightblue_upper = 90.0
    lightblue_lower = 70.0

    lightblue_sel_str = (
        selection
        + " & ! b < "
        + str(lightblue_lower)
        + " & ! b > "
        + str(lightblue_upper)
    )
    cmd.select("confident", lightblue_sel_str)
    cmd.set_color("lightblue_plddt", lightblue)
    cmd.color("lightblue_plddt", "confident")

    # select and colour yellow
    yellow_upper = 70.0
    yellow_lower = 50.0

    yellow_sel_str = (
        selection + " & ! b < " + str(yellow_lower) + " & ! b > " + str(yellow_upper)
    )
    cmd.select("low", yellow_sel_str)
    cmd.set_color("yellow_plddt", yellow)
    cmd.color("yellow_plddt", "low")

    # select and colour orange
    orange_upper = 50.0
    orange_lower = 0.0

    orange_sel_str = (
        selection + " & ! b < " + str(orange_lower) + " & ! b > " + str(orange_upper)
    )
    cmd.select("very_low", orange_sel_str)
    cmd.set_color("orange_plddt", orange)
    cmd.color("orange_plddt", "very_low")

    cmd.deselect()


cmd.extend("color_by_plddt", color_by_plddt)
