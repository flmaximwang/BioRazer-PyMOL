from pymol import cmd


def point(x, y, z, name="point", radius=0.1, color="white"):
    """
    Places a pseudoatom at the specified coordinates

    Parameters:
    x (float): X coordinate of the point
    y (float): Y coordinate of the point
    z (float): Z coordinate of the point
    name (str): Name of the pseudoatom object (default: "point")
    radius (float): Radius of the sphere representing the point (default: 0.1)
    color (str): Color of the sphere representing the point (default: "white")
    """

    cmd.delete(name)
    cmd.pseudoatom(name, pos=[x, y, z])
    cmd.show("spheres", name)
    cmd.set("sphere_scale", radius, name)
    cmd.color(color, name)


cmd.extend("point", point)
