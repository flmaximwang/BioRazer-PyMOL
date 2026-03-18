from pymol import cmd
import numpy as np
from pymol import cgo


def _arrow_between(
    point_start: np.ndarray,
    point_end: np.ndarray,
    name: str,
    cylinder_radius: float,
    cone_ratio: float,
    rgb_color: np.ndarray,
):
    direction = point_end - point_start

    length = np.linalg.norm(direction)
    if length == 0:
        raise ValueError("Start and end points cannot be the same")

    direction /= length

    cmd.delete(name)
    cone_length = min(length * cone_ratio, 0.6)
    obj = [
        cgo.CYLINDER,
        point_start[0],
        point_start[1],
        point_start[2],
        point_end[0],
        point_end[1],
        point_end[2],
        cylinder_radius,
        rgb_color[0],
        rgb_color[1],
        rgb_color[2],
        rgb_color[0],
        rgb_color[1],
        rgb_color[2],
        cgo.CONE,
        point_end[0],
        point_end[1],
        point_end[2],
        point_end[0] + cone_length * direction[0],
        point_end[1] + cone_length * direction[1],
        point_end[2] + cone_length * direction[2],
        cylinder_radius * 4,
        0.0,
        rgb_color[0],
        rgb_color[1],
        rgb_color[2],
        rgb_color[0],
        rgb_color[1],
        rgb_color[2],
        1.0,  # cone cap
        1.0,  # cone base
    ]
    cmd.load_cgo(obj, name)


def arrow_between(
    start_x,
    start_y,
    start_z,
    end_x,
    end_y,
    end_z,
    name="arrow",
    cylinder_radius=0.1,
    cone_ratio=0.2,
    r_color=1.0,
    g_color=1.0,
    b_color=0.0,
):
    _arrow_between(
        point_start=np.array([float(start_x), float(start_y), float(start_z)]),
        point_end=np.array([float(end_x), float(end_y), float(end_z)]),
        name=name,
        cylinder_radius=float(cylinder_radius),
        cone_ratio=float(cone_ratio),
        rgb_color=np.array([float(r_color), float(g_color), float(b_color)]),
    )


def arrow_pass(
    point_x,
    point_y,
    point_z,
    vector_x,
    vector_y,
    vector_z,
    length=1.0,
    name="arrow",
    cylinder_radius=0.1,
    cone_ratio=0.2,
    r_color=1.0,
    g_color=1.0,
    b_color=0.0,
):
    """
    DESCRIPTION

    Draw an arrow that passes through a point and follows a direction vector.
    The arrow is centered at (point_x, point_y, point_z) and extends by
    +/- length along the normalized direction vector.

    USAGE

    arrow_pass point_x, point_y, point_z, vector_x, vector_y, vector_z
    [, length [, name [, cylinder_radius [, cone_ratio
    [, r_color [, g_color [, b_color ]]]]]]]

    ARGUMENTS

    point_x, point_y, point_z = float
        Coordinates of the point that the arrow passes through.

    vector_x, vector_y, vector_z = float
        Direction vector of the arrow (cannot be all zeros).

    length = float, default=1.0
        Half-length from center point to each side of the arrow.

    name = str, default=arrow
        Name of the CGO object to create.

    cylinder_radius = float, default=0.1
        Radius of the cylinder body.

    cone_ratio = float, default=0.2
        Cone length as a fraction of the full arrow length.

    r_color, g_color, b_color = float, default=1.0, 1.0, 0.0
        RGB color values in the range [0.0, 1.0].

    EXAMPLE

    arrow_pass 0, 0, 0, 1, 0, 0, 4, x_axis_arrow, 0.12, 0.2, 1, 0.2, 0.2
    """
    length = float(length)
    cylinder_radius = float(cylinder_radius)
    cone_ratio = float(cone_ratio)
    r_color = float(r_color)
    g_color = float(g_color)
    b_color = float(b_color)

    point_passed = np.array([float(point_x), float(point_y), float(point_z)])
    vector_passed = np.array([float(vector_x), float(vector_y), float(vector_z)])
    vector_norm = np.linalg.norm(vector_passed)
    if vector_norm == 0:
        raise ValueError("Direction vector cannot be zero")
    direction = vector_passed / vector_norm
    start_point = point_passed - direction * length
    end_point = point_passed + direction * length
    _arrow_between(
        point_start=start_point,
        point_end=end_point,
        name=name,
        cylinder_radius=cylinder_radius,
        cone_ratio=cone_ratio,
        rgb_color=np.array([r_color, g_color, b_color]),
    )


cmd.extend("arrow_between", arrow_between)
cmd.extend("arrow_pass", arrow_pass)
