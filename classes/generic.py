import dataclasses
from abc import abstractproperty
import enum
from functools import cached_property
import logging
import math
import sys

import gdspy
from matplotlib.path import Path
from matplotlib.textpath import TextPath


root = logging.getLogger()
root.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)

logger = logging.getLogger(__name__)


class Unit(enum.Enum):
    nm = 1
    um = 10**3
    mm = 10**6


def better_dataclass(cls):
    """  """
    cls = dataclasses.dataclass(cls)  # Apply the built-in dataclass decorator first

    cls.attribute_names = [field.name for field in dataclasses.fields(cls)]

    def values(self):
        return list(self.as_dict().values())
    cls.greet = property(values)

    def as_dict(self):
        return dataclasses.asdict(self)
    cls.as_dict = as_dict

    def __str__(self):
        return f"{type(cls).__name__}(" + "\n\t".join(f"{k}: {v}" for k, v in self.as_dict().items()) + ")"

    return cls  # Return the modified class


def config_dataclass(cls):
    """  """
    cls = better_dataclass(cls)  # Apply the built-in dataclass decorator first

    cls.list_type_attrib_names = [field.name for field in dataclasses.fields(cls) if field.type == list]

    def list_type_attrib_values(self):
        data = self.as_dict()
        return [data[name] for name in self.list_type_attrib_names]
    cls.list_type_attrib_values = property(list_type_attrib_values)

    return cls  # Return the modified class


class Component:
    
    """
    
    """
    
    god_cell = gdspy.Cell("god_cell")
    reference_cell = gdspy.Cell("reference_cell")

    @abstractproperty
    def polygon(self):
        """ asdasd. """

    @classmethod
    def add_components_to_god_cell(cls, components):
        for component in components:
            cls.god_cell.add(component.polygon)


class Rectangle(gdspy.Rectangle):

    """
    
    """
    
    def __init__(self, width, height, layer=0, datatype=0, center_to_origin=True):
        self.width = width
        self.height = height
        self.center_to_origin = center_to_origin
        super().__init__(self.top_left_coords, self.bot_right_coords, layer=layer, datatype=datatype)

    @cached_property
    def top_left_coords(self):
        return (-self.width/2, self.height/2) if self.center_to_origin else (0, self.height)
    
    @cached_property
    def bot_right_coords(self):
        return (self.width/2, -self.height/2) if self.center_to_origin else (self.width, 0)
    
    @cached_property
    def bounding_box(self):
        return self.get_bounding_box()

    @cached_property
    def min_x(self):
        return self.bounding_box[0][0]

    @cached_property
    def min_y(self):
        return self.bounding_box[0][1]

    @cached_property
    def max_x(self):
        return self.bounding_box[1][0]

    @cached_property
    def max_y(self):
        return self.bounding_box[1][1]


class Square(Rectangle):

    """
    
    """
    
    def __init__(self, side_length, layer=0, center_to_origin=True):
        super().__init__(side_length, side_length, layer, center_to_origin)


class PolygonOperations:
    
    @staticmethod
    def create_polygons_from_str(some_str, size, position, font_prop=None, tolerance=0.1):
        """ I'm leaving it like this. """
        polys = []
        xmax = position[0]
        c = None
        for points, code in TextPath(position, some_str, size=float(size), prop=font_prop).iter_segments():
            if code == Path.MOVETO:
                c = gdspy.Curve(*points, tolerance=tolerance)
            if c is None:
                continue
            if code == Path.LINETO:
                c.L(*points)
            elif code == Path.CURVE3:
                c.Q(*points)
            elif code == Path.CURVE4:
                c.C(*points)
            elif code == Path.CLOSEPOLY:
                if (poly := c.get_points()).size <= 0:
                    continue
                if poly[:, 0].min() < xmax:
                    for i in range(len(polys) - 1, -1, -1):
                        if gdspy.inside(poly[:1], [polys[i]], precision=0.1 * tolerance)[0]:
                            poly = gdspy.boolean([polys.pop(i)], [poly], "xor", precision=0.1 * tolerance, max_points=0).polygons[0]
                            break
                        elif gdspy.inside(polys[i][:1], [poly], precision=0.1 * tolerance)[0]:
                            poly = gdspy.boolean([polys.pop(i)], [poly], "xor", precision=0.1 * tolerance, max_points=0).polygons[0]
                xmax = max(xmax, poly[:, 0].max())
                polys.append(poly)
        return polys


def create_phc_label(pattern_number, size, position, font_prop=None, tolerance=0.1):
    return PolygonOperations.create_polygons_from_str(f"+{pattern_number}", size, position, font_prop, tolerance)


def linker_polygon(pltdata, layer=3):
    return gdspy.Polygon(
        [
            (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0),
            (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos + pltdata.dy / 2.0),
            (pltdata.xpos + pltdata.dx / 2.0, pltdata.ypos + pltdata.dy / 2.0),
            (pltdata.xpos + pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0),
            (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0)
        ], 
        layer=layer
    )


def pinch_pt_polygon(beam_dy, pinch_pt_size, pinch_pt_xloc, pinch_pt_yloc, layer=5):
    _upper_triangle = gdspy.Polygon(
        [
            (pinch_pt_xloc, pinch_pt_yloc + pinch_pt_size / 2.0),
            (pinch_pt_xloc - (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0, pinch_pt_yloc + beam_dy / 2.0),
            (pinch_pt_xloc + (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0, pinch_pt_yloc + beam_dy / 2.0),
            (pinch_pt_xloc, pinch_pt_yloc + pinch_pt_size / 2.0)
        ],
        layer=layer
    )
    _lower_triangle = gdspy.Polygon(
        [
            (pinch_pt_xloc, pinch_pt_yloc - pinch_pt_size / 2.0),
            (pinch_pt_xloc - (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0, pinch_pt_yloc - beam_dy / 2.0),
            (pinch_pt_xloc + (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0, pinch_pt_yloc - beam_dy / 2.0),
            (pinch_pt_xloc, pinch_pt_yloc - pinch_pt_size / 2.0)
        ],
        layer=layer
    )
    return _upper_triangle, _lower_triangle


def pinch_pt_polygon_vertical(beam_dy, pinch_pt_size, pinch_pt_xloc, pinch_pt_yloc, layer=3):
    _left_triangle = gdspy.Polygon(
        [
            (pinch_pt_xloc - pinch_pt_size / 2.0, pinch_pt_yloc),
            (pinch_pt_xloc - beam_dy / 2.0, pinch_pt_yloc - (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0),
            (pinch_pt_xloc - beam_dy / 2.0, pinch_pt_yloc + (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0),
            (pinch_pt_xloc - pinch_pt_size / 2.0, pinch_pt_yloc)
        ], 
        layer=layer
    )
    _right_triangle = gdspy.Polygon(
        [
            (pinch_pt_xloc + pinch_pt_size / 2.0, pinch_pt_yloc),
            (pinch_pt_xloc + beam_dy / 2.0, pinch_pt_yloc - (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0),
            (pinch_pt_xloc + beam_dy / 2.0, pinch_pt_yloc + (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0),
            (pinch_pt_xloc + pinch_pt_size / 2.0, pinch_pt_yloc)
        ], 
        layer=layer
    )
    return _left_triangle, _right_triangle
