from abc import abstractproperty
from dataclasses import dataclass
import enum
from functools import cached_property
import gdspy
import numpy
import math
import copy


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

    def __init__(self, side_length, layer=0, center_to_origin=True):
        super().__init__(side_length, side_length, layer, center_to_origin)
