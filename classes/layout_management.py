import copy
from dataclasses import dataclass
from functools import cached_property

import gdspy
import numpy

from classes.generic import (Component, Rectangle, Square, Unit,
                             create_phc_label, linker_polygon,
                             pinch_pt_polygon, pinch_pt_polygon_vertical)
from classes.random_constants import (BEAM_SPACING, CORNER_BEND_PTS,
                                      CORNER_BEND_RAD, DEVICE_COL_COUNT,
                                      DEVICE_ROW_COUNT, EDGE_OFFSET,
                                      LINKER_EDGE_OFFSET, LINKER_WIDTH,
                                      PAIRS, PHC_GROUP_COUNT, SPACING, SPACING_Y,
                                      WRITE_FIELD_X_SIZE, WRITE_FIELD_Y_SIZE)

# Parameter for the support
SUPPORT_CONNECTOR_WIDTH = 1e3 #0.000001e3
SUPPORT_NOTCH_SIZE = 400
SUPPORT_BEAM_SIZE = 1e3
NUM_X_SUPPORT = 2
NUM_Y_SUPPORT = 10


class PhcGroup:

    """ 
    A group of Phcs.

    Phcs are written in groups to accommodate for any manufacturing errors. Not all Phcs within a group are neccesarily
    the same as test Phcs with no beam holes are created to test light transmittance.

    """


class SupportStructure:

    """
    Tabs used to cut the PhcGroup from the device layer.
    
    """

    OVERLAP_WIDTH = 0
    
    @dataclass
    class PLTdata:
        xpos: float
        ypos: float
        dx: float
        dy: float

    def __init__(self) -> None:
        pass

    @classmethod
    def write_support_region(cls, innerframe, layer=3):

        _support_box_outer = copy.copy(innerframe)
        _ymin = _support_box_outer.ypos - _support_box_outer.dy / 2.0
        _ymax = _support_box_outer.ypos + _support_box_outer.dy / 2.0
        _xmin = _support_box_outer.xpos - _support_box_outer.dx / 2.0
        _xmax = _support_box_outer.xpos + _support_box_outer.dx / 2.0
        _tmp_support_outer = linker_polygon(_support_box_outer, layer)

        _support_box_inner = copy.copy(_support_box_outer)
        _support_box_inner.dx = _support_box_inner.dx - 2 * cls.OVERLAP_WIDTH
        _support_box_inner.dy = _support_box_inner.dy - 2 * cls.OVERLAP_WIDTH
        _tmp_support_inner = linker_polygon(_support_box_inner, layer)

        _write_field = copy.copy(_support_box_outer)
        _write_field.dx = WRITE_FIELD_X_SIZE
        _write_field.dy = WRITE_FIELD_Y_SIZE
        _tmp_write_field = linker_polygon(_write_field, layer)

        _support_outline = copy.copy(_support_box_outer)
        _support_outline.dx = _support_outline.dx + 2 * SUPPORT_CONNECTOR_WIDTH
        _support_outline.dy = _support_outline.dy + 2 * SUPPORT_CONNECTOR_WIDTH
        _tmp_support_outline = linker_polygon(_support_outline, layer)

        _box_write_area = gdspy.fast_boolean(_tmp_support_outer, _tmp_support_inner, 'not', max_points=0, layer=layer)

        _tmp_notch_yarray = (_support_box_outer.dy - CORNER_BEND_RAD * 2 - SUPPORT_BEAM_SIZE) * numpy.linspace(0, 1.0, num=NUM_Y_SUPPORT)
        for i in range(NUM_Y_SUPPORT):
            _tmp_notch_ypos = _ymin * 1.15 + CORNER_BEND_RAD + SUPPORT_BEAM_SIZE / 2.0 + _tmp_notch_yarray[i]
            _tmp_beam_polygon_L = gdspy.Polygon(
                [
                    (_xmin - SUPPORT_CONNECTOR_WIDTH, _tmp_notch_ypos - SUPPORT_BEAM_SIZE / 2.0),
                    (_xmin - SUPPORT_CONNECTOR_WIDTH, _tmp_notch_ypos + SUPPORT_BEAM_SIZE / 2.0),
                    (_xmin, _tmp_notch_ypos + SUPPORT_BEAM_SIZE / 2.0),
                    (_xmin, _tmp_notch_ypos - SUPPORT_BEAM_SIZE / 2.0),
                    (_xmin - SUPPORT_CONNECTOR_WIDTH, _tmp_notch_ypos - SUPPORT_BEAM_SIZE / 2.0)
                ], 
                layer=layer
            )
            _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(
                SUPPORT_BEAM_SIZE, SUPPORT_NOTCH_SIZE,
                _xmin - SUPPORT_CONNECTOR_WIDTH / 2.0,
                _tmp_notch_ypos, layer
            )
            _tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L, _tmp_upper_triangle, 'not', max_points=0, layer=layer)
            _tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L, _tmp_lower_triangle, 'not', max_points=0, layer=layer)

            _tmp_beam_polygon_R = gdspy.Polygon(
                [
                    (_xmax, _tmp_notch_ypos - SUPPORT_BEAM_SIZE / 2.0),
                    (_xmax, _tmp_notch_ypos + SUPPORT_BEAM_SIZE / 2.0),
                    (_xmax + SUPPORT_CONNECTOR_WIDTH, _tmp_notch_ypos + SUPPORT_BEAM_SIZE / 2.0),
                    (_xmax + SUPPORT_CONNECTOR_WIDTH, _tmp_notch_ypos - SUPPORT_BEAM_SIZE / 2.0),
                    (_xmax, _tmp_notch_ypos - SUPPORT_BEAM_SIZE / 2.0)
                ], 
                layer=layer
            )
            _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(
                SUPPORT_BEAM_SIZE, SUPPORT_NOTCH_SIZE,
                _xmax + SUPPORT_CONNECTOR_WIDTH / 2.0,
                _tmp_notch_ypos, layer
            )
            _tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R, _tmp_upper_triangle, 'not', max_points=0, layer=layer)
            _tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R, _tmp_lower_triangle, 'not', max_points=0, layer=layer)

            _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_L, 'or', max_points=0, layer=layer)
            _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_R, 'or', max_points=0, layer=layer)

        _tmp_notch_xarray = (_support_box_outer.dx - CORNER_BEND_RAD * 2 - SUPPORT_BEAM_SIZE) * numpy.linspace(0, 1.0, num=NUM_X_SUPPORT)
        for i in range(NUM_X_SUPPORT):
            _tmp_notch_xpos = _xmin + CORNER_BEND_RAD + SUPPORT_BEAM_SIZE / 2.0 + _tmp_notch_xarray[i]
            _tmp_beam_polygon_T = gdspy.Polygon(
                [
                    (_tmp_notch_xpos - SUPPORT_BEAM_SIZE / 2.0, _ymax),
                    (_tmp_notch_xpos - SUPPORT_BEAM_SIZE / 2.0, _ymax + SUPPORT_CONNECTOR_WIDTH),
                    (_tmp_notch_xpos + SUPPORT_BEAM_SIZE / 2.0, _ymax + SUPPORT_CONNECTOR_WIDTH),
                    (_tmp_notch_xpos + SUPPORT_BEAM_SIZE / 2.0, _ymax),
                    (_tmp_notch_xpos - SUPPORT_BEAM_SIZE / 2.0, _ymax)
                ],
                layer=layer
            )
            _tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(
                SUPPORT_BEAM_SIZE, SUPPORT_NOTCH_SIZE,
                _tmp_notch_xpos,
                _ymax + SUPPORT_CONNECTOR_WIDTH / 2.0,
                layer
            )
            _tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T, _tmp_left_triangle, 'not', max_points=0, layer=layer)
            _tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T, _tmp_right_triangle, 'not', max_points=0, layer=layer)

            _tmp_beam_polygon_B = gdspy.Polygon(
                [
                    (_tmp_notch_xpos - SUPPORT_BEAM_SIZE / 2.0, _ymin - SUPPORT_CONNECTOR_WIDTH),
                    (_tmp_notch_xpos - SUPPORT_BEAM_SIZE / 2.0, _ymin),
                    (_tmp_notch_xpos + SUPPORT_BEAM_SIZE / 2.0, _ymin),
                    (_tmp_notch_xpos + SUPPORT_BEAM_SIZE / 2.0, _ymin - SUPPORT_CONNECTOR_WIDTH),
                    (_tmp_notch_xpos - SUPPORT_BEAM_SIZE / 2.0, _ymin - SUPPORT_CONNECTOR_WIDTH)
                ],
                layer=layer
            )
            _tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(
                SUPPORT_BEAM_SIZE, SUPPORT_NOTCH_SIZE,
                _tmp_notch_xpos,
                _ymin - SUPPORT_CONNECTOR_WIDTH / 2.0,
                layer
            )
            _tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B, _tmp_left_triangle, 'not', max_points=0, layer=layer)
            _tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B, _tmp_right_triangle, 'not', max_points=0, layer=layer)

            _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_T, 'or', max_points=0, layer=layer)
            _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_B, 'or', max_points=0, layer=layer)

        _support_combined = gdspy.fast_boolean(_tmp_write_field, _tmp_support_outline, 'not', max_points=0, layer=layer)
        _support_combined = gdspy.fast_boolean(_support_combined, _box_write_area, 'or', max_points=0, layer=layer)

        return _support_combined

    TEXT_HEIGHT = 8e3

    @classmethod
    def create_outer_frame_with_supports(cls, cell, beamdata, beam_dy_list, outer_box_bounding_box, pattern_number=None, reverse_tone=False, layer=3):

        _ymin = beamdata.ypos - beam_dy_list[0] / 2.0 - EDGE_OFFSET
        _ymax = _ymin + (PHC_GROUP_COUNT - 1) * BEAM_SPACING + numpy.sum(beam_dy_list) + EDGE_OFFSET * 2
        _tmp_beamdata = copy.copy(beamdata)
        _support_inner = cls.PLTdata(
            xpos=_tmp_beamdata.xpos, 
            ypos=(_ymin + _ymax) / 2.0, 
            dx=beamdata.dx + 2 * LINKER_WIDTH + 2 * LINKER_EDGE_OFFSET, 
            dy=_ymax - _ymin
        )

        _frame_write_area = cls.write_support_region(_support_inner, layer)

        if pattern_number is not None:
            _pattern_number_to_write = create_phc_label(pattern_number, cls.TEXT_HEIGHT, (_support_inner.xpos-8e3, _support_inner.ypos + _support_inner.dy / 2.0 + 2.5e3))
            _frame_write_area = gdspy.fast_boolean(_frame_write_area, _pattern_number_to_write, 'not', max_points=0, layer=layer)

        if not reverse_tone:
            bounds = _frame_write_area.get_bounding_box()

            shift = 5000

            box = gdspy.Rectangle(numpy.array(bounds[0]) - shift, numpy.array(bounds[1]) + shift, layer=layer)
            bounds_inner = outer_box_bounding_box
            box_inner = gdspy.Rectangle(bounds_inner[0], bounds_inner[1], layer=layer)
            box_inner = box_inner.fillet(CORNER_BEND_RAD, points_per_2pi=CORNER_BEND_PTS)

            _box_write_area_reverse = gdspy.fast_boolean(box, _frame_write_area, 'not', max_points=0, layer=layer)
            _box_write_area_reverse = gdspy.fast_boolean(_box_write_area_reverse, box_inner, 'not', max_points=0, layer=layer)
            cell.add(_box_write_area_reverse)
        else:
            cell.add(_frame_write_area)

class WorkArea(Component):

    """
    This defines the maximum area for the chip.

    """

    class PositionId:
        
        CELL = gdspy.Cell("position_ids_cell")
        LAYER = 0
        LENGTH = 20 * Unit.um.value
        DIMENSION_SCALAR = 1.2  # EPBG Alignment Marks (previously ebpg_mark_factor) | is this half of the width/height times the scalar? TODO

        def __init__(self, name: str, coords: set, add_orienter: bool) -> None:
            self.name = name
            self.add_orienter = add_orienter
            self.coords = coords  # no need to normalize to the origin because the square in the POSITION_ID_CELL is already centered

        def get_referenced_cells(self):
            referenced_cells = [gdspy.CellReference(self.CELL, self.coords)]
            if self.add_orienter:
                referenced_cells += self.get_orienters()
            return referenced_cells
        
        def get_orienters(self):
            x_sign = 1 if self.coords[0] > 0 else -1
            y_sign = 1 if self.coords[1] > 0 else -1
            return [
                gdspy.CellReference(self.CELL, (self.coords[0], self.coords[1] + 600 * Unit.um.value * y_sign)),
                gdspy.CellReference(self.CELL, (self.coords[0] + 600 * Unit.um.value * x_sign, self.coords[1])),
            ]
        
        @classmethod
        def create_position_id_rect(cls):
            return cls.CELL.add(Square(side_length=cls.LENGTH, layer=cls.LAYER))

        @classmethod
        def create_position_ids(cls, manufacturing_boundary, device_writing_frame):
            """ These are the small squares outside the corners of the device_writing_frame. """
            cls.create_position_id_rect()
            position_ids = [
                cls(name="top_right",    coords=(device_writing_frame.max_x * cls.DIMENSION_SCALAR, device_writing_frame.max_y * cls.DIMENSION_SCALAR), add_orienter=False),
                cls(name="bottom_right", coords=(device_writing_frame.max_x * cls.DIMENSION_SCALAR, device_writing_frame.min_y * cls.DIMENSION_SCALAR), add_orienter=False),
                cls(name="top_left",     coords=(device_writing_frame.min_x * cls.DIMENSION_SCALAR, device_writing_frame.max_y * cls.DIMENSION_SCALAR), add_orienter=False),
                cls(name="bottom_left",  coords=(device_writing_frame.min_x * cls.DIMENSION_SCALAR, device_writing_frame.min_y * cls.DIMENSION_SCALAR), add_orienter=True),
            ]
            for position_id in position_ids:
                for referenced_cell in position_id.get_referenced_cells():
                    manufacturing_boundary = gdspy.boolean(manufacturing_boundary, referenced_cell, "not", layer=cls.LAYER)
            return manufacturing_boundary
        
    class LithoMark:
        
        CELL = gdspy.Cell("litho_marks_cell")
        LAYER = 1
        HEIGHT = 20 * Unit.um.value  # previously mark_len
        WIDTH = 4 * Unit.um.value    # previously mark_cross
        DIMENSION_SCALAR = 1.4       # EPBG Alignment Marks (previously ebpg_mark_factor) | is this half of the width/height times the scalar? TODO
        SIDE_COUNT = 32              # The amount of litho marks there should be for each side

        def __init__(self) -> None:
            pass

        @classmethod
        def create_litho_mark_rects(cls):
            """ Just like in create_position_id_rect, we reference the same geometry in some reference cell. """
            cls.CELL.add(
                [
                    Rectangle(width=cls.WIDTH, height=cls.HEIGHT, layer=cls.LAYER),  # h_rect
                    Rectangle(width=cls.HEIGHT, height=cls.WIDTH, layer=cls.LAYER),  # v_rect
                ]
            )

        @classmethod
        def create_litho_marks(cls, manufacturing_boundary, device_writing_frame):
            """ Technically we don't need manufacturing_boundary... TODO (previously litho_marks). """
            cls.create_litho_mark_rects()
            x_coords = numpy.linspace(device_writing_frame.min_x * cls.DIMENSION_SCALAR, device_writing_frame.max_x * cls.DIMENSION_SCALAR, num=cls.SIDE_COUNT)
            y_coords = numpy.linspace(device_writing_frame.min_y * cls.DIMENSION_SCALAR, device_writing_frame.max_y * cls.DIMENSION_SCALAR, num=cls.SIDE_COUNT)
            referenced_cells = []
            for y_i, y_coord in enumerate(y_coords):
                for x_i, x_coord in enumerate(x_coords):
                    if 0 < y_i < len(y_coords) -1 and 0 < x_i < len(x_coords) -1:
                        continue
                    referenced_cells.append(gdspy.CellReference(cls.CELL, (x_coord, y_coord)))
                    WorkArea.god_cell.add(referenced_cells[-1])
            return manufacturing_boundary

    def __init__(self, x_length, y_length, layer, add_position_ids=True, add_litho_marks=True) -> None:
        super().__init__()
        self.x_length = x_length
        self.y_length = y_length
        self.layer = layer
        self.add_position_ids = add_position_ids
        self.add_litho_marks = add_litho_marks

    @cached_property
    def manufacturing_boundary(self):
        """ The maximum chip area (previously: outer_frame). """
        return Rectangle(self.x_length, self.y_length, layer=4)
    
    @cached_property
    def device_writing_frame(self):
        """ The main area the devices are written (previously: inner_frame). """
        return Rectangle(SPACING * Unit.um.value * DEVICE_ROW_COUNT, SPACING_Y * Unit.um.value * DEVICE_COL_COUNT * PAIRS, layer=self.layer)

    @cached_property
    def polygon(self):
        """ # The main area the devices are written (previously: frame). """
        polygon = gdspy.boolean(self.manufacturing_boundary, self.device_writing_frame, "not", layer=self.layer)
        if self.add_position_ids:
            polygon = self.PositionId.create_position_ids(polygon, self.device_writing_frame)
        if self.add_litho_marks:
            polygon = self.LithoMark.create_litho_marks(polygon, self.device_writing_frame)
        return polygon
