import copy
import math
from dataclasses import dataclass
from functools import cached_property

import gdspy
import numpy

from classes.generic import (Component, Polygon, Rectangle, Unit, better_dataclass,
                             create_phc_label, linker_polygon,
                             pinch_pt_polygon, pinch_pt_polygon_vertical)
from classes.layout_management import SupportStructure
from classes.random_constants import (BEAM_SPACING, CORNER_BEND_PTS,
                                      CORNER_BEND_RAD, DEVICE_COL_COUNT,
                                      DEVICE_ROW_COUNT, EDGE_OFFSET, GAP,
                                      HEIGHT, LENGTH, LINKER_BEAM_SIZE,
                                      LINKER_EDGE_OFFSET, LINKER_NOTCH_SIZE,
                                      LINKER_WIDTH, LINKER_X_NOTCHES, PAIRS,
                                      PHC_GROUP_COUNT, SPACING,
                                      SPACING_Y, SUPPORT_PINCH_LENGTH,
                                      SUPPORT_PINCH_WIDTH, SUPPORT_PITCH,
                                      SUPPORT_WIDTH, WRITE_FIELD_X_SIZE,
                                      WRITE_FIELD_Y_SIZE)

CELL_MAIN = gdspy.Cell("main_cell")
GRATING_CELL = gdspy.Cell("grating")
REFLECTOR_CELL = gdspy.Cell("reflector")
GRATING_TAPER_CELL = gdspy.Cell("grating_taper")
REFLECTOR_TAPER_CELL = gdspy.Cell("reflector_taper")
DEVS = 6
SHOT_PITCH = 0.1
NUM_CIRC_POINTS = 181

# Parameters for nanobeam
HOLE_POS_OFFSET = 0

# For grating spacer
SPACER = 0  # Spacer size

UNDER_EDGE_SIZE = 12

# Parameters for the linker
LINKER_CONNECTOR_WIDTH = 1e3
LINKER_NOTCH_SIZE = 500
LINKER_BEAM_SIZE = 500
LINKER_X_NOTCHES = 0
LINKER_Y_NOTCHES = 0

# Parameters for circular grating structure
PERIOD = 777e-9
DUTY_CYCLE = 0.376
GRATING_LINE_WIDTH = numpy.round(PERIOD * (1 - DUTY_CYCLE) / 1e-9, 1)
GRATING_SPACING = numpy.round(PERIOD * DUTY_CYCLE / 1e-9, 1)
NUM_CIRC_GRATING_POINTS = 181
NUM_GRATING = 2
CIRC_GRATING_BASE = 1500
GRATING_ANGLE = 10 * numpy.pi / 180  # 15 degree in radian
ODD_SUPPORT_ANGLE = 75 * numpy.pi / 180  # 55 degree in radian
EVEN_SUPPORT_ANGLE = numpy.pi / 2.0  # 90 degree in radian

SUPPORT_ANGLE_WIDTH = 2.3 * numpy.pi / 180  # 10 degree in radian


class Hole:
    
    @better_dataclass
    class Config:
        hole_width_variations             : list  # previously hx_list
        hole_height_variations            : list  # previously hy_list
        hole_spacing_variations           : list  # previously acav
        hole_group_padding_variations     : list  # previously amir
        hole_modification_func_variations : any   # previously 
        hole_count_variations             : list  # previously mirror_list
        hole_count_with_acav_spacing      : list  # number of holes with acav spacing (previously cavity_list)
        taper_hole_count                  : list 
        taper_hole_modification_func      : any 

    """
    This is the cavity that the photons bounce around in.

    """


class NanoBeam(Polygon):

    """
    The area that houses the holes.

    Note:    
        This is a 1 dimmensional beam. There can be more in the future. TODO
    """

    @better_dataclass
    class Config:
        width_variations  : float  # previously wy_list
        height_variations : float  # previously 


    def __init__(self, points, layer=0, datatype=0):
        super().__init__(points, layer, datatype)

    @staticmethod
    def write_beam_single(cell, pltdata, layer=1):
        cell.add(gdspy.Polygon([
            (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0),
            (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos + pltdata.dy / 2.0),
            (pltdata.xpos + pltdata.dx / 2.0, pltdata.ypos + pltdata.dy / 2.0),
            (pltdata.xpos + pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0),
            (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0)
        ], layer=layer))

    @classmethod
    def exponential_cavity_period_list(cls, amir, acav, num_cav, exponent):
        N = int(num_cav / 2)
        cavity_period_list = numpy.zeros(N)
        a = (amir - acav) / ((N - 1) ** exponent)
        b = acav

        for i in range(N):
            cavity_period_list[i] = a * (i ** exponent) + b

        cavity_period_list = numpy.append(numpy.flipud(cavity_period_list), cavity_period_list)
        return cavity_period_list

    @staticmethod
    def get_tapers(num_end_taper, end_period, amir):
        taper_period_list = numpy.linspace(amir, end_period, num_end_taper + 1)
        return taper_period_list[1:]

    @staticmethod
    def write_hole_single(pltdata, layer=2):
        _philist = numpy.linspace(0, 2 * numpy.pi, NUM_CIRC_POINTS)
        _circ_pts = [
            ((round(pltdata.xpos / SHOT_PITCH) + round(pltdata.dx / 2 * numpy.cos(phi) / SHOT_PITCH)) * SHOT_PITCH,
            (round(pltdata.ypos / SHOT_PITCH) + round(pltdata.dy / 2 * numpy.sin(phi) / SHOT_PITCH)) * SHOT_PITCH)
            for phi in _philist
        ]
        return gdspy.Polygon(_circ_pts, layer=layer)

    DO_CUSTOM_CAVITY = True
    SEGMENTS = numpy.array([1.57385141e-07, 1.58156287e-07, 1.59014195e-07, 1.59789244e-07, 1.61512865e-07, 1.63789232e-07, 1.65369591e-07, 1.67000125e-07]) * 1e9
    CUSTOM_CAVITY_LIST = numpy.append(numpy.flip(SEGMENTS), SEGMENTS)

    @classmethod
    def write_hole_1D(cls, cell, beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, amir, end_period,
                    end_taper_L=False, end_taper_R=False, num_end_taper=0, reverse_tone=False, layer=2):
        n = 1.5  # Exponent for polynomial

        _one_side_holes = (num_cav_hole - 1) / 2.0
        _idx_list = numpy.linspace(-_one_side_holes, _one_side_holes, int(num_cav_hole))
        X1 = _idx_list[int(num_cav_hole / 2)]  # First Hole
        XN = _idx_list[-1]  # Last Hole

        _acav_list = cls.exponential_cavity_period_list(amir, acav, num_cav_hole, exponent=2)

        if cls.DO_CUSTOM_CAVITY:
            _acav_list = cls.CUSTOM_CAVITY_LIST

        if end_taper_L == True:
            taper_periods = cls.get_tapers(num_end_taper, end_period, amir=amir)

        _amir_list_L = numpy.append(numpy.flipud(taper_periods), [amir] * int(num_mir_hole_L))

        _amir_list_R = [amir] * int(num_mir_hole_R)
        _amir_list_R.extend(taper_periods)
        _aper_list = list(copy.copy(_amir_list_L))

        _aper_list.extend(_acav_list)
        _aper_list.extend(_amir_list_R)

        _hole_write = copy.copy(holedata)
        _hole_write.xpos = _hole_write.xpos - numpy.sum(numpy.array(_acav_list)) / 2.0 - numpy.sum(
            numpy.array(_amir_list_L))
        _taper_scale_list = []
        if num_end_taper > 0:
            _taper_scale_list = numpy.linspace(0, 1.0, num_end_taper + 2)
            _taper_scale_list_L = _taper_scale_list[1:-1]
            _taper_scale_list_R = numpy.flipud(_taper_scale_list[1:-1])

        for i in range(len(_aper_list)):
            _hole_write.xpos = _hole_write.xpos + _aper_list[i] / 2.0
            if i < num_end_taper * end_taper_L:
                _hole_write.dx = holedata.dx * 1
                _hole_write.dy = holedata.dy * 1
            else:
                _hole_write.dx = holedata.dx
                _hole_write.dy = holedata.dy
            if i >= (len(_aper_list) - num_end_taper * end_taper_R):
                _hole_write.dx = holedata.dx * 1
                _hole_write.dy = holedata.dy * 1
            _single_hole_polygon = cls.write_hole_single(_hole_write)

            if i == 0:
                _hole_polygon_combined = _single_hole_polygon
            else:
                _hole_polygon_combined = gdspy.fast_boolean(_hole_polygon_combined, _single_hole_polygon, 'or', max_points=0, layer=layer)
            _hole_write.xpos = _hole_write.xpos + _aper_list[i] / 2.0

        if reverse_tone:
            _tmp_beam_polygon = gdspy.Polygon(
                [
                    (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
                    (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
                    (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
                    (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
                    (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
                ],
                layer=layer
            )
            _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _hole_polygon_combined, 'not', max_points=0, layer=layer)
            send = _tmp_beam_polygon
        else:
            send = _hole_polygon_combined
        return send

    UNDER_RING_SIZE = 12
    OVER_RING_SIZE = 12
    OVER_EDGE_SIZE = 12

    @classmethod
    def create_nano_beam(
            cls, i, _initial_ypos, _tmp_beamdata, cell, holedata, beam_dy_list, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, amir,
            end_period, guides, end_taper_L=False, end_taper_R=False, num_end_taper=0, reverse_tone=False, edge_overdose=False,
            ring_overdose=False, ring_underdose=True, edge_underdose=True
        ):
        """ (previously write_hole_2D). """
        holedata.ypos = _initial_ypos + i * BEAM_SPACING + numpy.sum(beam_dy_list[:i]) + beam_dy_list[i] / 2.0
        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0
        _tmp_beamdata.dy = beam_dy_list[i]

        if guides == 1 and i == 0:
            cls.write_beam_single(cell, _tmp_beamdata, layer=2)
        elif guides == 2 and (i == 0 or i == (PHC_GROUP_COUNT-1)):
            cls.write_beam_single(cell, _tmp_beamdata, layer=2)
        else:
            ord_holes = cls.write_hole_1D(
                cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, 
                amir, end_period, end_taper_L, end_taper_R, num_end_taper, reverse_tone=reverse_tone
            )
            if ring_overdose or edge_overdose is True:
                if reverse_tone == True:
                    holedata.dx = holedata.dx + cls.OVER_RING_SIZE * ring_overdose * 2
                    holedata.dy = holedata.dy + cls.OVER_RING_SIZE * ring_overdose * 2
                    _tmp_beamdata.dy = _tmp_beamdata.dy - cls.OVER_EDGE_SIZE * 2 * edge_overdose
                    overdose_holes = cls.write_hole_1D(
                        cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, 
                        amir, end_period, end_taper_L, end_taper_R, num_end_taper, reverse_tone=reverse_tone
                    )
                    region = gdspy.fast_boolean(ord_holes, overdose_holes, 'not', max_points=0, layer=11)
                    cell.add(region)
                    cell.add(overdose_holes)
                    holedata.dx = holedata.dx - cls.OVER_RING_SIZE * ring_overdose * 2
                    holedata.dy = holedata.dy - cls.OVER_RING_SIZE * ring_overdose * 2
            elif ring_underdose or edge_underdose is True:
                if reverse_tone is False:
                    ref_ord_edge = cls.write_hole_1D(
                        cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, 
                        amir, end_period, end_taper_L, end_taper_R, num_end_taper, reverse_tone=True
                    )
                    holedata.dx = holedata.dx - cls.UNDER_RING_SIZE * ring_underdose
                    holedata.dy = holedata.dy - cls.UNDER_RING_SIZE * ring_underdose
                    _tmp_beamdata.dy = _tmp_beamdata.dy + UNDER_EDGE_SIZE * 2 * edge_underdose
                    underdose_holes = cls.write_hole_1D(
                        cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, 
                        amir, end_period, end_taper_L, end_taper_R, num_end_taper, reverse_tone=reverse_tone
                    )
                    _tmp_beamdata.dx = _tmp_beamdata.dx - SPACER * 2
                    underdose_edge = cls.write_hole_1D(
                        cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav,
                        amir, end_period, end_taper_L, end_taper_R, num_end_taper, reverse_tone=True
                    )
                    region = gdspy.fast_boolean(ord_holes, underdose_holes, 'not', max_points=0, layer=12)
                    region2 = gdspy.fast_boolean(underdose_edge, ref_ord_edge, 'not', max_points=0, layer=12)
                    # underdose_holes = gdspy.fast_boolean(underdose_edge, region2, 'not', max_points=0, layer=12)
                    cell.add(underdose_holes)
                    cell.add(region)
                    cell.add(region2)
            else:
                cell.add(ord_holes)
        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0 + BEAM_SPACING


class GratingCoupler:

    """
    Where the photons are inputted.
    
    """

    def __init__(self, beam_position, layer) -> None:

        self.beam_position = beam_position
        self.layer = layer

        self.left_shockwave = self.create_left_circ_grating(self.beam_position, layer=layer)
        self.right_shockwave = self.create_right_circ_grating(self.beam_position, layer=layer)
    
    @classmethod
    def create_shockwaves(cls, beamdata, is_left, layer):
        """ idk what to call this, but this is a better name (was write_left_circ_grating/write_right_circ_grating). """

        direction_mult = -1 if is_left else 1
        some_angle = (numpy.pi * 3.0) if is_left else numpy.pi

        phi_list = numpy.append(
            numpy.linspace(
                -1 * direction_mult * numpy.pi / 2.0, 
                some_angle / 2.0, 
                NUM_CIRC_GRATING_POINTS - 1
            ), 
            -1 * direction_mult * numpy.pi / 2.0
        )
        _ini_pt = [(beamdata.xpos + (direction_mult * beamdata.dx / 2), beamdata.ypos)]
        _tmp_beamdata = copy.copy(beamdata)
        _radius_inner = CIRC_GRATING_BASE

        grating = None
        for i in range(NUM_GRATING + 1):
            grating_inner = gdspy.Polygon(
                [
                    (
                        (round((_tmp_beamdata.xpos + (direction_mult * _tmp_beamdata.dx / 2)) / SHOT_PITCH) + round(_radius_inner * numpy.cos(phi) / SHOT_PITCH)) * SHOT_PITCH,
                        (round(_tmp_beamdata.ypos / SHOT_PITCH) + round(_radius_inner * numpy.sin(phi) / SHOT_PITCH)) * SHOT_PITCH
                    )
                    for phi in phi_list
                ], 
                layer=layer
            )

            _radius_outer = _radius_inner + GRATING_SPACING
            grating_outer = gdspy.Polygon(
                [
                    (
                        (round((_tmp_beamdata.xpos + (direction_mult * _tmp_beamdata.dx / 2)) / SHOT_PITCH) + round(_radius_outer * numpy.cos(phi) / SHOT_PITCH)) * SHOT_PITCH,
                        (round(_tmp_beamdata.ypos / SHOT_PITCH) + round(_radius_outer * numpy.sin(phi) / SHOT_PITCH)) * SHOT_PITCH
                    )
                    for phi in phi_list
                ],
                layer=layer
            )
            grating_tmp = gdspy.fast_boolean(grating_outer, grating_inner, 'not', max_points=0, layer=layer)
            _radius_outer = _radius_outer + 10
            support_angle = ODD_SUPPORT_ANGLE if (i % 2 == 0) else EVEN_SUPPORT_ANGLE  # this is stupid
            _philist_support = numpy.linspace(
                -1 * direction_mult * numpy.pi / 2.0 + support_angle - SUPPORT_ANGLE_WIDTH / 2.0,
                -1 * direction_mult * numpy.pi / 2.0 + support_angle + SUPPORT_ANGLE_WIDTH / 2.0,
                NUM_CIRC_GRATING_POINTS - 1
            )
            _support_pt = [
                (
                    (round((_tmp_beamdata.xpos + (direction_mult *_tmp_beamdata.dx / 2)) / SHOT_PITCH) + round(_radius_outer * numpy.cos(phi) / SHOT_PITCH)) * SHOT_PITCH,
                    (round(_tmp_beamdata.ypos / SHOT_PITCH) + round(_radius_outer * numpy.sin(phi) / SHOT_PITCH)) * SHOT_PITCH
                )
                for phi in _philist_support
            ]
            _support_pt_combined = copy.copy(_ini_pt)
            _support_pt_combined.extend(_support_pt)
            _support_pt_combined.extend(_ini_pt)
            _support_frame = gdspy.Polygon(_support_pt_combined, layer=layer)
            grating_tmp = gdspy.fast_boolean(grating_tmp, _support_frame, 'not', max_points=0, layer=layer)
            if (i % 2 == 0):  # even
                _philist_support = numpy.linspace(
                    some_angle / 2.0 - ODD_SUPPORT_ANGLE - SUPPORT_ANGLE_WIDTH / 2.0,
                    some_angle / 2.0 - ODD_SUPPORT_ANGLE + SUPPORT_ANGLE_WIDTH / 2.0,
                    NUM_CIRC_GRATING_POINTS - 1
                )
                _support_pt = [
                    (
                        (round((_tmp_beamdata.xpos + (direction_mult * _tmp_beamdata.dx / 2)) / SHOT_PITCH) + round(_radius_outer * numpy.cos(phi) / SHOT_PITCH)) * SHOT_PITCH,
                        (round(_tmp_beamdata.ypos / SHOT_PITCH) + round(_radius_outer * numpy.sin(phi) / SHOT_PITCH)) * SHOT_PITCH
                    )
                    for phi in _philist_support
                ]
                _support_pt_combined = copy.copy(_ini_pt)
                _support_pt_combined.extend(_support_pt)
                _support_pt_combined.extend(_ini_pt)
                _support_frame = gdspy.Polygon(_support_pt_combined, layer=layer)
                grating_tmp = gdspy.fast_boolean(grating_tmp, _support_frame, 'not', max_points=0, layer=layer)
            if grating is None:
                grating = grating_tmp
            else:
                grating = gdspy.fast_boolean(grating, grating_tmp, 'or', max_points=0, layer=layer)
            _radius_inner = _radius_outer + GRATING_LINE_WIDTH
            _philist_frame = numpy.linspace(-1 * direction_mult * numpy.pi / 2.0 + GRATING_ANGLE, some_angle / 2.0 - GRATING_ANGLE, NUM_CIRC_GRATING_POINTS - 1)
            frame_pt = [
                (
                    (round((_tmp_beamdata.xpos + (direction_mult * _tmp_beamdata.dx / 2)) / SHOT_PITCH) + round(_radius_outer * numpy.cos(phi) / SHOT_PITCH)) * SHOT_PITCH,
                    (round(_tmp_beamdata.ypos / SHOT_PITCH) + round(_radius_outer * numpy.sin(phi) / SHOT_PITCH)) * SHOT_PITCH
                )
                for phi in _philist_frame
            ]
        frame_pt_combined = copy.copy(_ini_pt)
        frame_pt_combined.extend(frame_pt)
        frame_pt_combined.extend(_ini_pt)
        frame = gdspy.Polygon(frame_pt_combined, layer=layer)
        grating = gdspy.fast_boolean(grating, frame, 'and', max_points=0, layer=layer)
        return grating

    @classmethod
    def create_left_circ_grating(cls, beamdata, layer):
        return cls.create_shockwaves(beamdata, is_left=True, layer=layer)

    @classmethod
    def create_right_circ_grating(cls, beamdata, layer):
        return cls.create_shockwaves(beamdata, is_left=False, layer=layer)

    @classmethod
    def create_circ_grating(cls, _circ_grating_combined, _tmp_beamdata, i, beam_dy_list, layer=5):
        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0
        _tmp_beamdata.dy = beam_dy_list[i]

        _circ_grating_L = cls.create_left_circ_grating(_tmp_beamdata, layer=layer)
        _circ_grating_R = cls.create_right_circ_grating(_tmp_beamdata, layer=layer)
        if i == 0:
            _circ_grating_combined = _circ_grating_L
            _circ_grating_combined = gdspy.fast_boolean(_circ_grating_combined, _circ_grating_R, 'or', max_points=0, layer=layer)
        else:
            _circ_grating_combined = gdspy.fast_boolean(_circ_grating_combined, _circ_grating_L, 'or', max_points=0, layer=layer)
            _circ_grating_combined = gdspy.fast_boolean(_circ_grating_combined, _circ_grating_R, 'or', max_points=0, layer=layer)
        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0 + BEAM_SPACING
        return _circ_grating_combined


class LinearPattern:

    """
    This is for defining how many copies of a Phc you want for compensating for manufacturing errors.

    """

    def __init__(self, x_count=1, y_count=1) -> None:
        self.x_count = x_count
        self.y_count = y_count


class Phc(Component):

    """
    The photonic crystal cavity.

    """

    SUPPORTS_MASK_BOOL = False

    CELL_MEMBRANE_ROW = gdspy.Cell('membrane_width_sweep')  # Single instance of a sweep over the membrane widths
    CELL_MEMBRANE_ROW_NOHOLE = gdspy.Cell('membrane_width_sweep_nohole')  # Single instance of a sweep over the membrane widths
    REMOVE_INNER_MEM = True

    # General Pattern Layout
    WRITE_QFC = False
    WRITE_PHC = True
    PHC_SWEEPER = not WRITE_PHC
    WRITE_MEMS = False

    DEV_LIST = numpy.arange(0, DEVICE_ROW_COUNT * DEVICE_COL_COUNT, 1)

    ACAV = [155.5]
    AMIR = [172.1]

    WY_LIST = [470]
    WY_END_TAPER_LIST = []

    HX_LIST = [63, 66, 69, 73, 76, 79]
    HY_LIST = [180, 184, 188, 192, 196, 200]

    MIRROR_LIST = [10]
    CAVITY_LIST = [16]
    TAPER_HOLES = [0]

    PARAM_SWEEP = list(
        zip(
            DEV_LIST, 
            numpy.array(numpy.meshgrid(ACAV, AMIR, WY_LIST, MIRROR_LIST, CAVITY_LIST, HX_LIST, HY_LIST, TAPER_HOLES)).T.reshape(-1, 8).astype(int)
        )
    )

    @dataclass
    class PLTdata:
        xpos: float
        ypos: float
        dx: float
        dy: float

    phc_instances = {}

    BOUNDING_RECTANGLE_SCALAR = 1.2
    END_PERIOD = numpy.round(955 / (3.1 * 2), 1)
    NUM_GUIDES = 0  # Define the number of blank waveguides for control measurements (0,1,2)

    # For grating spacer
    GRATING_SPACER = False
    # Tone
    REVERSE_TONE = True

    RING_R = 25
    COUPLING_L = 25
    GRATING_OFFSET = 10
    RING_POS_Y_OFFSET = 20
    
    PHC_SWEEP = gdspy.Cell('phc_sweep')
    PHC_SWEEP_FRAME_CELL = gdspy.Cell('phc_sweep_frame')

    # Membrane parameters
    PERFORATIONS = False
    TEXT = True

    # QFC Device Parameters
    WIDTHS = numpy.linspace(0, 1, DEVS)
    GAPS = numpy.linspace(0, 1, DEVS)

    PHC_Y_OFFSET = (35 - 0.275) * Unit.um.value  # An arbitrary parameter for centering the devices at the origin

    def __init__(self, device_id, row_i, col_i, beam, grating_coupler, enable_bounding_rectangle=True) -> None:
        super().__init__()

        self.device_id = device_id
        self.row_i = row_i
        self.col_i = col_i

        self.key = f"phc_{len(self.phc_instances)}"
        self.beam = beam                        # should be referenced as ReferenecCells
        self.grating_coupler = grating_coupler  # should be referenced as ReferenecCells
        self.enable_bounding_rectangle = enable_bounding_rectangle
        self.phc_instances[self.key] = self

        # ---- attributes --------------------
        name = DEVICE_COL_COUNT * row_i + col_i
        gap = self.GAPS[col_i]
        params = self.PARAM_SWEEP[name]
        width = self.WIDTHS[row_i]
        # ------------------------------------

        xpos = (row_i * SPACING - SPACING * (len(self.WIDTHS) - 1) / 2) * Unit.um.value
        ypos = (SPACING_Y / 2 + col_i * SPACING_Y - SPACING_Y * len(self.GAPS) * PAIRS / 2) * Unit.um.value
        # QFC Device
        if self.WRITE_QFC:
            device = self.defineDevice(width, self.RING_R, self.COUPLING_L, gap, self.GRATING_OFFSET, g)  # what is g? TODO
            device_pos = gdspy.CellReference(device, (xpos, ypos + self.RING_POS_Y_OFFSET))
            self.god_cell.add(device_pos)

        if self.PHC_SWEEPER and not self.WRITE_PHC and device_id == 0:
            phc_pos = gdspy.CellReference(self.PHC_SWEEP, (xpos, ypos - self.PHC_Y_OFFSET))
            self.god_cell.add(phc_pos)

        if self.WRITE_PHC:
            phc = self.create_phc_group(params, end_period=self.END_PERIOD, blank_guides=self.NUM_GUIDES, create_label=self.TEXT)
            phc_pos = gdspy.CellReference(phc, (xpos, ypos - self.PHC_Y_OFFSET))
            self.god_cell.add(phc_pos)

        if self.WRITE_MEMS and device_id == 0:
            membrane = self.defineMembrane(f"{col_i}{row_i}", SPACING, SPACING_Y, LENGTH, HEIGHT, layer=5)

            phc_sweep_frame = gdspy.Rectangle((-LENGTH/2, -HEIGHT/2), (LENGTH/2, HEIGHT/2), layer=5)
            self.PHC_SWEEP_FRAME_CELL.add(phc_sweep_frame)
            phc_sweep_frame_cell_pos = gdspy.CellReference(self.PHC_SWEEP_FRAME_CELL, (xpos, ypos), magnification=1000)
            membrane_pos = gdspy.CellReference(membrane, (xpos, ypos), magnification=1000)

            if self.REMOVE_INNER_MEM:
                membrane_pos = gdspy.boolean(membrane_pos, phc_sweep_frame_cell_pos, "not", layer=5)
            self.god_cell.add(membrane_pos)

        # # Supports Mask
        if self.SUPPORTS_MASK_BOOL and row_i == 0 and col_i == 0:
            self.supportsMask(name, SPACING)

    @cached_property
    def bounding_box_cell(self):
        return gdspy.Cell(f"{self.key}_bounding_box_cell")

    @cached_property
    def bounding_rectangle(self):
        rect = Rectangle(
            width=self.grating_coupler.width * self.BOUNDING_RECTANGLE_SCALAR, 
            height=self.grating_coupler.height * self.BOUNDING_RECTANGLE_SCALAR,
        )
        self.bounding_box_cell.add(rect)
        return rect
        
    @classmethod
    def write_linker_region(cls, beamdata, xmin, xmax, ymin, ymax, round_corner=False, layer=6):

        _tmp_beamdata = copy.copy(beamdata)

        _linkerdata_inner = cls.PLTdata(xpos=_tmp_beamdata.xpos, ypos=(ymin + ymax) / 2.0, dx=_tmp_beamdata.dx, dy=ymax - ymin - 2.0 * LINKER_EDGE_OFFSET - 2.0 * LINKER_CONNECTOR_WIDTH)
        _tmp_linker_inner = linker_polygon(_linkerdata_inner, layer)

        _linkerdata_outer = cls.PLTdata(xpos=_tmp_beamdata.xpos, ypos=(ymin + ymax) / 2.0, dx=_tmp_beamdata.dx + 2.0 * LINKER_WIDTH, dy=ymax - ymin - 2.0 * LINKER_EDGE_OFFSET)
        _tmp_linker_outer = linker_polygon(_linkerdata_outer, layer)

        if round_corner is True:
            # _tmp_linker_inner = _tmp_linker_inner.fillet(corner_bend_rad, points_per_2pi=corner_bend_pts)
            _tmp_linker_outer = _tmp_linker_outer.fillet(CORNER_BEND_RAD, points_per_2pi=CORNER_BEND_PTS)

        _tmp_notch_yarray = (_linkerdata_outer.dy - CORNER_BEND_RAD * 2 - LINKER_BEAM_SIZE) * numpy.linspace(0, 1.0,
                                                                                                            num=LINKER_Y_NOTCHES)
        for i in range(LINKER_Y_NOTCHES):
            _tmp_notch_ypos = ymin + LINKER_EDGE_OFFSET + CORNER_BEND_RAD + LINKER_BEAM_SIZE / 2.0 + _tmp_notch_yarray[i]
            _tmp_beam_polygon_L = gdspy.Polygon([
                (xmin, _tmp_notch_ypos - LINKER_BEAM_SIZE / 2.0),
                (xmin, _tmp_notch_ypos + LINKER_BEAM_SIZE / 2.0),
                (xmin + LINKER_EDGE_OFFSET, _tmp_notch_ypos + LINKER_BEAM_SIZE / 2.0),
                (xmin + LINKER_EDGE_OFFSET, _tmp_notch_ypos - LINKER_BEAM_SIZE / 2.0),
                (xmin, _tmp_notch_ypos - LINKER_BEAM_SIZE / 2.0)
            ], layer=layer)
            _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(LINKER_BEAM_SIZE, LINKER_NOTCH_SIZE, xmin + LINKER_EDGE_OFFSET / 2.0, _tmp_notch_ypos, layer)
            _tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L, _tmp_upper_triangle, 'not', max_points=0, layer=layer)
            _tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L, _tmp_lower_triangle, 'not', max_points=0, layer=layer)

            _tmp_beam_polygon_R = gdspy.Polygon(
                [
                    (xmax - LINKER_EDGE_OFFSET, _tmp_notch_ypos - LINKER_BEAM_SIZE / 2.0),
                    (xmax - LINKER_EDGE_OFFSET, _tmp_notch_ypos + LINKER_BEAM_SIZE / 2.0),
                    (xmax, _tmp_notch_ypos + LINKER_BEAM_SIZE / 2.0),
                    (xmax, _tmp_notch_ypos - LINKER_BEAM_SIZE / 2.0),
                    (xmax - LINKER_EDGE_OFFSET, _tmp_notch_ypos - LINKER_BEAM_SIZE / 2.0)
                ],
                layer=layer
            )
            _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(LINKER_BEAM_SIZE, LINKER_NOTCH_SIZE, xmax - LINKER_EDGE_OFFSET / 2.0, _tmp_notch_ypos, layer)
            _tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R, _tmp_upper_triangle, 'not', max_points=0, layer=layer)
            _tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R, _tmp_lower_triangle, 'not', max_points=0, layer=layer)

            _tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer, _tmp_beam_polygon_L, 'or', max_points=0, layer=layer)
            _tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer, _tmp_beam_polygon_R, 'or', max_points=0, layer=layer)

        _tmp_notch_xarray = (_linkerdata_outer.dx - CORNER_BEND_RAD * 2 - LINKER_BEAM_SIZE) * numpy.linspace(0, 1.0, num=LINKER_X_NOTCHES)
        for i in range(LINKER_X_NOTCHES):
            _tmp_notch_xpos = xmin + LINKER_EDGE_OFFSET + CORNER_BEND_RAD + LINKER_BEAM_SIZE / 2.0 + _tmp_notch_xarray[i]
            _tmp_beam_polygon_T = gdspy.Polygon(
                [
                    (_tmp_notch_xpos - LINKER_BEAM_SIZE / 2.0, ymax),
                    (_tmp_notch_xpos - LINKER_BEAM_SIZE / 2.0, ymax - LINKER_EDGE_OFFSET),
                    (_tmp_notch_xpos + LINKER_BEAM_SIZE / 2.0, ymax - LINKER_EDGE_OFFSET),
                    (_tmp_notch_xpos + LINKER_BEAM_SIZE / 2.0, ymax),
                    (_tmp_notch_xpos - LINKER_BEAM_SIZE / 2.0, ymax)
                ],
                layer=layer
            )
            _tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(LINKER_BEAM_SIZE, LINKER_NOTCH_SIZE, _tmp_notch_xpos, ymax - LINKER_EDGE_OFFSET / 2.0, layer)
            _tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T, _tmp_left_triangle, 'not', max_points=0, layer=layer)
            _tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T, _tmp_right_triangle, 'not', max_points=0, layer=layer)

            _tmp_beam_polygon_B = gdspy.Polygon(
                [
                    (_tmp_notch_xpos - LINKER_BEAM_SIZE / 2.0, ymin),
                    (_tmp_notch_xpos - LINKER_BEAM_SIZE / 2.0, ymin + LINKER_EDGE_OFFSET),
                    (_tmp_notch_xpos + LINKER_BEAM_SIZE / 2.0, ymin + LINKER_EDGE_OFFSET),
                    (_tmp_notch_xpos + LINKER_BEAM_SIZE / 2.0, ymin),
                    (_tmp_notch_xpos - LINKER_BEAM_SIZE / 2.0, ymin)
                ],
                layer=layer
            )
            _tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(LINKER_BEAM_SIZE, LINKER_NOTCH_SIZE, _tmp_notch_xpos, ymin + LINKER_EDGE_OFFSET / 2.0, layer)
            _tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B, _tmp_left_triangle, 'not', max_points=0, layer=layer)
            _tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B, _tmp_right_triangle, 'not', max_points=0, layer=layer)

            _tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer, _tmp_beam_polygon_T, 'or', max_points=0, layer=layer)
            _tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer, _tmp_beam_polygon_B, 'or', max_points=0, layer=layer)

        _linker_combined = gdspy.fast_boolean(_tmp_linker_outer, _tmp_linker_inner, 'not', max_points=0, layer=layer)

        return _linker_combined

    # Parameters for alignment marks
    MARK_THICKNESS = 1e3
    MARK_LENGTH = 5e3

    @classmethod
    def write_alignment_mark(cls, cell, alignment_xpos, alignment_ypos, layer=3):
        cell.add(
            gdspy.Polygon(
                [
                    (alignment_xpos - cls.MARK_LENGTH / 2.0, alignment_ypos - cls.MARK_THICKNESS / 2.0),
                    (alignment_xpos - cls.MARK_LENGTH / 2.0, alignment_ypos + cls.MARK_THICKNESS / 2.0),
                    (alignment_xpos + cls.MARK_LENGTH / 2.0, alignment_ypos + cls.MARK_THICKNESS / 2.0),
                    (alignment_xpos + cls.MARK_LENGTH / 2.0, alignment_ypos - cls.MARK_THICKNESS / 2.0),
                    (alignment_xpos - cls.MARK_LENGTH / 2.0, alignment_ypos - cls.MARK_THICKNESS / 2.0),
                ], 
                layer=layer
            )
        )
        cell.add(
            gdspy.Polygon(
                [
                    (alignment_xpos - cls.MARK_THICKNESS / 2.0, alignment_ypos - cls.MARK_LENGTH / 2.0),
                    (alignment_xpos - cls.MARK_THICKNESS / 2.0, alignment_ypos + cls.MARK_LENGTH / 2.0),
                    (alignment_xpos + cls.MARK_THICKNESS / 2.0, alignment_ypos + cls.MARK_LENGTH / 2.0),
                    (alignment_xpos + cls.MARK_THICKNESS / 2.0, alignment_ypos - cls.MARK_LENGTH / 2.0),
                    (alignment_xpos - cls.MARK_THICKNESS / 2.0, alignment_ypos - cls.MARK_LENGTH / 2.0),
                ], 
                layer=layer
            )
        )

    # def get_beam_surface():
    #     return gdspy.Polygon(
    #         [
    #             (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
    #             (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy / 2.0),
    #             (_tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy / 2.0),
    #             (_tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
    #         ], 
    #         layer=layer
    #     )

    TEXT_DIST_TO_TOP = 6e3
    # Underdose Regions
    EDGE_UNDERDOSE = False

    @classmethod
    def get_grating_couplers(cls, beam_position, beam_dy_list):
        beam_position = copy.copy(beam_position)
        beam_position.ypos = beam_position.ypos - beam_dy_list[0] / 2.0
        grating_couplers = []
        for i in range(PHC_GROUP_COUNT):
            beam_position.ypos = beam_position.ypos + beam_dy_list[i] / 2.0
            beam_position.dy = beam_dy_list[i]
            grating_couplers.append(
                GratingCoupler(
                    beam_position=beam_position,
                    layer=5,
                )
            )
            beam_position.ypos = beam_position.ypos + beam_dy_list[i] / 2.0 + BEAM_SPACING
        return grating_couplers

    @classmethod
    def get_box_write_area(cls, _outer_box, _tmp_beam_polygon_combined, layer, write_linker, circ_grating, pattern_number, beamdata, _xmin, _xmax, _ymin, _ymax, round_corner, beam_dy_list, cell):
        _box_write_area = gdspy.fast_boolean(_outer_box, _tmp_beam_polygon_combined, 'not', max_points=0, layer=layer)
        if write_linker:
            _linker_combined = cls.write_linker_region(beamdata, _xmin, _xmax, _ymin, _ymax, round_corner, layer=layer)
            _box_write_area = gdspy.fast_boolean(_box_write_area, _linker_combined, 'not', max_points=0, layer=layer)
        if circ_grating:
            grating_couplers = cls.get_grating_couplers(beamdata, beam_dy_list)
            all_grating_coupler_geometry = Polygon.fast_combine(*([gc.left_shockwave for gc in grating_couplers] + [gc.right_shockwave for gc in grating_couplers]), layer=layer)
            _box_write_area = gdspy.fast_boolean(_box_write_area, all_grating_coupler_geometry, 'or', max_points=0, layer=layer)
        if pattern_number is not None:
            _left_pattern_number = create_phc_label(cell, _xmin + LINKER_EDGE_OFFSET + LINKER_WIDTH / 2.0, _ymax - LINKER_EDGE_OFFSET - cls.TEXT_DIST_TO_TOP, pattern_number)
            _right_pattern_number = create_phc_label(cell, _xmax - LINKER_EDGE_OFFSET - LINKER_WIDTH / 2.0, _ymax - LINKER_EDGE_OFFSET - cls.TEXT_DIST_TO_TOP, pattern_number)
            _pattern_number_combined = gdspy.fast_boolean(_left_pattern_number, _right_pattern_number, 'or', max_points=0, layer=layer)
            _box_write_area = gdspy.fast_boolean(_box_write_area, _pattern_number_combined, 'xor', max_points=0, layer=layer)
        return _box_write_area

    @classmethod
    def write_outer_box(
            cls, cell, beamdata, beam_dy_list, grating_spacer=False, round_corner=False, direct_write_area=False, alignment_mark=False,
            write_linker=False, pinch_pt_L=False, pinch_pt_R=False, pinch_pt_L_offset=0, pinch_pt_R_offset=0, pinch_pt_size=0, 
            circ_grating=False, pattern_number=None, reverse_tone=False, layer=3
        ):
        _xmin = beamdata.xpos - beamdata.dx / 2.0 - int(write_linker) * (LINKER_WIDTH + LINKER_EDGE_OFFSET)
        _xmax = beamdata.xpos + beamdata.dx / 2.0 + int(write_linker) * (LINKER_WIDTH + LINKER_EDGE_OFFSET)
        _ymin = beamdata.ypos - beam_dy_list[0] / 2.0 - EDGE_OFFSET
        _ymax = _ymin + (PHC_GROUP_COUNT - 1) * BEAM_SPACING + numpy.sum(beam_dy_list) + EDGE_OFFSET * 2
        _outer_box = gdspy.Polygon([(_xmin, _ymin), (_xmin, _ymax), (_xmax, _ymax), (_xmax, _ymin), (_xmin, _ymin)], layer=layer)
        if round_corner is True:
            _outer_box = _outer_box.fillet(CORNER_BEND_RAD, points_per_2pi=CORNER_BEND_PTS)

        if direct_write_area is True:
            _tmp_beamdata = copy.copy(beamdata)
            _tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0] / 2.0
            for i in range(PHC_GROUP_COUNT - 5):
                _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0
                if cls.EDGE_UNDERDOSE and not reverse_tone:
                    _tmp_beamdata.dy = beam_dy_list[i] + UNDER_EDGE_SIZE * 2
                else:
                    _tmp_beamdata.dy = beam_dy_list[i]

                beam = NanoBeam(
                    [
                        (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
                        (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy / 2.0),
                        (_tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy / 2.0),
                        (_tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
                        (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
                    ], 
                    layer=layer
                )
                


                _tmp_beam_polygon = gdspy.Polygon(
                    [
                        (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
                        (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy / 2.0),
                        (_tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy / 2.0),
                        (_tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
                        (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
                    ], 
                    layer=layer
                )
                if pinch_pt_L:
                    _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(
                        _tmp_beamdata.dy, pinch_pt_size, 
                        _tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0 + pinch_pt_L_offset, _tmp_beamdata.ypos, layer
                    )
                    _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _tmp_upper_triangle, 'not', max_points=0, layer=layer)
                    _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _tmp_lower_triangle, 'not', max_points=0, layer=layer)
                if pinch_pt_R:
                    _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(
                        _tmp_beamdata.dy, pinch_pt_size, 
                        _tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0 - pinch_pt_R_offset, 
                        _tmp_beamdata.ypos, layer
                    )
                    _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _tmp_upper_triangle, 'not', max_points=0, layer=layer)
                    _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _tmp_lower_triangle, 'not', max_points=0, layer=layer)

                _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0 + BEAM_SPACING
                if i == 0:
                    _tmp_beam_polygon_combined = _tmp_beam_polygon
                else:
                    _tmp_beam_polygon_combined = gdspy.fast_boolean(_tmp_beam_polygon_combined, _tmp_beam_polygon, 'or', max_points=0, layer=layer)
            _box_write_area = cls.get_box_write_area(_outer_box, _tmp_beam_polygon_combined, layer, write_linker, circ_grating, pattern_number, beamdata, _xmin, _xmax, _ymin, _ymax, round_corner, beam_dy_list, cell)

            rough = gdspy.Rectangle(
                (beamdata.xpos - beamdata.dx / 2, beamdata.ypos - EDGE_OFFSET - beamdata.dy / 2), 
                (beamdata.xpos - beamdata.dx / 2 + SPACER, beamdata.ypos + EDGE_OFFSET + beamdata.dy / 2),
                layer=layer
            )
            rough2 = gdspy.Rectangle(
                (beamdata.xpos + beamdata.dx / 2, beamdata.ypos - EDGE_OFFSET - beamdata.dy / 2),
                (beamdata.xpos + beamdata.dx / 2 - SPACER, beamdata.ypos + EDGE_OFFSET + beamdata.dy / 2),
                layer=layer
            )
            rough2 = gdspy.fast_boolean(rough2, _tmp_beam_polygon_combined, 'not', max_points=0, layer=layer)
            rough = gdspy.fast_boolean(rough, _tmp_beam_polygon_combined, 'not', max_points=0, layer=layer)

        if reverse_tone:
            _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_combined, 'or', max_points=0, layer=layer)
            _box_write_area_reverse = gdspy.fast_boolean(_outer_box, _box_write_area, 'not', max_points=0, layer=layer)
            if grating_spacer:
                cell.add(rough)
                cell.add(rough2)
            cell.add(_box_write_area_reverse)
        else:
            if grating_spacer:
                _box_write_area = gdspy.fast_boolean(_box_write_area, rough, 'not', max_points=0, layer=layer)
                _box_write_area = gdspy.fast_boolean(_box_write_area, rough2, 'not', max_points=0, layer=layer)
            cell.add(_box_write_area)

        if alignment_mark:
            _yoffset = 10000
            _alignment_x = beamdata.xpos
            _alignment_y = _ymax + _yoffset
            cls.write_alignment_mark(cell, _alignment_x, _alignment_y)

        return _outer_box

    def write_beam_array(self, cell, first_beam, width_list=None):

        if width_list is None:
            _dy_list = first_beam.dy * numpy.ones(PHC_GROUP_COUNT, dtype=int)
        else:
            _dy_list = width_list

        _beam_write = copy.copy(first_beam)
        _initial_ypos = first_beam.ypos - _dy_list[0] / 2.0
        for i in range(PHC_GROUP_COUNT):
            _beam_write.dy = _dy_list[i]
            _beam_write.ypos = _initial_ypos + i * BEAM_SPACING + numpy.sum(_dy_list[:i]) + _beam_write.dy / 2
        # write_beam_single(cell, _beam_write)
        return _dy_list

    # Overdose Regions
    RING_OVERDOSE = False
    EDGE_OVERDOSE = False

    # Underdose Regions
    RING_UNDERDOSE = False

    def create_phc_group(self, param_sweep, end_period, blank_guides, create_label=True):

        acav = param_sweep[1][0]
        amir = param_sweep[1][1]
        hx = param_sweep[1][5]
        hy = param_sweep[1][6]
        wy = param_sweep[1][2]
        num_cav = param_sweep[1][4]
        num_tap = param_sweep[1][7]
        num_mirr = param_sweep[1][3] * 2
        beams = gdspy.Cell(str(param_sweep[0]))

        beamdata_rowicolj = self.PLTdata(xpos=0, ypos=0, dx=10e3, dy=wy)
        beam_dy_list = self.write_beam_array(beams, beamdata_rowicolj)

        holedata_circ_rowicolj = self.PLTdata(xpos=beamdata_rowicolj.xpos - HOLE_POS_OFFSET, ypos=beamdata_rowicolj.ypos, dx=hx, dy=hy)

        _initial_ypos = holedata_circ_rowicolj.ypos - beam_dy_list[0] / 2.0
        _tmp_beamdata = copy.copy(beamdata_rowicolj)
        _tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0] / 2.0

        for i in range(PHC_GROUP_COUNT):
            NanoBeam.create_nano_beam(
                i              = i,
                _initial_ypos  = _initial_ypos,
                _tmp_beamdata  = _tmp_beamdata,
                cell           = beams, 
                holedata       = holedata_circ_rowicolj, 
                beam_dy_list   = beam_dy_list, 
                num_cav_hole   = num_cav,
                num_mir_hole_L = num_mirr / 2,
                num_mir_hole_R = num_mirr / 2, 
                acav           = acav, 
                amir           = amir, 
                end_period     = end_period, 
                guides         = blank_guides,
                end_taper_L    = True, 
                end_taper_R    = True, 
                num_end_taper  = num_tap, 
                reverse_tone   = self.REVERSE_TONE, 
                edge_overdose  = self.EDGE_OVERDOSE,
                ring_overdose  = self.RING_OVERDOSE,
                ring_underdose = self.RING_UNDERDOSE, 
                edge_underdose = self.RING_UNDERDOSE
            )

        direct_write_area = True
        circ_grating = True

        outer_box = self.write_outer_box(
            cell                 = beams,  # this cell seems to have all of the grating coupler references and the bounding box
            # gdspy.Cell(f"temptemptemp-{param_sweep[0]}"), 
            beamdata             = beamdata_rowicolj, 
            beam_dy_list         = beam_dy_list, 
            grating_spacer       = self.GRATING_SPACER, 
            round_corner         = False, 
            direct_write_area    = direct_write_area, 
            alignment_mark       = False,
            write_linker         = True, 
            pinch_pt_L           = False, 
            pinch_pt_R           = False, 
            pinch_pt_L_offset    = 0, 
            pinch_pt_R_offset    = 0, 
            pinch_pt_size        = 0, 
            circ_grating         = circ_grating, 
            pattern_number       = None, 
            reverse_tone         = self.REVERSE_TONE, 
            layer                = 3
        )

        SupportStructure.create_outer_frame_with_supports(
            beams,
            beamdata_rowicolj, 
            beam_dy_list,
            outer_box.get_bounding_box(), 
            pattern_number=param_sweep[0] if create_label else None, 
            reverse_tone=self.REVERSE_TONE
        )
        return beams

    def supportsMask(self, name, spacing):
        x_frame = spacing*6 * Unit.um.value /2
        y_frame = ((LINKER_WIDTH/1e3)*2+self.WY_LIST[0]/1e3) * Unit.um.value * 12/2

        phc_frame = gdspy.Rectangle((-x_frame, -y_frame), (x_frame, y_frame), layer=4)

        support_mask_inner = 0.82
        spacing = 440
        support_mask_out = gdspy.Rectangle(
            (-spacing * support_mask_inner * Unit.um.value / 2, spacing * support_mask_inner * Unit.um.value / 2),
            (spacing * support_mask_inner * Unit.um.value / 2, -spacing * support_mask_inner * Unit.um.value / 2)
        )

        holder = gdspy.Cell('support' + str(name))

        support_mask_frame = gdspy.boolean(support_mask_out, phc_frame, "not", layer=6)
        holder.add(support_mask_frame)

        support_mask_pos = gdspy.CellReference(holder, (x_frame/2, 3300))
        CELL_MAIN.add(support_mask_pos)

    # Taper for fibre and reflector
    GRATING_TAPER_LENGTH = 5
    REFLECTOR_TAPER_LENGTH = 5

    # Fibre-coupling grating
    GRATING_PERIOD = 2.05
    GRATING_FF = 0.5
    GRATING_W = 2
    N_GRATINGS = 8

    GRATING_FW = GRATING_PERIOD * GRATING_FF
    GRATING_EW = GRATING_PERIOD * (1 - GRATING_FF)

    # Reflector grating
    REFLECTOR_PERIOD = 1.2
    REFLECTOR_FF = 0.775
    REFLECTOR_W = 2
    REFLECTOR_N_GRATINGS = 12

    REFLECTOR_FW = REFLECTOR_PERIOD * REFLECTOR_FF
    REFLECTOR_EW = REFLECTOR_PERIOD * (1 - REFLECTOR_FF)

    def defineDevice(self, wg_w, ring_r, coupling_l, coupling_gap, grating_offset, g):
        # Create the grating couplers
        for i in range(self.N_GRATINGS):
            grating_post = gdspy.Rectangle((self.GRATING_EW + i * self.GRATING_PERIOD, - self.GRATING_W / 2), (self.GRATING_EW + i * self.GRATING_PERIOD + self.GRATING_FW, self.GRATING_W / 2))
            GRATING_CELL.add(grating_post)

        # Create the reflector grating
        for i in range(self.REFLECTOR_N_GRATINGS):
            reflector_post = gdspy.Rectangle(
                (self.REFLECTOR_EW + i * self.REFLECTOR_PERIOD, - self.REFLECTOR_W / 2),
                (self.REFLECTOR_EW + i * self.REFLECTOR_PERIOD + self.REFLECTOR_FW, self.REFLECTOR_W / 2)
            )
            REFLECTOR_CELL.add(reflector_post)

        # Create the grating coupler taper
        grating_taper_points = [(0, -wg_w / 2), (0, wg_w / 2), (self.GRATING_TAPER_LENGTH, self.GRATING_W / 2), (self.GRATING_TAPER_LENGTH, -self.GRATING_W / 2)]
        grating_taper_poly = gdspy.Polygon(grating_taper_points)
        GRATING_TAPER_CELL.add(grating_taper_poly)

        # Create the reflector taper
        reflector_taper_points = [(0, -wg_w / 2), (0, wg_w / 2), (self.REFLECTOR_TAPER_LENGTH, self.REFLECTOR_W / 2), (self.REFLECTOR_TAPER_LENGTH, -self.REFLECTOR_W / 2)]
        reflector_taper_poly = gdspy.Polygon(reflector_taper_points)
        REFLECTOR_TAPER_CELL.add(reflector_taper_poly)

        # Create the ring resonator
        ring = gdspy.Path(wg_w)
        ring.segment(coupling_l / 2, "+x")
        ring.arc(ring_r, math.pi / 2, -math.pi / 2)
        ring.segment(coupling_l / 2, "-x")

        ring2 = gdspy.copy(ring)
        ring2.mirror((0, 0), (0, 100))

        # Create the coupling region
        coupler = gdspy.Rectangle((-coupling_l / 2 - grating_offset, coupling_gap + wg_w / 2), (coupling_l / 2 + grating_offset, coupling_gap + wg_w / 2 + wg_w))
        grating_taper = gdspy.CellReference(GRATING_TAPER_CELL, (-coupling_l / 2 - grating_offset, wg_w + coupling_gap), 180)
        coupler_grating = gdspy.CellReference(GRATING_CELL, (-coupling_l / 2 - grating_offset - self.GRATING_TAPER_LENGTH, wg_w + coupling_gap), 180)
        reflector_taper = gdspy.CellReference(REFLECTOR_TAPER_CELL, (coupling_l / 2 + grating_offset, wg_w + coupling_gap))
        reflector_grating = gdspy.CellReference(REFLECTOR_CELL, (coupling_l / 2 + grating_offset + self.REFLECTOR_TAPER_LENGTH, wg_w + coupling_gap))

        # Add a grating-only region
        # grating_test_wg = gdspy.Rectangle((-coupling_l/2 - grating_offset,-2*ring_r-5), (coupling_l/2 + grating_offset,-2*ring_r-5-wg_w))
        # grating_test_in = gdspy.CellReference(grating_cell,(-coupling_l/2 - grating_offset, -2*ring_r-5-wg_w/2),180)
        # grating_test_out = gdspy.CellReference(grating_cell,(coupling_l/2 + grating_offset, -2*ring_r-5-wg_w/2),0)

        device = gdspy.Cell(str(g) + str(wg_w) + str(coupling_gap) + 'device')
        device.add(gdspy.Text(str(int(wg_w * 1000)), 4, (-2.5, -30)))
        device.add(gdspy.Text(str(int(coupling_gap * 1000)), 4, (-2.5, -20)))
        device.add([ring, ring2, coupler, grating_taper, coupler_grating, reflector_taper, reflector_grating])
        return device
    
    # Define the cell reference for the support structure
    CELL_SUPPORT = gdspy.Cell('support_structure')  # Single instance of the support structure for the membrane layer
    SUPPORT_VERTICES = [
        (-GAP / 2, SUPPORT_WIDTH / 2), 
        (-SUPPORT_PINCH_LENGTH / 2, SUPPORT_PINCH_WIDTH / 2),
        (0, SUPPORT_PINCH_WIDTH / 2), (0, 0), (-GAP / 2, 0)
    ]
    SUPPORT_1 = gdspy.Polygon(SUPPORT_VERTICES, layer=5)
    SUPPORT_2 = gdspy.Polygon(SUPPORT_VERTICES, layer=5)
    SUPPORT_3 = gdspy.Polygon(SUPPORT_VERTICES, layer=5)
    SUPPORT_4 = gdspy.Polygon(SUPPORT_VERTICES, layer=5)
    SUPPORT_2.mirror((0, -1), (0, 1))
    SUPPORT_3.mirror((-1, 0), (1, 0))
    SUPPORT_4.mirror((0, -1), (0, 1))
    SUPPORT_4.mirror((-1, 0), (1, 0))
    CELL_SUPPORT.add(SUPPORT_1)
    CELL_SUPPORT.add(SUPPORT_2)
    CELL_SUPPORT.add(SUPPORT_3)
    CELL_SUPPORT.add(SUPPORT_4)

    def defineMembrane(self, identifier, spacing, spacing_y, length, height, layer = 5):
        # Create a membrane cell for each desired membrane size
        cell_membrane_nohole = gdspy.Cell('membrane' + identifier)
        # Define the frame for the membrane
        membrane_outer_frame = gdspy.Rectangle((-spacing / 2, spacing_y / 2), (spacing / 2, -spacing_y / 2), layer=layer)
        membrane_inner_frame = gdspy.Rectangle((-length / 2 - GAP, height / 2 + GAP), (length / 2 + GAP, -height / 2 - GAP), layer=layer)
        membrane_frame = gdspy.boolean(membrane_outer_frame, membrane_inner_frame, "not", layer=layer)
        # Define the membrane itself
        membrane = gdspy.Rectangle((-length / 2, height / 2), (length / 2, -height / 2), layer=layer)
        # Add to the membrane cell
        cell_membrane_nohole.add(membrane_frame)

        if self.PERFORATIONS:
            factor = 0.5
            perfPos = length / 2 * factor
            perfSize = 3
            perfPoints = [(1, -1), (1, 1), (-1, 1), (-1, -1)]
            for i in perfPoints:
                x = i[0] * perfPos
                y = i[1] * perfPos
                perf = gdspy.Rectangle((x, y), (x + perfSize, y + perfSize), layer=layer)
                membrane = gdspy.boolean(membrane, perf, "not", layer=layer)
            perf = gdspy.Rectangle((-perfSize / 2, -perfSize / 2), (perfSize / 2, perfSize / 2), layer=layer)
            perf_cut = gdspy.boolean(membrane, perf, "not", layer=layer)
            cell_membrane_nohole.add(perf_cut)
        else:
            cell_membrane_nohole.add(membrane)

        # Add the support structures
        support_period = SUPPORT_WIDTH / SUPPORT_PITCH
        air_width = support_period * (1 - SUPPORT_PITCH)

        num_supports_x = math.floor(length / support_period + SUPPORT_PITCH)
        num_supports_y = math.floor(height / support_period + SUPPORT_PITCH)

        if (num_supports_y % 2) == 0:
            # Number of supports is even
            support_length_y = num_supports_y * support_period - support_period / 2
            support_offset_y = (length - support_length_y) / 2 + SUPPORT_WIDTH / 2
            for i in range(0, num_supports_y):
                ref1 = gdspy.CellReference(self.CELL_SUPPORT, (
                    -length / 2 - GAP / 2, -length / 2 + support_offset_y + i * support_period))
                ref2 = gdspy.CellReference(self.CELL_SUPPORT, (
                    length / 2 + GAP / 2, -length / 2 + support_offset_y + i * support_period))
                cell_membrane_nohole.add([ref1, ref2])
        else:
            # Number of supports is odd
            for i in range(-math.floor(num_supports_y / 2), math.floor(num_supports_y / 2) + 1):
                ref1 = gdspy.CellReference(self.CELL_SUPPORT, (-length / 2 - GAP / 2, i * support_period))
                ref2 = gdspy.CellReference(self.CELL_SUPPORT, (length / 2 + GAP / 2, i * support_period))
                cell_membrane_nohole.add([ref1, ref2])
        # Membrane Sides For Rectangles
        if (num_supports_x % 2) == 0:
            # For rectangles
            support_length_x = num_supports_x * support_period - support_period / 2
            support_offset_x = (height - support_length_x) / 2 + SUPPORT_WIDTH / 2
            for i in range(0, num_supports_x):
                # Rectangle
                ref3 = gdspy.CellReference(self.CELL_SUPPORT, (-height / 2 + support_offset_x + i * support_period, height / 2 + GAP / 2), rotation=90)
                ref4 = gdspy.CellReference(self.CELL_SUPPORT, (-height / 2 + support_offset_x + i * support_period, -height / 2 - GAP / 2), rotation=90)
                cell_membrane_nohole.add([ref3, ref4])
        else:
            # Number of supports is odd
            for i in range(-math.floor(num_supports_x / 2), math.floor(num_supports_x / 2) + 1):
                ref3 = gdspy.CellReference(self.CELL_SUPPORT, (i * support_period, height / 2 + GAP / 2), rotation=90)
                ref4 = gdspy.CellReference(self.CELL_SUPPORT, (i * support_period, -height / 2 - GAP / 2), rotation=90)
                cell_membrane_nohole.add([ref3, ref4])

        cell_membrane_ref_nohole = gdspy.CellReference(cell_membrane_nohole, (0, 0))
        self.CELL_MEMBRANE_ROW_NOHOLE.add(cell_membrane_ref_nohole)
        return self.CELL_MEMBRANE_ROW_NOHOLE
