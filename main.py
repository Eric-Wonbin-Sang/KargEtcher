# def write_hole_1D_mir_sweep(cell, beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, amir,
#                             end_taper_L=False, end_taper_R=False, num_end_taper=0, reverse_tone=False, layer=2):
#     _one_side_holes = 0
#     _amir_list_L = [amir for x in range(int(num_mir_hole_L + num_end_taper * end_taper_L))]
#     _amir_list_R = [amir for x in range(int(num_mir_hole_R + num_end_taper * end_taper_R))]
#     _aper_list = copy.copy(_amir_list_L)
#     _aper_list.extend(_amir_list_R)
#     if len(_aper_list) > 0:
#         _hole_write = copy.copy(holedata)
#         _hole_write.xpos = _hole_write.xpos - numpy.sum(numpy.array(_amir_list_L))
#         _taper_scale_list = []
#         if num_end_taper > 0:
#             _taper_scale_list = numpy.linspace(0, 1.0, num_end_taper + 2)
#             _taper_scale_list_L = _taper_scale_list[1:-1]
#             _taper_scale_list_R = numpy.flipud(_taper_scale_list[1:-1])

#         for i in range(len(_aper_list)):
#             _hole_write.xpos = _hole_write.xpos + _aper_list[i] / 2.0
#             if i < num_end_taper * end_taper_L:
#                 _hole_write.dx = holedata.dx * _taper_scale_list_L[i]
#                 _hole_write.dy = holedata.dy * _taper_scale_list_L[i]
#             else:
#                 _hole_write.dx = holedata.dx
#                 _hole_write.dy = holedata.dy
#             if i >= (len(_aper_list) - num_end_taper * end_taper_R):
#                 _hole_write.dx = holedata.dx * _taper_scale_list_R[i - len(_aper_list) + num_end_taper]
#                 _hole_write.dy = holedata.dy * _taper_scale_list_R[i - len(_aper_list) + num_end_taper]
#             _single_hole_polygon = write_hole_single(cell, _hole_write)

#             if i == 0:
#                 _hole_polygon_combined = _single_hole_polygon
#             else:
#                 _hole_polygon_combined = gdspy.fast_boolean(_hole_polygon_combined, _single_hole_polygon, 'or',
#                                                             max_points=0, layer=layer)

#             _hole_write.xpos = _hole_write.xpos + _aper_list[i] / 2.0

#         if reverse_tone is True:
#             _tmp_beam_polygon = gdspy.Polygon([
#                 (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
#                 (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
#                 (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
#                 (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
#                 (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0)], layer=layer)
#             _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _hole_polygon_combined, 'not', max_points=0,
#                                                    layer=layer)
#             cell.add(_tmp_beam_polygon)
#         else:
#             cell.add(_hole_polygon_combined)
#     else:
#         _tmp_beam_polygon = gdspy.Polygon([
#             (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
#             (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
#             (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
#             (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
#             (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0)], layer=layer)
#         cell.add(_tmp_beam_polygon)


# def write_hole_2D_mir_sweep(cell, beamdata, holedata, beam_dy_list, num_cav_hole, num_mir_hole_L, num_mir_hole_R,
#                             acav, amir, guides,
#                             end_taper_L=False, end_taper_R=False, num_end_taper=0, reverse_tone=False,
#                             edge_overdose=False, ring_overdose=False):
#     _initial_ypos = holedata.ypos - beam_dy_list[0] / 2.0
#     _tmp_beamdata = copy.copy(beamdata)
#     _tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0] / 2.0

#     for num_mir in range(num_beam):
#         holedata.ypos = _initial_ypos + num_mir * beam_spacing + numpy.sum(beam_dy_list[:num_mir]) + beam_dy_list[num_mir] / 2.0
#         _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[num_mir] / 2.0
#         _tmp_beamdata.dy = beam_dy_list[num_mir]
#         write_hole_1D_mir_sweep(cell, _tmp_beamdata, holedata, num_cav_hole, num_mir, num_mir, acav, amir, end_taper_L, end_taper_R, num_end_taper, reverse_tone=reverse_tone)
#         _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[num_mir] / 2.0 + beam_spacing

from datetime import datetime
from functools import cached_property
import os
import gdspy
from classes.base_classes import Component

from classes.photonic_components import WorkArea, Unit, Hole, Beam, GratingCoupler, Phc, logger


class Creator:
    
    """
    EBPG: Electron Beam Projection Lithography (Raith EBPG system)
    Heidelberg Laser Writer: high-precision laser lithography system

    https://gdsfactory.github.io/gdsfactory/index.html
    https://en.wikipedia.org/wiki/Avalanche_photodiode

    A cell that is referenced must have the referenced geometry exist on some cell.
    """

    SAVE_DIR = "C:\\Users\\ericw\\playground\\gds_files"
    UNIT = 1.0e-9
    PRESICION = 1.0e-11

    def __init__(self, name) -> None:

        self.name = name
        
        self.work_area = WorkArea(
            x_length=6 * Unit.mm.value,
            y_length=6 * Unit.mm.value,
            layer=4,
            add_position_ids=True,
            add_litho_marks=True,
        )

        self.phc = Phc(
            beam=Beam(width=60 * Unit.um.value, height= 30 * Unit.um.value),
            grating_coupler=GratingCoupler(),
        )

    @cached_property
    def filename(self):
        """ Ex: 2023-11-18-M12.gds """
        return f"{datetime.now().strftime('%Y-%m-%d')}-{self.name}.gds"

    @cached_property
    def filepath(self):
        return os.path.join(self.SAVE_DIR, self.filename)

    @property
    def all_components(self):
        return [
            self.work_area
        ]

    def preprocess_geometries(self):
        Component.add_components_to_god_cell(self.all_components)

    def save_file(self, unit=UNIT, precision=PRESICION):
        self.preprocess_geometries()
        gdspy.write_gds(self.filepath, unit=unit, precision=precision)
        if os.path.exists(self.filepath):
            logger.info(f"{self.filepath} saved!")
        else:
            logger.error(f"{self.filepath} could not be saved!")


def main():
    creator = Creator(name="test")
    creator.save_file()


if __name__ == "__main__":
    main()
