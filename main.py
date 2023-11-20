# Import requisite python packages
import gdspy
import numpy
import math
import copy

# Define the fixed variables of the structure
import numpy as np


nm = 1
um = 1000

# General Pattern Layout
writeMems = False
writePhC = True
PhCSweeper = not writePhC
writeQFC = False
EBPG_markers = True
litho_markers = True
supportsMaskBool = False
removeInnerMem = True


# Define the key parameters of the membrane and support structure
chip_size = 6000 * um  # The size of our chip - assumed to be 6mm x 6mm - reduce slightly in case of variations

invert = False

# Create the cells of the cad file
cell_main = gdspy.Cell("main_cell")
grating_cell = gdspy.Cell("grating")
reflector_cell = gdspy.Cell("reflector")
grating_taper_cell = gdspy.Cell("grating_taper")
reflector_taper_cell = gdspy.Cell("reflector_taper")
device_pos = gdspy.Cell('device_pos')
support_mask = gdspy.Cell('support_mask')
phc_frame = gdspy.Cell('phc_frame')
phc_sweep = gdspy.Cell('phc_sweep')
EBPG_marks = gdspy.Cell("ebpg_marks")
heid_alignment_mark = gdspy.Cell('heid_alignment mark')
phc_sweep_frame_cell = gdspy.Cell('phc_sweep_frame')

devs = 6

num_rows = 6# Number of times to write each row of membranes
num_cols = 6

dev_list = np.arange(0, num_rows * num_cols, 1)



# Unit for all parameters: nm
# Parameters for exposure
write_field_size = 100e3
write_field_x_size = 50e3
write_field_y_size = 100e3
shot_pitch = 0.1
pattern_x_sep = 80e3
pattern_y_sep = 115e3

# Parameters for nanobeam
num_beam = 11
pinch_pt_size = [100]
beam_spacing = 6.5e3
num_circ_points = 181
hole_pos_offset = 0
edge_offset = 5e3
corner_bend_rad = 3e3
corner_bend_pts = 41

num_rows_phc = 12
num_cols_phc = 12
spacing_phc = 2*(edge_offset/1e3)+10

acav = [155.5]
amir = [172.1]

wy_list = [470]
wy_end_taper_list = []

hx_list = [63, 66, 69, 73, 76, 79]
hy_list = [180, 184, 188, 192, 196, 200]

mirror_list = [10]
cavity_list = [16]
taper_holes = [0]


phc_dev_list = np.arange(0, num_rows_phc * num_cols_phc, 1)
param_sweep = np.array(np.meshgrid(acav, amir, wy_list, mirror_list, cavity_list, hx_list, hy_list, taper_holes)).T.reshape(-1, 8).astype(int)
param_sweep = list(zip(dev_list, param_sweep))


custom_cavity = True
seg = np.array([1.57385141e-07, 1.58156287e-07, 1.59014195e-07, 1.59789244e-07,
                1.61512865e-07, 1.63789232e-07, 1.65369591e-07, 1.67000125e-07]) * 1e9
custom_cavity_list = np.append(np.flip(seg), seg)
num_guides = 0  # Define the number of blank waveguides for control measurements (0,1,2)

taper_neff = 1.8
end_period = np.round(955 / (3.1 * 2), 1)

# For grating spacer
grating_spacer = False
spacer = 0  # Spacer size

# Tone
reverse_tone = True
# Overdose Regions
ring_overdose = False
edge_overdose = False

over_ring_size = 12
over_edge_size = 12

# Underdose Regions
ring_underdose = False
edge_underdose = False

under_ring_size = 12
under_edge_size = 12

# Parameters for alignment marks
mark_thickness = 1e3
mark_length = 5e3

# Parameters for the linker
linker_edgeoffset = 0
linker_width = 6e3
linker_connector_width = 1e3
linker_notch_size = 500
linker_beam_size = 500
linker_xnotches = 0
linker_ynotches = 0

# Parameter for the support
support_connector_width = 1e3 #0.000001e3
support_notch_size = 400
support_beam_size = 1e3
num_xsupport = 2
num_ysupport = 10
overlap_width = 0

# Parameters for circular grating structure
period = 777e-9
duty_cycle = 0.376
grating_linewidth = np.round(period * (1 - duty_cycle) / 1e-9, 1)
grating_spacing = np.round(period * (duty_cycle) / 1e-9, 1)

circ_grating_support = [200]
num_circ_grating_points = 181

# grating_spacing = 280
# grating_linewidth = 160  # design width = 160nm

num_grating = 2
circ_grating_base = 1500
grating_angle = 10 * numpy.pi / 180  # 15 degree in radian
odd_support_angle = 75 * numpy.pi / 180  # 55 degree in radian
even_support_angle = numpy.pi / 2.0  # 90 degree in radian

support_angle_width = 2.3 * numpy.pi / 180  # 10 degree in radian

# Parameters for the text
text = True
textheight = 8e3
textwidth = textheight * 20.0 / 9.0
textsep = textheight * 8.0 / 9.0
text_dist_to_top = 6e3

matrix_x_size = len(hy_list)
matrix_y_size = len(mirror_list)
block_x_size = len(hx_list)
block_y_size = len(cavity_list)
block_x_sep = pattern_x_sep * matrix_x_size
block_y_sep = (pattern_y_sep) * matrix_y_size
k = 0

origin_x = 0
origin_y = 0

# Membrane Parameters

gap = 10  # 10um gap between the membrane and support structure
support_width = 10  # The width of the support for the suspended layer
support_pinch_width = 3  # The width of the pinch point for stamping
support_pinch_length = 2  # The length of the pinch point region
support_pitch = 0.4  # Ratio of support structure to exposed area

# Membrane parameters
perforations = False

# length = spacing_phc*12
# height = ((linker_width / 1e3) * 2 + wy_list[0] / 1e3)*12

length = 150
height = 200


safety = 60 # Space between membranes
spacing = length + safety + support_width * 2
spacing_y = height + safety + support_width * 2


phc_y_offset = (35 - 0.275) * um  # An arbitrary parameter for centering the devices at the origin

column_size = num_rows * spacing
pairs = 1
# PhC Device Parameters

# QFC Device Parameters
widths = np.linspace(0, 1, devs)
gaps = np.linspace(0, 1, devs)

wg_w = 0.36

ring_r = 25
coupling_l = 25
coupling_gap = 0.15
grating_offset = 10

ringPos_offset_y = 20

# Fibre-coupling grating
grating_period = 2.05
grating_ff = 0.5
grating_w = 2
n_gratings = 8

grating_fw = grating_period * grating_ff
grating_ew = grating_period * (1 - grating_ff)

# Reflector grating
reflector_period = 1.2
reflector_ff = 0.775
reflector_w = 2
reflector_n_gratings = 12

reflector_fw = reflector_period * reflector_ff
reflector_ew = reflector_period * (1 - reflector_ff)

# Taper for fibre and reflector
grating_taper_length = 5
reflector_taper_length = 5

# Define Alignment Marks

write_size = (length + safety + gap * 2) * devs
mark_pos = write_size * um

outer_frame = gdspy.Rectangle((-chip_size / 2, chip_size / 2), (chip_size / 2, -chip_size / 2), layer=4)
inner_frame = gdspy.Rectangle((-spacing * um * num_rows / 2, spacing_y * um * num_cols * pairs / 2),
                              (spacing * um * num_rows / 2, -spacing_y * um * num_cols * pairs / 2), layer=4)
frame = gdspy.boolean(outer_frame, inner_frame, "not", layer=4)




def executeWriting():

    for k in range(num_rows_phc):

        # Write multiple rows
        for j in range(num_cols_phc):
            identifier = str(j) + str(k)

            if PhCSweeper == True:
                xpos = (k * spacing_phc - spacing_phc * (num_rows_phc - 1) / 2) * um
                ypos = ((linker_width / 1e3) * 2 + wy_list[0] / 1e3) / spacing_phc * (
                            spacing_phc / 2 + j * spacing_phc - spacing_phc * num_cols_phc * pairs / 2) * um


                name = num_cols_phc * k + j
                phc = PhC_Writer(param_sweep[name], end_period=end_period, blank_guides=num_guides, text=text)

                # phc_sweep_y_offset = -34.84
                phc_pos = gdspy.CellReference(phc, (xpos, ypos + phc_y_offset))
                phc_sweep.add(phc_pos)


    i=0

    for k in range(num_rows):

        # Write multiple rows
        for j in range(num_cols):
            identifier = str(j) + str(k)

            xpos = (k * spacing - spacing * (len(widths) - 1) / 2) * um
            ypos = (spacing_y / 2 + j * spacing_y - spacing_y * len(gaps) * pairs / 2) * um
            # QFC Device
            if writeQFC == True:
                device = defineDevice(widths[k], ring_r, coupling_l, gaps[j], grating_offset, g)
                device_pos = gdspy.CellReference(device, (xpos, ypos + ringPos_offset_y))
                cell_main.add(device_pos)

            if PhCSweeper == True and writePhC == False and i==0:
                name = num_cols * k + j

                phc_pos = gdspy.CellReference(phc_sweep, (xpos, ypos - phc_y_offset))
                cell_main.add(phc_pos)

            if writePhC == True:
                name = num_cols * k + j
                phc = PhC_Writer(param_sweep[name], end_period=end_period, blank_guides=num_guides, text=text)

                phc_pos = gdspy.CellReference(phc, (xpos, ypos - phc_y_offset))
                cell_main.add(phc_pos)

            if writeMems == True and i==0:
                membrane = defineMembrane(identifier, spacing, spacing_y, length, height, layer=5)

                phc_sweep_frame = gdspy.Rectangle((-length/2, -height/2), (length/2, height/2), layer=5)
                phc_sweep_frame_cell.add(phc_sweep_frame)
                phc_sweep_frame_cell_pos = gdspy.CellReference(phc_sweep_frame_cell, (xpos, ypos), magnification=1000)
                membrane_pos = gdspy.CellReference(membrane, (xpos, ypos), magnification=1000)

                if removeInnerMem == True:
                    membrane_pos = gdspy.boolean(membrane_pos, phc_sweep_frame_cell_pos, "not", layer=5)
                cell_main.add(membrane_pos)

            # # Supports Mask

            if supportsMaskBool == True:
                if k == 0 and j == 0:
                    supportsMask(name, xpos, ypos, spacing)

            if EBPG_markers == True:
                ebeam_marks(frame)
            if litho_markers == True:
                litho_marks()
            i+=1

def defineDevice(wg_w, ring_r, coupling_l, coupling_gap, grating_offset, g):
    # Create the grating couplers
    for i in range(n_gratings):
        grating_post = gdspy.Rectangle((grating_ew + i * grating_period, - grating_w / 2),
                                       (grating_ew + i * grating_period + grating_fw, grating_w / 2))
        grating_cell.add(grating_post)

    # Create the reflector grating
    for i in range(reflector_n_gratings):
        reflector_post = gdspy.Rectangle((reflector_ew + i * reflector_period, - reflector_w / 2),
                                         (reflector_ew + i * reflector_period + reflector_fw, reflector_w / 2))
        reflector_cell.add(reflector_post)

    # Create the grating coupler taper
    grating_taper_points = [(0, -wg_w / 2), (0, wg_w / 2), (grating_taper_length, grating_w / 2),
                            (grating_taper_length, -grating_w / 2)]
    grating_taper_poly = gdspy.Polygon(grating_taper_points)
    grating_taper_cell.add(grating_taper_poly)

    # Create the reflector taper
    reflector_taper_points = [(0, -wg_w / 2), (0, wg_w / 2), (reflector_taper_length, reflector_w / 2),
                              (reflector_taper_length, -reflector_w / 2)]
    reflector_taper_poly = gdspy.Polygon(reflector_taper_points)
    reflector_taper_cell.add(reflector_taper_poly)

    # Create the ring resonator
    ring = gdspy.Path(wg_w)
    ring.segment(coupling_l / 2, "+x")
    ring.arc(ring_r, math.pi / 2, -math.pi / 2)
    ring.segment(coupling_l / 2, "-x")

    ring2 = gdspy.copy(ring)
    ring2.mirror((0, 0), (0, 100))

    # Create the coupling region
    coupler = gdspy.Rectangle((-coupling_l / 2 - grating_offset, coupling_gap + wg_w / 2),
                              (coupling_l / 2 + grating_offset, coupling_gap + wg_w / 2 + wg_w))
    grating_taper = gdspy.CellReference(grating_taper_cell, (-coupling_l / 2 - grating_offset, wg_w + coupling_gap),
                                        180)
    coupler_grating = gdspy.CellReference(grating_cell, (
        -coupling_l / 2 - grating_offset - grating_taper_length, wg_w + coupling_gap), 180)
    reflector_taper = gdspy.CellReference(reflector_taper_cell, (coupling_l / 2 + grating_offset, wg_w + coupling_gap))
    reflector_grating = gdspy.CellReference(reflector_cell, (
        coupling_l / 2 + grating_offset + reflector_taper_length, wg_w + coupling_gap))

    # Add a grating-only region
    # grating_test_wg = gdspy.Rectangle((-coupling_l/2 - grating_offset,-2*ring_r-5), (coupling_l/2 + grating_offset,-2*ring_r-5-wg_w))
    # grating_test_in = gdspy.CellReference(grating_cell,(-coupling_l/2 - grating_offset, -2*ring_r-5-wg_w/2),180)
    # grating_test_out = gdspy.CellReference(grating_cell,(coupling_l/2 + grating_offset, -2*ring_r-5-wg_w/2),0)

    # Add identifying text

    device = gdspy.Cell(str(g) + str(wg_w) + str(coupling_gap) + 'device')

    text = gdspy.Text(str(int(wg_w * 1000)), 4, (-2.5, -30))
    device.add(text)

    text = gdspy.Text(str(int(coupling_gap * 1000)), 4, (-2.5, -20))
    device.add(text)

    device.add(ring)
    device.add(ring2)
    device.add(coupler)
    device.add(grating_taper)
    device.add(coupler_grating)
    device.add(reflector_taper)
    device.add(reflector_grating)

    return device


cell_membrane_row = gdspy.Cell('membrane_width_sweep')  # Single instance of a sweep over the membrane widths
cell_membrane_row_nohole = gdspy.Cell(
    'membrane_width_sweep_nohole')  # Single instance of a sweep over the membrane widths
cell_support = gdspy.Cell('support_structure')  # Single instance of the support structure for the membrane layer

# center_hole_size = 10 #Size of square holes in the center of the membrane


num_membranes = num_cols * 1  # Number of membranes in a given row - repeat each membrane three times
row_length = num_membranes * spacing  # Length of a row
column_size = num_rows * spacing  # Length of a column

# Define the cell reference for the support structure
support_vertices = [(-gap / 2, support_width / 2), (-support_pinch_length / 2, support_pinch_width / 2),
                    (0, support_pinch_width / 2), (0, 0), (-gap / 2, 0)]
support_1 = gdspy.Polygon(support_vertices, layer=5)
support_2 = gdspy.Polygon(support_vertices, layer=5)
support_3 = gdspy.Polygon(support_vertices, layer=5)
support_4 = gdspy.Polygon(support_vertices, layer=5)
support_2.mirror((0, -1), (0, 1))
support_3.mirror((-1, 0), (1, 0))
support_4.mirror((0, -1), (0, 1))
support_4.mirror((-1, 0), (1, 0))
cell_support.add(support_1)
cell_support.add(support_2)
cell_support.add(support_3)
cell_support.add(support_4)


# Define the square membrane hole cut-out
# membrane_hole = gdspy.Rectangle((-center_hole_size/2, center_hole_size/2), (center_hole_size/2, -center_hole_size/2))
# cell_membrane_hole.add(membrane_hole)


def defineMembrane(identifier, spacing, spacing_y, length, height, layer = 5):
    # Create a membrane cell for each desired membrane size

    cell_membrane_nohole = gdspy.Cell('membrane' + identifier)

    # Define the frame for the membrane
    membrane_outer_frame = gdspy.Rectangle((-spacing / 2, spacing_y / 2),
                                           (spacing / 2, -spacing_y / 2), layer=layer)
    membrane_inner_frame = gdspy.Rectangle((-length / 2 - gap, height / 2 + gap),
                                           (length / 2 + gap, -height / 2 - gap), layer=layer)
    membrane_frame = gdspy.boolean(membrane_outer_frame, membrane_inner_frame, "not", layer=layer)
    # Define the membrane itself
    membrane = gdspy.Rectangle((-length / 2, height / 2), (length / 2, -height / 2), layer=layer)

    # Add to the membrane cell
    cell_membrane_nohole.add(membrane_frame)

    # Perforations
    if perforations == True:
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
    support_period = support_width / support_pitch
    air_width = support_period * (1 - support_pitch)

    num_supports_x = math.floor(length / support_period + support_pitch)
    num_supports_y = math.floor(height / support_period + support_pitch)

    if (num_supports_y % 2) == 0:

        # Number of supports is even
        support_length_y = num_supports_y * support_period - support_period / 2
        support_offset_y = (length - support_length_y) / 2 + support_width / 2

        for i in range(0, num_supports_y):
            ref1 = gdspy.CellReference(cell_support, (
                -length / 2 - gap / 2, -length / 2 + support_offset_y + i * support_period))
            ref2 = gdspy.CellReference(cell_support, (
                length / 2 + gap / 2, -length / 2 + support_offset_y + i * support_period))
            cell_membrane_nohole.add([ref1, ref2])
    else:
        # Number of supports is odd
        for i in range(-math.floor(num_supports_y / 2), math.floor(num_supports_y / 2) + 1):
            ref1 = gdspy.CellReference(cell_support, (-length / 2 - gap / 2, i * support_period))
            ref2 = gdspy.CellReference(cell_support, (length / 2 + gap / 2, i * support_period))
            cell_membrane_nohole.add([ref1, ref2])

    # Membrane Sides For Rectangles
    if (num_supports_x % 2) == 0:

        # For rectangles
        support_length_x = num_supports_x * support_period - support_period / 2
        support_offset_x = (height - support_length_x) / 2 + support_width / 2

        for i in range(0, num_supports_x):
            # Rectangle
            ref3 = gdspy.CellReference(cell_support, (
                -height / 2 + support_offset_x + i * support_period, height / 2 + gap / 2), rotation=90)
            ref4 = gdspy.CellReference(cell_support, (
                -height / 2 + support_offset_x + i * support_period, -height / 2 - gap / 2),
                                       rotation=90)
            cell_membrane_nohole.add([ref3, ref4])

    else:
        # Number of supports is odd
        for i in range(-math.floor(num_supports_x / 2), math.floor(num_supports_x / 2) + 1):
            ref3 = gdspy.CellReference(cell_support, (i * support_period, height / 2 + gap / 2), rotation=90)
            ref4 = gdspy.CellReference(cell_support, (i * support_period, -height / 2 - gap / 2), rotation=90)
            cell_membrane_nohole.add([ref3, ref4])

    cell_membrane_ref_nohole = gdspy.CellReference(cell_membrane_nohole, (0, 0))
    cell_membrane_row_nohole.add(cell_membrane_ref_nohole)

    return cell_membrane_row_nohole


########################################
""" PhC Code"""


def write_beam_single(cell, pltdata, layer=1):
    cell.add(gdspy.Polygon([
        (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0),
        (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos + pltdata.dy / 2.0),
        (pltdata.xpos + pltdata.dx / 2.0, pltdata.ypos + pltdata.dy / 2.0),
        (pltdata.xpos + pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0),
        (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0)
    ], layer=layer))


def write_beam_array(cell, first_beam, width_list=None):
    global num_beam, beam_spacing

    if width_list is None:
        _dy_list = first_beam.dy * numpy.ones(num_beam, dtype=int)
    else:
        _dy_list = width_list

    _beam_write = copy.copy(first_beam)
    _initial_ypos = first_beam.ypos - _dy_list[0] / 2.0
    for i in range(num_beam):
        _beam_write.dy = _dy_list[i]
        _beam_write.ypos = _initial_ypos + i * beam_spacing + numpy.sum(_dy_list[:i]) + _beam_write.dy / 2

    # write_beam_single(cell, _beam_write)

    return _dy_list


def write_hole_single(cell, pltdata, layer=2):
    global num_circ_points, shot_pitch

    _philist = numpy.linspace(0, 2 * numpy.pi, num_circ_points)

    _circ_pts = [((round(pltdata.xpos / shot_pitch) + round(pltdata.dx / 2 * numpy.cos(phi) / shot_pitch)) * shot_pitch,
                  (round(pltdata.ypos / shot_pitch) + round(pltdata.dy / 2 * numpy.sin(phi) / shot_pitch)) * shot_pitch)
                 for phi in _philist]

    return gdspy.Polygon(_circ_pts, layer=layer)


def get_tapers(num_end_taper, end_period, amir):
    taper_period_list = np.linspace(amir, end_period, num_end_taper + 1)
    return taper_period_list[1:]


def exponential_cavity_period_list(amir, acav, num_cav, exponent):
    N = int(num_cav / 2)
    cavity_period_list = np.zeros(N)
    a = (amir - acav) / ((N - 1) ** exponent)
    b = acav

    for i in range(N):
        cavity_period_list[i] = a * (i ** exponent) + b

    cavity_period_list = np.append(np.flipud(cavity_period_list), cavity_period_list)
    return cavity_period_list


def write_hole_1D(cell, beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, amir, end_period,
                  end_taper_L=False, end_taper_R=False, num_end_taper=0, reverse_tone=False, layer=2):
    n = 1.5  # Exponent for polynomial

    _one_side_holes = (num_cav_hole - 1) / 2.0
    _idx_list = numpy.linspace(-_one_side_holes, _one_side_holes, int(num_cav_hole))
    X1 = _idx_list[int(num_cav_hole / 2)]  # First Hole
    XN = _idx_list[-1]  # Last Hole

    _acav_list = exponential_cavity_period_list(amir, acav, num_cav_hole, exponent=2)

    if custom_cavity == True:
        _acav_list = custom_cavity_list

    if end_taper_L == True:
        taper_periods = get_tapers(num_end_taper, end_period, amir=amir)

    _amir_list_L = np.append(np.flipud(taper_periods), [amir for x in range(int(num_mir_hole_L))])

    _amir_list_R = [amir for x in range(int(num_mir_hole_R))]
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
        _single_hole_polygon = write_hole_single(cell, _hole_write)

        if i == 0:
            _hole_polygon_combined = _single_hole_polygon
        else:
            _hole_polygon_combined = gdspy.fast_boolean(_hole_polygon_combined, _single_hole_polygon, 'or',
                                                        max_points=0, layer=layer)

        _hole_write.xpos = _hole_write.xpos + _aper_list[i] / 2.0

    if reverse_tone is True:
        _tmp_beam_polygon = gdspy.Polygon([
            (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
            (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
            (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
            (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
            (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0)
        ], layer=layer)
        _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _hole_polygon_combined, 'not', max_points=0,
                                               layer=layer)
        send = _tmp_beam_polygon
    else:
        send = _hole_polygon_combined

    return send


def write_hole_1D_mir_sweep(cell, beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, amir,
                            end_taper_L=False, end_taper_R=False, num_end_taper=0, reverse_tone=False, layer=2):
    _one_side_holes = 0
    _amir_list_L = [amir for x in range(int(num_mir_hole_L + num_end_taper * end_taper_L))]
    _amir_list_R = [amir for x in range(int(num_mir_hole_R + num_end_taper * end_taper_R))]
    _aper_list = copy.copy(_amir_list_L)
    _aper_list.extend(_amir_list_R)
    if len(_aper_list) > 0:
        _hole_write = copy.copy(holedata)
        _hole_write.xpos = _hole_write.xpos - numpy.sum(numpy.array(_amir_list_L))
        _taper_scale_list = []
        if num_end_taper > 0:
            _taper_scale_list = numpy.linspace(0, 1.0, num_end_taper + 2)
            _taper_scale_list_L = _taper_scale_list[1:-1]
            _taper_scale_list_R = numpy.flipud(_taper_scale_list[1:-1])

        for i in range(len(_aper_list)):
            _hole_write.xpos = _hole_write.xpos + _aper_list[i] / 2.0
            if i < num_end_taper * end_taper_L:
                _hole_write.dx = holedata.dx * _taper_scale_list_L[i]
                _hole_write.dy = holedata.dy * _taper_scale_list_L[i]
            else:
                _hole_write.dx = holedata.dx
                _hole_write.dy = holedata.dy
            if i >= (len(_aper_list) - num_end_taper * end_taper_R):
                _hole_write.dx = holedata.dx * _taper_scale_list_R[i - len(_aper_list) + num_end_taper]
                _hole_write.dy = holedata.dy * _taper_scale_list_R[i - len(_aper_list) + num_end_taper]
            _single_hole_polygon = write_hole_single(cell, _hole_write)

            if i == 0:
                _hole_polygon_combined = _single_hole_polygon
            else:
                _hole_polygon_combined = gdspy.fast_boolean(_hole_polygon_combined, _single_hole_polygon, 'or',
                                                            max_points=0, layer=layer)

            _hole_write.xpos = _hole_write.xpos + _aper_list[i] / 2.0

        if reverse_tone is True:
            _tmp_beam_polygon = gdspy.Polygon([
                (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
                (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
                (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
                (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
                (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0)], layer=layer)
            _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _hole_polygon_combined, 'not', max_points=0,
                                                   layer=layer)
            cell.add(_tmp_beam_polygon)
        else:
            cell.add(_hole_polygon_combined)
    else:
        _tmp_beam_polygon = gdspy.Polygon([
            (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
            (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
            (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos + beamdata.dy / 2.0),
            (beamdata.xpos + beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0),
            (beamdata.xpos - beamdata.dx / 2.0, beamdata.ypos - beamdata.dy / 2.0)], layer=layer)
        cell.add(_tmp_beam_polygon)


def write_hole_2D(cell, beamdata, holedata, beam_dy_list, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav, amir,
                  end_period, guides,
                  end_taper_L=False, end_taper_R=False, num_end_taper=0, reverse_tone=False, edge_overdose=False,
                  ring_overdose=False,
                  ring_underdose=True, edge_underdose=True):
    global num_beam, beam_spacing, over_edge_size, over_ring_size, spacer

    _initial_ypos = holedata.ypos - beam_dy_list[0] / 2.0
    _tmp_beamdata = copy.copy(beamdata)
    _tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0] / 2.0

    for i in range(num_beam):
        holedata.ypos = _initial_ypos + i * beam_spacing + numpy.sum(beam_dy_list[:i]) + beam_dy_list[i] / 2.0
        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0
        _tmp_beamdata.dy = beam_dy_list[i]


        if guides == 1 and i == 0:
            write_beam_single(cell, _tmp_beamdata, layer=2)
        elif guides == 2 and (i == 0 or i == (num_beam-1)):
            write_beam_single(cell, _tmp_beamdata, layer=2)
        else:

            ord_holes = write_hole_1D(cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L, num_mir_hole_R, acav,
                                      amir,
                                      end_period,
                                      end_taper_L, end_taper_R, num_end_taper, reverse_tone=reverse_tone)

            if ring_overdose or edge_overdose is True:
                if reverse_tone == True:
                    holedata.dx = holedata.dx + over_ring_size * ring_overdose * 2
                    holedata.dy = holedata.dy + over_ring_size * ring_overdose * 2

                    _tmp_beamdata.dy = _tmp_beamdata.dy - over_edge_size * 2 * edge_overdose

                    overdose_holes = write_hole_1D(cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L,
                                                   num_mir_hole_R, acav, amir,
                                                   end_period,
                                                   end_taper_L, end_taper_R, num_end_taper, reverse_tone=reverse_tone)

                    region = gdspy.fast_boolean(ord_holes, overdose_holes, 'not', max_points=0, layer=11)
                    cell.add(region)
                    cell.add(overdose_holes)

                    holedata.dx = holedata.dx - over_ring_size * ring_overdose * 2
                    holedata.dy = holedata.dy - over_ring_size * ring_overdose * 2

            elif ring_underdose or edge_underdose is True:
                if reverse_tone is False:
                    ref_ord_edge = write_hole_1D(cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L,
                                                 num_mir_hole_R, acav, amir,
                                                 end_period,
                                                 end_taper_L, end_taper_R, num_end_taper, reverse_tone=True)

                    holedata.dx = holedata.dx - under_ring_size * ring_underdose
                    holedata.dy = holedata.dy - under_ring_size * ring_underdose
                    _tmp_beamdata.dy = _tmp_beamdata.dy + under_edge_size * 2 * edge_underdose

                    underdose_holes = write_hole_1D(cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L,
                                                    num_mir_hole_R, acav, amir,
                                                    end_period,
                                                    end_taper_L, end_taper_R, num_end_taper, reverse_tone=reverse_tone)

                    _tmp_beamdata.dx = _tmp_beamdata.dx - spacer * 2
                    underdose_edge = write_hole_1D(cell, _tmp_beamdata, holedata, num_cav_hole, num_mir_hole_L,
                                                   num_mir_hole_R, acav, amir,
                                                   end_period,
                                                   end_taper_L, end_taper_R, num_end_taper, reverse_tone=True)

                    region = gdspy.fast_boolean(ord_holes, underdose_holes, 'not', max_points=0, layer=12)
                    region2 = gdspy.fast_boolean(underdose_edge, ref_ord_edge, 'not', max_points=0, layer=12)

                    # underdose_holes = gdspy.fast_boolean(underdose_edge, region2, 'not', max_points=0, layer=12)

                    cell.add(underdose_holes)
                    cell.add(region)
                    cell.add(region2)
            else:
                cell.add(ord_holes)

        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0 + beam_spacing


def write_hole_2D_mir_sweep(cell, beamdata, holedata, beam_dy_list, num_cav_hole, num_mir_hole_L, num_mir_hole_R,
                            acav, amir, guides,
                            end_taper_L=False, end_taper_R=False, num_end_taper=0, reverse_tone=False,
                            edge_overdose=False, ring_overdose=False):
    global num_beam, beam_spacing

    _initial_ypos = holedata.ypos - beam_dy_list[0] / 2.0
    _tmp_beamdata = copy.copy(beamdata)
    _tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0] / 2.0

    for i in range(num_beam):
        num_mir = i
        holedata.ypos = _initial_ypos + i * beam_spacing + numpy.sum(beam_dy_list[:i]) + beam_dy_list[i] / 2.0
        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0
        _tmp_beamdata.dy = beam_dy_list[i]
        write_hole_1D_mir_sweep(cell, _tmp_beamdata, holedata, num_cav_hole, num_mir, num_mir, acav, amir,
                                end_taper_L, end_taper_R, num_end_taper, reverse_tone=reverse_tone)
        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0 + beam_spacing


def write_alignment_mark(cell, alignment_xpos, alignment_ypos, layer=3):
    global mark_thickness, mark_length

    cell.add(gdspy.Polygon([
        (alignment_xpos - mark_length / 2.0, alignment_ypos - mark_thickness / 2.0),
        (alignment_xpos - mark_length / 2.0, alignment_ypos + mark_thickness / 2.0),
        (alignment_xpos + mark_length / 2.0, alignment_ypos + mark_thickness / 2.0),
        (alignment_xpos + mark_length / 2.0, alignment_ypos - mark_thickness / 2.0),
        (alignment_xpos - mark_length / 2.0, alignment_ypos - mark_thickness / 2.0)
    ], layer=layer))
    cell.add(gdspy.Polygon([
        (alignment_xpos - mark_thickness / 2.0, alignment_ypos - mark_length / 2.0),
        (alignment_xpos - mark_thickness / 2.0, alignment_ypos + mark_length / 2.0),
        (alignment_xpos + mark_thickness / 2.0, alignment_ypos + mark_length / 2.0),
        (alignment_xpos + mark_thickness / 2.0, alignment_ypos - mark_length / 2.0),
        (alignment_xpos - mark_thickness / 2.0, alignment_ypos - mark_length / 2.0)
    ], layer=layer))


def pinch_pt_polygon(beam_dy, pinch_pt_size, pinch_pt_xloc, pinch_pt_yloc, layer=5):
    _upper_triangle = gdspy.Polygon([
        (pinch_pt_xloc, pinch_pt_yloc + pinch_pt_size / 2.0),
        (pinch_pt_xloc - (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0, pinch_pt_yloc + beam_dy / 2.0),
        (pinch_pt_xloc + (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0, pinch_pt_yloc + beam_dy / 2.0),
        (pinch_pt_xloc, pinch_pt_yloc + pinch_pt_size / 2.0)
    ], layer=layer)
    _lower_triangle = gdspy.Polygon([
        (pinch_pt_xloc, pinch_pt_yloc - pinch_pt_size / 2.0),
        (pinch_pt_xloc - (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0, pinch_pt_yloc - beam_dy / 2.0),
        (pinch_pt_xloc + (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0, pinch_pt_yloc - beam_dy / 2.0),
        (pinch_pt_xloc, pinch_pt_yloc - pinch_pt_size / 2.0)
    ], layer=layer)

    return _upper_triangle, _lower_triangle


def pinch_pt_polygon_vertical(beam_dy, pinch_pt_size, pinch_pt_xloc, pinch_pt_yloc, layer=3):
    _left_triangle = gdspy.Polygon([
        (pinch_pt_xloc - pinch_pt_size / 2.0, pinch_pt_yloc),
        (pinch_pt_xloc - beam_dy / 2.0, pinch_pt_yloc - (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0),
        (pinch_pt_xloc - beam_dy / 2.0, pinch_pt_yloc + (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0),
        (pinch_pt_xloc - pinch_pt_size / 2.0, pinch_pt_yloc)
    ], layer=layer)
    _right_triangle = gdspy.Polygon([
        (pinch_pt_xloc + pinch_pt_size / 2.0, pinch_pt_yloc),
        (pinch_pt_xloc + beam_dy / 2.0, pinch_pt_yloc - (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0),
        (pinch_pt_xloc + beam_dy / 2.0, pinch_pt_yloc + (beam_dy - pinch_pt_size) * math.sqrt(3) / 2.0),
        (pinch_pt_xloc + pinch_pt_size / 2.0, pinch_pt_yloc)
    ], layer=layer)

    return _left_triangle, _right_triangle


def linker_polygon(pltdata, layer=3):
    _linker = gdspy.Polygon([
        (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0),
        (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos + pltdata.dy / 2.0),
        (pltdata.xpos + pltdata.dx / 2.0, pltdata.ypos + pltdata.dy / 2.0),
        (pltdata.xpos + pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0),
        (pltdata.xpos - pltdata.dx / 2.0, pltdata.ypos - pltdata.dy / 2.0)
    ], layer=layer)

    return _linker


def write_linker_region(beamdata, xmin, xmax, ymin, ymax, round_corner=False, layer=6):
    global corner_bend_rad, corner_bend_pts, linker_edgeoffset, linker_width, linker_connector_width
    global linker_notch_size, linker_beam_size, linker_xnotches, linker_ynotches

    _tmp_beamdata = copy.copy(beamdata)
    _linkerdata_inner = PLTdata()
    _linkerdata_inner.xpos = _tmp_beamdata.xpos
    _linkerdata_inner.ypos = (ymin + ymax) / 2.0
    _linkerdata_inner.dx = _tmp_beamdata.dx
    _linkerdata_inner.dy = ymax - ymin - 2.0 * linker_edgeoffset - 2.0 * linker_connector_width
    _tmp_linker_inner = linker_polygon(_linkerdata_inner, layer)

    _linkerdata_outer = PLTdata()
    _linkerdata_outer.xpos = _tmp_beamdata.xpos
    _linkerdata_outer.ypos = (ymin + ymax) / 2.0
    _linkerdata_outer.dx = _tmp_beamdata.dx + 2.0 * linker_width
    _linkerdata_outer.dy = ymax - ymin - 2.0 * linker_edgeoffset
    _tmp_linker_outer = linker_polygon(_linkerdata_outer, layer)

    if round_corner is True:
        # _tmp_linker_inner = _tmp_linker_inner.fillet(corner_bend_rad, points_per_2pi=corner_bend_pts)
        _tmp_linker_outer = _tmp_linker_outer.fillet(corner_bend_rad, points_per_2pi=corner_bend_pts)

    _tmp_notch_yarray = (_linkerdata_outer.dy - corner_bend_rad * 2 - linker_beam_size) * numpy.linspace(0, 1.0,
                                                                                                         num=linker_ynotches)
    for i in range(linker_ynotches):
        _tmp_notch_ypos = ymin + linker_edgeoffset + corner_bend_rad + linker_beam_size / 2.0 + _tmp_notch_yarray[i]
        _tmp_beam_polygon_L = gdspy.Polygon([
            (xmin, _tmp_notch_ypos - linker_beam_size / 2.0),
            (xmin, _tmp_notch_ypos + linker_beam_size / 2.0),
            (xmin + linker_edgeoffset, _tmp_notch_ypos + linker_beam_size / 2.0),
            (xmin + linker_edgeoffset, _tmp_notch_ypos - linker_beam_size / 2.0),
            (xmin, _tmp_notch_ypos - linker_beam_size / 2.0)
        ], layer=layer)
        _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(linker_beam_size, linker_notch_size,
                                                                    xmin + linker_edgeoffset / 2.0, _tmp_notch_ypos,
                                                                    layer)
        _tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L, _tmp_upper_triangle, 'not', max_points=0,
                                                 layer=layer)
        _tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L, _tmp_lower_triangle, 'not', max_points=0,
                                                 layer=layer)

        _tmp_beam_polygon_R = gdspy.Polygon([
            (xmax - linker_edgeoffset, _tmp_notch_ypos - linker_beam_size / 2.0),
            (xmax - linker_edgeoffset, _tmp_notch_ypos + linker_beam_size / 2.0),
            (xmax, _tmp_notch_ypos + linker_beam_size / 2.0),
            (xmax, _tmp_notch_ypos - linker_beam_size / 2.0),
            (xmax - linker_edgeoffset, _tmp_notch_ypos - linker_beam_size / 2.0)
        ], layer=layer)
        _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(linker_beam_size, linker_notch_size,
                                                                    xmax - linker_edgeoffset / 2.0, _tmp_notch_ypos,
                                                                    layer)
        _tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R, _tmp_upper_triangle, 'not', max_points=0,
                                                 layer=layer)
        _tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R, _tmp_lower_triangle, 'not', max_points=0,
                                                 layer=layer)

        _tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer, _tmp_beam_polygon_L, 'or', max_points=0, layer=layer)
        _tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer, _tmp_beam_polygon_R, 'or', max_points=0, layer=layer)

    _tmp_notch_xarray = (_linkerdata_outer.dx - corner_bend_rad * 2 - linker_beam_size) * numpy.linspace(0, 1.0,
                                                                                                         num=linker_xnotches)
    for i in range(linker_xnotches):
        _tmp_notch_xpos = xmin + linker_edgeoffset + corner_bend_rad + linker_beam_size / 2.0 + _tmp_notch_xarray[i]
        _tmp_beam_polygon_T = gdspy.Polygon([
            (_tmp_notch_xpos - linker_beam_size / 2.0, ymax),
            (_tmp_notch_xpos - linker_beam_size / 2.0, ymax - linker_edgeoffset),
            (_tmp_notch_xpos + linker_beam_size / 2.0, ymax - linker_edgeoffset),
            (_tmp_notch_xpos + linker_beam_size / 2.0, ymax),
            (_tmp_notch_xpos - linker_beam_size / 2.0, ymax)
        ], layer=layer)
        _tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(linker_beam_size, linker_notch_size,
                                                                            _tmp_notch_xpos,
                                                                            ymax - linker_edgeoffset / 2.0, layer)
        _tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T, _tmp_left_triangle, 'not', max_points=0,
                                                 layer=layer)
        _tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T, _tmp_right_triangle, 'not', max_points=0,
                                                 layer=layer)

        _tmp_beam_polygon_B = gdspy.Polygon([
            (_tmp_notch_xpos - linker_beam_size / 2.0, ymin),
            (_tmp_notch_xpos - linker_beam_size / 2.0, ymin + linker_edgeoffset),
            (_tmp_notch_xpos + linker_beam_size / 2.0, ymin + linker_edgeoffset),
            (_tmp_notch_xpos + linker_beam_size / 2.0, ymin),
            (_tmp_notch_xpos - linker_beam_size / 2.0, ymin)
        ], layer=layer)
        _tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(linker_beam_size, linker_notch_size,
                                                                            _tmp_notch_xpos,
                                                                            ymin + linker_edgeoffset / 2.0, layer)
        _tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B, _tmp_left_triangle, 'not', max_points=0,
                                                 layer=layer)
        _tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B, _tmp_right_triangle, 'not', max_points=0,
                                                 layer=layer)

        _tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer, _tmp_beam_polygon_T, 'or', max_points=0, layer=layer)
        _tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer, _tmp_beam_polygon_B, 'or', max_points=0, layer=layer)

    _linker_combined = gdspy.fast_boolean(_tmp_linker_outer, _tmp_linker_inner, 'not', max_points=0, layer=layer)

    return _linker_combined


def write_left_circ_grating(beamdata, layer=4):
    global shot_pitch, num_circ_grating_points, grating_spacing, grating_linewidth, num_grating
    global circ_grating_base, grating_angle, odd_support_angle, even_support_angle, support_angle_width

    _philist_L = numpy.linspace(numpy.pi / 2.0, numpy.pi * 3.0 / 2.0, num_circ_grating_points - 1)
    _philist_L = numpy.append(_philist_L, numpy.pi / 2.0)
    _tmp_beamdata = copy.copy(beamdata)
    _ini_pt = [(_tmp_beamdata.xpos - _tmp_beamdata.dx / 2, _tmp_beamdata.ypos)]
    #
    # rough = gdspy.Rectangle((_tmp_beamdata.xpos - _tmp_beamdata.dx / 2, _tmp_beamdata.ypos - 4000),
    #                         (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2 + 500, _tmp_beamdata.ypos - 4000), layer=9)

    _radius_inner = circ_grating_base

    for i in range(num_grating + 1):
        _left_grating_inner = gdspy.Polygon([
            ((round((_tmp_beamdata.xpos - _tmp_beamdata.dx / 2) / shot_pitch) + round(
                _radius_inner * numpy.cos(phi) / shot_pitch)) * shot_pitch,
             (round(_tmp_beamdata.ypos / shot_pitch) + round(_radius_inner * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
            for phi in _philist_L], layer=layer)

        _radius_outer = _radius_inner + grating_spacing
        _left_grating_outer = gdspy.Polygon([
            ((round((_tmp_beamdata.xpos - _tmp_beamdata.dx / 2) / shot_pitch) + round(
                _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
             (round(_tmp_beamdata.ypos / shot_pitch) + round(_radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
            for phi in _philist_L], layer=layer)

        _left_grating_tmp = gdspy.fast_boolean(_left_grating_outer, _left_grating_inner, 'not', max_points=0,
                                               layer=layer)

        if (i % 2 == 0):
            _radius_outer = _radius_outer + 10
            _philist_support = numpy.linspace(numpy.pi / 2.0 + odd_support_angle - support_angle_width / 2.0,
                                              numpy.pi / 2.0 + odd_support_angle + support_angle_width / 2.0,
                                              num_circ_grating_points - 1)
            _support_pt = [
                ((round((_tmp_beamdata.xpos - _tmp_beamdata.dx / 2) / shot_pitch) + round(
                    _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
                 (round(_tmp_beamdata.ypos / shot_pitch) + round(
                     _radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
                for phi in _philist_support]
            _support_pt_combined = copy.copy(_ini_pt)
            _support_pt_combined.extend(_support_pt)
            _support_pt_combined.extend(_ini_pt)
            _support_frame = gdspy.Polygon(_support_pt_combined, layer=layer)
            _left_grating_tmp = gdspy.fast_boolean(_left_grating_tmp, _support_frame, 'not', max_points=0, layer=layer)

            _philist_support = numpy.linspace(numpy.pi * 3.0 / 2.0 - odd_support_angle - support_angle_width / 2.0,
                                              numpy.pi * 3.0 / 2.0 - odd_support_angle + support_angle_width / 2.0,
                                              num_circ_grating_points - 1)
            _support_pt = [
                ((round((_tmp_beamdata.xpos - _tmp_beamdata.dx / 2) / shot_pitch) + round(
                    _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
                 (round(_tmp_beamdata.ypos / shot_pitch) + round(
                     _radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
                for phi in _philist_support]
            _support_pt_combined = copy.copy(_ini_pt)
            _support_pt_combined.extend(_support_pt)
            _support_pt_combined.extend(_ini_pt)
            _support_frame = gdspy.Polygon(_support_pt_combined, layer=layer)
            _left_grating_tmp = gdspy.fast_boolean(_left_grating_tmp, _support_frame, 'not', max_points=0, layer=layer)
        else:
            _radius_outer = _radius_outer + 10
            _philist_support = numpy.linspace(numpy.pi / 2.0 + even_support_angle - support_angle_width / 2.0,
                                              numpy.pi / 2.0 + even_support_angle + support_angle_width / 2.0,
                                              num_circ_grating_points - 1)
            _support_pt = [
                ((round((_tmp_beamdata.xpos - _tmp_beamdata.dx / 2) / shot_pitch) + round(
                    _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
                 (round(_tmp_beamdata.ypos / shot_pitch) + round(
                     _radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
                for phi in _philist_support]
            _support_pt_combined = copy.copy(_ini_pt)
            _support_pt_combined.extend(_support_pt)
            _support_pt_combined.extend(_ini_pt)
            _support_frame = gdspy.Polygon(_support_pt_combined, layer=layer)
            _left_grating_tmp = gdspy.fast_boolean(_left_grating_tmp, _support_frame, 'not', max_points=0, layer=layer)

        if i == 0:
            _left_grating = _left_grating_tmp
        else:
            _left_grating = gdspy.fast_boolean(_left_grating, _left_grating_tmp, 'or', max_points=0, layer=layer)

        _radius_inner = _radius_outer + grating_linewidth

    _philist_frame = numpy.linspace(numpy.pi / 2.0 + grating_angle, numpy.pi * 3.0 / 2.0 - grating_angle,
                                    num_circ_grating_points - 1)
    _grating_frame_pt = [
        ((round((_tmp_beamdata.xpos - _tmp_beamdata.dx / 2) / shot_pitch) + round(
            _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
         (round(_tmp_beamdata.ypos / shot_pitch) + round(_radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
        for phi in _philist_frame]

    _grating_frame_pt_combined = copy.copy(_ini_pt)
    _grating_frame_pt_combined.extend(_grating_frame_pt)
    _grating_frame_pt_combined.extend(_ini_pt)

    _grating_frame = gdspy.Polygon(_grating_frame_pt_combined, layer=layer)

    _left_grating = gdspy.fast_boolean(_left_grating, _grating_frame, 'and', max_points=0, layer=layer)

    return _left_grating


def write_right_circ_grating(beamdata, layer=5):
    global shot_pitch, num_circ_grating_points, grating_spacing, grating_linewidth, num_grating
    global circ_grating_base, grating_angle, odd_support_angle, even_support_angle, support_angle_width

    _philist_R = numpy.linspace(-1 * numpy.pi / 2.0, numpy.pi / 2.0, num_circ_grating_points - 1)
    _philist_R = numpy.append(_philist_R, -1 * numpy.pi / 2.0)
    _tmp_beamdata = copy.copy(beamdata)
    _ini_pt = [(_tmp_beamdata.xpos + _tmp_beamdata.dx / 2, _tmp_beamdata.ypos)]

    _radius_inner = circ_grating_base

    for i in range(num_grating + 1):
        _right_grating_inner = gdspy.Polygon([
            ((round((_tmp_beamdata.xpos + _tmp_beamdata.dx / 2) / shot_pitch) + round(
                _radius_inner * numpy.cos(phi) / shot_pitch)) * shot_pitch,
             (round(_tmp_beamdata.ypos / shot_pitch) + round(_radius_inner * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
            for phi in _philist_R], layer=layer)

        _radius_outer = _radius_inner + grating_spacing
        _right_grating_outer = gdspy.Polygon([
            ((round((_tmp_beamdata.xpos + _tmp_beamdata.dx / 2) / shot_pitch) + round(
                _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
             (round(_tmp_beamdata.ypos / shot_pitch) + round(_radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
            for phi in _philist_R], layer=layer)

        _right_grating_tmp = gdspy.fast_boolean(_right_grating_outer, _right_grating_inner, 'not', max_points=0,
                                                layer=layer)

        if (i % 2 == 0):
            _radius_outer = _radius_outer + 10
            _philist_support = numpy.linspace(-1 * numpy.pi / 2.0 + odd_support_angle - support_angle_width / 2.0,
                                              -1 * numpy.pi / 2.0 + odd_support_angle + support_angle_width / 2.0,
                                              num_circ_grating_points - 1)
            _support_pt = [
                ((round((_tmp_beamdata.xpos + _tmp_beamdata.dx / 2) / shot_pitch) + round(
                    _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
                 (round(_tmp_beamdata.ypos / shot_pitch) + round(
                     _radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
                for phi in _philist_support]
            _support_pt_combined = copy.copy(_ini_pt)
            _support_pt_combined.extend(_support_pt)
            _support_pt_combined.extend(_ini_pt)
            _support_frame = gdspy.Polygon(_support_pt_combined, layer=layer)
            _right_grating_tmp = gdspy.fast_boolean(_right_grating_tmp, _support_frame, 'not', max_points=0,
                                                    layer=layer)

            _philist_support = numpy.linspace(numpy.pi / 2.0 - odd_support_angle - support_angle_width / 2.0,
                                              numpy.pi / 2.0 - odd_support_angle + support_angle_width / 2.0,
                                              num_circ_grating_points - 1)
            _support_pt = [
                ((round((_tmp_beamdata.xpos + _tmp_beamdata.dx / 2) / shot_pitch) + round(
                    _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
                 (round(_tmp_beamdata.ypos / shot_pitch) + round(
                     _radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
                for phi in _philist_support]
            _support_pt_combined = copy.copy(_ini_pt)
            _support_pt_combined.extend(_support_pt)
            _support_pt_combined.extend(_ini_pt)
            _support_frame = gdspy.Polygon(_support_pt_combined, layer=layer)
            _right_grating_tmp = gdspy.fast_boolean(_right_grating_tmp, _support_frame, 'not', max_points=0,
                                                    layer=layer)
        else:
            _radius_outer = _radius_outer + 10
            _philist_support = numpy.linspace(-1 * numpy.pi / 2.0 + even_support_angle - support_angle_width / 2.0,
                                              -1 * numpy.pi / 2.0 + even_support_angle + support_angle_width / 2.0,
                                              num_circ_grating_points - 1)
            _support_pt = [
                ((round((_tmp_beamdata.xpos + _tmp_beamdata.dx / 2) / shot_pitch) + round(
                    _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
                 (round(_tmp_beamdata.ypos / shot_pitch) + round(
                     _radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
                for phi in _philist_support]
            _support_pt_combined = copy.copy(_ini_pt)
            _support_pt_combined.extend(_support_pt)
            _support_pt_combined.extend(_ini_pt)
            _support_frame = gdspy.Polygon(_support_pt_combined, layer=layer)
            _right_grating_tmp = gdspy.fast_boolean(_right_grating_tmp, _support_frame, 'not', max_points=0,
                                                    layer=layer)

        if i == 0:
            _right_grating = _right_grating_tmp
        else:
            _right_grating = gdspy.fast_boolean(_right_grating, _right_grating_tmp, 'or', max_points=0, layer=layer)

        _radius_inner = _radius_outer + grating_linewidth

    _philist_frame = numpy.linspace(-1 * numpy.pi / 2.0 + grating_angle, numpy.pi / 2.0 - grating_angle,
                                    num_circ_grating_points - 1)
    _grating_frame_pt = [
        ((round((_tmp_beamdata.xpos + _tmp_beamdata.dx / 2) / shot_pitch) + round(
            _radius_outer * numpy.cos(phi) / shot_pitch)) * shot_pitch,
         (round(_tmp_beamdata.ypos / shot_pitch) + round(_radius_outer * numpy.sin(phi) / shot_pitch)) * shot_pitch) \
        for phi in _philist_frame]

    _grating_frame_pt_combined = copy.copy(_ini_pt)
    _grating_frame_pt_combined.extend(_grating_frame_pt)
    _grating_frame_pt_combined.extend(_ini_pt)

    _grating_frame = gdspy.Polygon(_grating_frame_pt_combined, layer=layer)

    _right_grating = gdspy.fast_boolean(_right_grating, _grating_frame, 'and', max_points=0, layer=layer)

    return _right_grating


def write_circ_grating(beamdata, beam_dy_list, circ_grating_support=None, layer=5):
    global num_beam, beam_spacing

    _tmp_beamdata = copy.copy(beamdata)
    _tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0] / 2.0
    for i in range(num_beam):
        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0
        _tmp_beamdata.dy = beam_dy_list[i]

        _circ_grating_L = write_left_circ_grating(_tmp_beamdata, layer=layer)
        _circ_grating_R = write_right_circ_grating(_tmp_beamdata, layer=layer)

        if i == 0:
            _circ_grating_combined = _circ_grating_L
            _circ_grating_combined = gdspy.fast_boolean(_circ_grating_combined, _circ_grating_R, 'or', max_points=0,
                                                        layer=layer)
        else:
            _circ_grating_combined = gdspy.fast_boolean(_circ_grating_combined, _circ_grating_L, 'or', max_points=0,
                                                        layer=layer)
            _circ_grating_combined = gdspy.fast_boolean(_circ_grating_combined, _circ_grating_R, 'or', max_points=0,
                                                        layer=layer)

        _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0 + beam_spacing

    return _circ_grating_combined


from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath


def write_pattern_number(pattern_number, size, position, font_prop=None, tolerance=0.1):
    text = "+" + str(pattern_number)
    path = TextPath(position, text, size=float(size), prop=font_prop)
    polys = []
    xmax = position[0]
    for points, code in path.iter_segments():
        if code == path.MOVETO:
            c = gdspy.Curve(*points, tolerance=tolerance)
        elif code == path.LINETO:
            c.L(*points)
        elif code == path.CURVE3:
            c.Q(*points)
        elif code == path.CURVE4:
            c.C(*points)
        elif code == path.CLOSEPOLY:
            poly = c.get_points()
            if poly.size > 0:
                if poly[:, 0].min() < xmax:
                    i = len(polys) - 1
                    while i >= 0:
                        if gdspy.inside(
                            poly[:1], [polys[i]], precision=0.1 * tolerance
                        )[0]:
                            p = polys.pop(i)
                            poly = gdspy.boolean(
                                [p],
                                [poly],
                                "xor",
                                precision=0.1 * tolerance,
                                max_points=0,
                            ).polygons[0]
                            break
                        elif gdspy.inside(
                            polys[i][:1], [poly], precision=0.1 * tolerance
                        )[0]:
                            p = polys.pop(i)
                            poly = gdspy.boolean(
                                [p],
                                [poly],
                                "xor",
                                precision=0.1 * tolerance,
                                max_points=0,
                            ).polygons[0]
                        i -= 1
                xmax = max(xmax, poly[:, 0].max())
                polys.append(poly)
    return polys


def write_outer_box(cell, beamdata, beam_dy_list, grating_spacer=False, round_corner=False, direct_write_area=False,
                    alignment_mark=False,
                    write_linker=False, pinch_pt_L=False, pinch_pt_R=False, pinch_pt_L_offset=0, pinch_pt_R_offset=0,
                    pinch_pt_size=0, circ_grating=False, grating_support_size=None, pattern_number=None,
                    reverse_tone=False, layer=3):
    global num_beam, beam_spacing, edge_offset, corner_bend_rad, corner_bend_pts, linker_edgeoffset, spacer
    global linker_width, text_dist_to_top, under_edge_size

    _xmin = beamdata.xpos - beamdata.dx / 2.0 - int(write_linker) * (linker_width + linker_edgeoffset)
    _xmax = beamdata.xpos + beamdata.dx / 2.0 + int(write_linker) * (linker_width + linker_edgeoffset)
    _ymin = beamdata.ypos - beam_dy_list[0] / 2.0 - edge_offset
    _ymax = _ymin + (num_beam - 1) * beam_spacing + numpy.sum(beam_dy_list) + edge_offset * 2

    _outer_box = gdspy.Polygon([
        (_xmin, _ymin), (_xmin, _ymax), (_xmax, _ymax), (_xmax, _ymin), (_xmin, _ymin)
    ], layer=layer)

    if round_corner is True:
        _outer_box = _outer_box.fillet(corner_bend_rad, points_per_2pi=corner_bend_pts)

    if direct_write_area is True:
        _tmp_beamdata = copy.copy(beamdata)
        _tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0] / 2.0
        for i in range(num_beam):
            _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0
            if edge_underdose is True and reverse_tone is False:
                _tmp_beamdata.dy = beam_dy_list[i] + under_edge_size * 2
            else:
                _tmp_beamdata.dy = beam_dy_list[i]
            _tmp_beam_polygon = gdspy.Polygon([
                (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
                (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy / 2.0),
                (_tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy / 2.0),
                (_tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0),
                (_tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy / 2.0)
            ], layer=layer)

            if pinch_pt_L is True:
                _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(_tmp_beamdata.dy, pinch_pt_size,
                                                                            _tmp_beamdata.xpos - _tmp_beamdata.dx / 2.0 + pinch_pt_L_offset,
                                                                            _tmp_beamdata.ypos, layer)
                _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _tmp_upper_triangle, 'not', max_points=0,
                                                       layer=layer)
                _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _tmp_lower_triangle, 'not', max_points=0,
                                                       layer=layer)
            if pinch_pt_R is True:
                _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(_tmp_beamdata.dy, pinch_pt_size,
                                                                            _tmp_beamdata.xpos + _tmp_beamdata.dx / 2.0 - pinch_pt_R_offset,
                                                                            _tmp_beamdata.ypos, layer)
                _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _tmp_upper_triangle, 'not', max_points=0,
                                                       layer=layer)
                _tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon, _tmp_lower_triangle, 'not', max_points=0,
                                                       layer=layer)

            _tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i] / 2.0 + beam_spacing
            if i == 0:
                _tmp_beam_polygon_combined = _tmp_beam_polygon
            else:
                _tmp_beam_polygon_combined = gdspy.fast_boolean(_tmp_beam_polygon_combined, _tmp_beam_polygon, 'or',
                                                                max_points=0, layer=layer)

        _box_write_area = gdspy.fast_boolean(_outer_box, _tmp_beam_polygon_combined, 'not', max_points=0, layer=layer)

        if write_linker is True:
            _linker_combined = write_linker_region(beamdata, _xmin, _xmax, _ymin, _ymax, round_corner, layer=layer)
            _box_write_area = gdspy.fast_boolean(_box_write_area, _linker_combined, 'not', max_points=0, layer=layer)

        if circ_grating is True:
            _circ_grating_combined = write_circ_grating(beamdata, beam_dy_list, grating_support_size, layer=5)
            _box_write_area = gdspy.fast_boolean(_box_write_area, _circ_grating_combined, 'or', max_points=0,
                                                 layer=layer)

        if pattern_number is not None:
            _left_pattern_number = write_pattern_number(cell, _xmin + linker_edgeoffset + linker_width / 2.0,
                                                        _ymax - linker_edgeoffset - text_dist_to_top, pattern_number)
            _right_pattern_number = write_pattern_number(cell, _xmax - linker_edgeoffset - linker_width / 2.0,
                                                         _ymax - linker_edgeoffset - text_dist_to_top, pattern_number)
            _pattern_number_combined = gdspy.fast_boolean(_left_pattern_number, _right_pattern_number, 'or',
                                                          max_points=0, layer=layer)
            _box_write_area = gdspy.fast_boolean(_box_write_area, _pattern_number_combined, 'xor', max_points=0,
                                                 layer=layer)

    rough = gdspy.Rectangle((beamdata.xpos - beamdata.dx / 2, beamdata.ypos - edge_offset - beamdata.dy / 2),
                            (beamdata.xpos - beamdata.dx / 2 + spacer, beamdata.ypos + edge_offset + beamdata.dy / 2),
                            layer=layer)

    rough2 = gdspy.Rectangle((beamdata.xpos + beamdata.dx / 2, beamdata.ypos - edge_offset - beamdata.dy / 2),
                             (beamdata.xpos + beamdata.dx / 2 - spacer, beamdata.ypos + edge_offset + beamdata.dy / 2),
                             layer=layer)

    rough2 = gdspy.fast_boolean(rough2, _tmp_beam_polygon_combined, 'not', max_points=0, layer=layer)
    rough = gdspy.fast_boolean(rough, _tmp_beam_polygon_combined, 'not', max_points=0, layer=layer)

    if reverse_tone is True:
        _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_combined, 'or', max_points=0,
                                             layer=layer)

        _box_write_area_reverse = gdspy.fast_boolean(_outer_box, _box_write_area, 'not', max_points=0, layer=layer)

        if grating_spacer == True:
            cell.add(rough)
            cell.add(rough2)

        cell.add(_box_write_area_reverse)

    else:
        if grating_spacer == True:
            _box_write_area = gdspy.fast_boolean(_box_write_area, rough, 'not', max_points=0, layer=layer)
            _box_write_area = gdspy.fast_boolean(_box_write_area, rough2, 'not', max_points=0, layer=layer)
        cell.add(_box_write_area)

    if alignment_mark:
        _yoffset = 10000
        _alignment_x = beamdata.xpos
        _alignment_y = _ymax + _yoffset
        write_alignment_mark(cell, _alignment_x, _alignment_y)

    return _outer_box


def write_support_region(innerframe, layer=3):
    global corner_bend_rad, support_connector_width, support_notch_size, support_beam_size
    global num_xsupport, num_ysupport, write_field_x_size, write_field_y_size

    _support_box_outer = copy.copy(innerframe)
    _ymin = _support_box_outer.ypos - _support_box_outer.dy / 2.0
    _ymax = _support_box_outer.ypos + _support_box_outer.dy / 2.0
    _xmin = _support_box_outer.xpos - _support_box_outer.dx / 2.0
    _xmax = _support_box_outer.xpos + _support_box_outer.dx / 2.0
    _tmp_support_outer = linker_polygon(_support_box_outer, layer)

    _support_box_inner = copy.copy(_support_box_outer)
    _support_box_inner.dx = _support_box_inner.dx - 2 * overlap_width
    _support_box_inner.dy = _support_box_inner.dy - 2 * overlap_width
    _tmp_support_inner = linker_polygon(_support_box_inner, layer)

    _write_field = copy.copy(_support_box_outer)
    _write_field.dx = write_field_x_size
    _write_field.dy = write_field_y_size
    _tmp_write_field = linker_polygon(_write_field, layer)

    _support_outline = copy.copy(_support_box_outer)
    _support_outline.dx = _support_outline.dx + 2 * support_connector_width
    _support_outline.dy = _support_outline.dy + 2 * support_connector_width
    _tmp_support_outline = linker_polygon(_support_outline, layer)

    _box_write_area = gdspy.fast_boolean(_tmp_support_outer, _tmp_support_inner, 'not', max_points=0, layer=layer)

    _tmp_notch_yarray = (_support_box_outer.dy - corner_bend_rad * 2 - support_beam_size) * numpy.linspace(0, 1.0,
                                                                                                           num=num_ysupport)
    for i in range(num_ysupport):
        _tmp_notch_ypos = _ymin * 1.15 + corner_bend_rad + support_beam_size / 2.0 + _tmp_notch_yarray[i]
        _tmp_beam_polygon_L = gdspy.Polygon([
            (_xmin - support_connector_width, _tmp_notch_ypos - support_beam_size / 2.0),
            (_xmin - support_connector_width, _tmp_notch_ypos + support_beam_size / 2.0),
            (_xmin, _tmp_notch_ypos + support_beam_size / 2.0),
            (_xmin, _tmp_notch_ypos - support_beam_size / 2.0),
            (_xmin - support_connector_width, _tmp_notch_ypos - support_beam_size / 2.0)
        ], layer=layer)
        _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(support_beam_size, support_notch_size,
                                                                    _xmin - support_connector_width / 2.0,
                                                                    _tmp_notch_ypos, layer)
        _tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L, _tmp_upper_triangle, 'not', max_points=0,
                                                 layer=layer)
        _tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L, _tmp_lower_triangle, 'not', max_points=0,
                                                 layer=layer)

        _tmp_beam_polygon_R = gdspy.Polygon([
            (_xmax, _tmp_notch_ypos - support_beam_size / 2.0),
            (_xmax, _tmp_notch_ypos + support_beam_size / 2.0),
            (_xmax + support_connector_width, _tmp_notch_ypos + support_beam_size / 2.0),
            (_xmax + support_connector_width, _tmp_notch_ypos - support_beam_size / 2.0),
            (_xmax, _tmp_notch_ypos - support_beam_size / 2.0)
        ], layer=layer)
        _tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(support_beam_size, support_notch_size,
                                                                    _xmax + support_connector_width / 2.0,
                                                                    _tmp_notch_ypos, layer)
        _tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R, _tmp_upper_triangle, 'not', max_points=0,
                                                 layer=layer)
        _tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R, _tmp_lower_triangle, 'not', max_points=0,
                                                 layer=layer)

        _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_L, 'or', max_points=0, layer=layer)
        _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_R, 'or', max_points=0, layer=layer)

    _tmp_notch_xarray = (_support_box_outer.dx - corner_bend_rad * 2 - support_beam_size) * numpy.linspace(0, 1.0,
                                                                                                           num=num_xsupport)
    for i in range(num_xsupport):
        _tmp_notch_xpos = _xmin + corner_bend_rad + support_beam_size / 2.0 + _tmp_notch_xarray[i]
        _tmp_beam_polygon_T = gdspy.Polygon([
            (_tmp_notch_xpos - support_beam_size / 2.0, _ymax),
            (_tmp_notch_xpos - support_beam_size / 2.0, _ymax + support_connector_width),
            (_tmp_notch_xpos + support_beam_size / 2.0, _ymax + support_connector_width),
            (_tmp_notch_xpos + support_beam_size / 2.0, _ymax),
            (_tmp_notch_xpos - support_beam_size / 2.0, _ymax)
        ], layer=layer)
        _tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(support_beam_size, support_notch_size,
                                                                            _tmp_notch_xpos,
                                                                            _ymax + support_connector_width / 2.0,
                                                                            layer)
        _tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T, _tmp_left_triangle, 'not', max_points=0,
                                                 layer=layer)
        _tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T, _tmp_right_triangle, 'not', max_points=0,
                                                 layer=layer)

        _tmp_beam_polygon_B = gdspy.Polygon([
            (_tmp_notch_xpos - support_beam_size / 2.0, _ymin - support_connector_width),
            (_tmp_notch_xpos - support_beam_size / 2.0, _ymin),
            (_tmp_notch_xpos + support_beam_size / 2.0, _ymin),
            (_tmp_notch_xpos + support_beam_size / 2.0, _ymin - support_connector_width),
            (_tmp_notch_xpos - support_beam_size / 2.0, _ymin - support_connector_width)
        ], layer=layer)
        _tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(support_beam_size, support_notch_size,
                                                                            _tmp_notch_xpos,
                                                                            _ymin - support_connector_width / 2.0,
                                                                            layer)
        _tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B, _tmp_left_triangle, 'not', max_points=0,
                                                 layer=layer)
        _tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B, _tmp_right_triangle, 'not', max_points=0,
                                                 layer=layer)

        _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_T, 'or', max_points=0, layer=layer)
        _box_write_area = gdspy.fast_boolean(_box_write_area, _tmp_beam_polygon_B, 'or', max_points=0, layer=layer)

    _support_combined = gdspy.fast_boolean(_tmp_write_field, _tmp_support_outline, 'not', max_points=0, layer=layer)
    _support_combined = gdspy.fast_boolean(_support_combined, _box_write_area, 'or', max_points=0, layer=layer)

    return _support_combined


def write_outer_frame(cell, beamdata, beam_dy_list, outer_box, pattern_number=None, reverse_tone=False, layer=3):
    global num_beam, beam_spacing, edge_offset, linker_edgeoffset, linker_width, text_dist_to_top

    _ymin = beamdata.ypos - beam_dy_list[0] / 2.0 - edge_offset
    _ymax = _ymin + (num_beam - 1) * beam_spacing + numpy.sum(beam_dy_list) + edge_offset * 2

    _tmp_beamdata = copy.copy(beamdata)
    _support_inner = PLTdata()
    _support_inner.xpos = _tmp_beamdata.xpos
    _support_inner.ypos = (_ymin + _ymax) / 2.0
    _support_inner.dx = beamdata.dx + 2 * linker_width + 2 * linker_edgeoffset
    _support_inner.dy = _ymax - _ymin

    _frame_write_area = write_support_region(_support_inner, layer)

    if pattern_number is not None:
        _pattern_number_to_write = write_pattern_number(pattern_number, textheight, (_support_inner.xpos-8e3,
                                                        _support_inner.ypos + _support_inner.dy / 2.0 + 2.5e3))

        # write_pattern_number(pattern_number, size=textheight, position=(xloc, yloc), font_prop=None, tolerance=0.1):


        _frame_write_area = gdspy.fast_boolean(_frame_write_area, _pattern_number_to_write, 'not', max_points=0,
                                               layer=layer)

    if reverse_tone is False:
        bounds = _frame_write_area.get_bounding_box()

        shift = 5000

        box = gdspy.Rectangle(np.array(bounds[0]) - shift, np.array(bounds[1]) + shift, layer=layer)
        bounds_inner = outer_box.get_bounding_box()
        box_inner = gdspy.Rectangle(bounds_inner[0], bounds_inner[1], layer=layer)
        box_inner = box_inner.fillet(corner_bend_rad, points_per_2pi=corner_bend_pts)

        _box_write_area_reverse = gdspy.fast_boolean(box, _frame_write_area, 'not', max_points=0, layer=layer)
        _box_write_area_reverse = gdspy.fast_boolean(_box_write_area_reverse, box_inner, 'not', max_points=0,
                                                     layer=layer)
        cell.add(_box_write_area_reverse)

    else:
        cell.add(_frame_write_area)


class PLTdata:
    def __init__(self):
        self.xpos = None
        self.ypos = None
        self.dx = None
        self.dy = None


def PhC_Writer(param_sweep, end_period=end_period, blank_guides=num_guides, text=text):
    acav = param_sweep[1][0]
    amir = param_sweep[1][1]
    hx = param_sweep[1][5]
    hy = param_sweep[1][6]
    wy = param_sweep[1][2]
    num_cav = param_sweep[1][4]
    num_tap = param_sweep[1][7]
    num_mirr = param_sweep[1][3] * 2
    beams = gdspy.Cell(str(param_sweep[0]))

    rowicolj = PLTdata()
    rowicolj.dx = 10e3
    rowicolj.dy = wy

    rowicolj.xpos = 0
    rowicolj.ypos = 0
    dy_list = write_beam_array(beams, rowicolj)

    circ_rowicolj = PLTdata()
    circ_rowicolj.xpos = rowicolj.xpos - hole_pos_offset
    circ_rowicolj.ypos = rowicolj.ypos
    circ_rowicolj.dx = hx
    circ_rowicolj.dy = hy

    write_hole_2D(beams, rowicolj, circ_rowicolj, dy_list, num_cav_hole=num_cav,
                  num_mir_hole_L=num_mirr / 2,
                  num_mir_hole_R=num_mirr / 2, acav=acav, amir=amir, end_period=end_period, num_end_taper=num_tap,
                  guides=blank_guides,
                  end_taper_L=True,
                  end_taper_R=True, reverse_tone=reverse_tone, edge_overdose=edge_overdose, ring_overdose=ring_overdose,
                  ring_underdose=ring_underdose, edge_underdose=ring_underdose)

    outer_box = write_outer_box(beams, rowicolj, dy_list, grating_spacer=grating_spacer, round_corner=False,
                                direct_write_area=True, write_linker=True,
                                circ_grating=True, reverse_tone=reverse_tone)

    if text == False:
        write_outer_frame(beams, rowicolj, dy_list, outer_box, pattern_number=None, reverse_tone=reverse_tone)
    else:
        write_outer_frame(beams, rowicolj, dy_list, outer_box, pattern_number=param_sweep[0], reverse_tone=reverse_tone)

    return beams


def litho_marks():
    # Define the Heidelberg alignment marks
    heid_mark_factor = 0.9
    heid_mark_pos = mark_pos * heid_mark_factor

    mark_len = 20 * um
    mark_cross = 4 * um
    markersSpacing = 60 * um

    horz = gdspy.Rectangle((-mark_len / 2, -mark_cross / 2), (mark_len / 2, mark_cross / 2), layer=7)
    vert = gdspy.Rectangle((-mark_cross / 2, -mark_len / 2), (mark_cross / 2, mark_len / 2), layer=7)
    heid_alignment_mark.add([horz, vert])

    side = np.linspace(-heid_mark_pos, heid_mark_pos, int(heid_mark_pos * 2 / (markersSpacing + mark_len)))
    corners = [heid_mark_pos, -heid_mark_pos]

    blank = np.empty(len(side))
    for i in corners:
        blank.fill(i)

        side1 = list(zip(blank, side))
        side2 = list(zip(side, blank))

        sides = side1 + side2
        for coord in sides:
            mark = gdspy.CellReference(heid_alignment_mark, coord)
            cell_main.add(mark)


def ebeam_marks(frame):  # Marks defined in the Heidelberg for EBPG alignment
    # # EPBG Alignment Marks
    ebpg_mark_factor = 0.75
    ebpg_size = 20 * um
    corners = [mark_pos * ebpg_mark_factor, -mark_pos * ebpg_mark_factor]

    markHeidelberg = gdspy.Rectangle((0, 0), (ebpg_size, ebpg_size), layer=4)
    EBPG_marks.add(markHeidelberg)

    for i in corners:
        square_pos = gdspy.CellReference(EBPG_marks, (i - ebpg_size / 2, i - ebpg_size / 2))
        frame = gdspy.boolean(frame, square_pos, "not", layer=4)

        square_pos = gdspy.CellReference(EBPG_marks, (-i - ebpg_size / 2, i - ebpg_size / 2))
        frame = gdspy.boolean(frame, square_pos, "not", layer=4)

        # Additional markers to differentiate rotation.
        bot_left_pos = gdspy.CellReference(EBPG_marks,
                                           (-corners[0] - ebpg_size / 2, -corners[0] - ebpg_size / 2 - 600 * um))
        frame = gdspy.boolean(frame, bot_left_pos, "not", layer=4)

        bot_left_pos = gdspy.CellReference(EBPG_marks,
                                           (-corners[0] - ebpg_size / 2 - 600 * um, -corners[0] - ebpg_size / 2))
        frame = gdspy.boolean(frame, bot_left_pos, "not", layer=4)

    cell_main.add(frame)


def supportsMask(name, xpos, ypos, spacing):
    x_frame = write_field_x_size / 2 * 0.9
    y_frame = write_field_y_size / 2 * 0.9

    x_frame = spacing*6 * um /2
    y_frame = ((linker_width/1e3)*2+wy_list[0]/1e3) * um * 12/2


    phc_frame = gdspy.Rectangle((-x_frame, -y_frame), (x_frame, y_frame), layer=4)


    support_mask_inner = 0.82
    spacing = 440
    support_mask_out = gdspy.Rectangle((-spacing * support_mask_inner * um / 2, spacing * support_mask_inner * um / 2),
                                       (spacing * support_mask_inner * um / 2, -spacing * support_mask_inner * um / 2))

    holder = gdspy.Cell('support' + str(name))

    support_mask_frame = gdspy.boolean(support_mask_out, phc_frame, "not", layer=6)
    holder.add(support_mask_frame)

    support_mask_pos = gdspy.CellReference(holder, (x_frame/2, 3300))
    cell_main.add(support_mask_pos)



executeWriting()


if invert == True:
    frame = gdspy.fast_boolean(outer_frame, cell_main, 'not')
    cell_main.add(frame)


def save_file(filename):
    gdspy.write_gds(filename, unit=1.0e-9, precision=1.0e-11)


save_file("2023-11-18-M12.gds")