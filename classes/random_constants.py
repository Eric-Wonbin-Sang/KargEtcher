# Membrane Parameters
GAP = 10  # 10um gap between the membrane and support structure
SUPPORT_WIDTH = 10  # The width of the support for the suspended layer
SUPPORT_PINCH_WIDTH = 3  # The width of the pinch point for stamping
SUPPORT_PINCH_LENGTH = 2  # The length of the pinch point region
SUPPORT_PITCH = 0.4  # Ratio of support structure to exposed area

# Membrane parameters
LENGTH = 150
HEIGHT = 200

DEVICE_ROW_COUNT = 6 # Number of times to write each row of membranes
DEVICE_COL_COUNT = 6

# Unit for all parameters: nm
# Parameters for exposure
WRITE_FIELD_X_SIZE = 50e3  # frame surrounding the phc
WRITE_FIELD_Y_SIZE = 100e3

# Parameters for nanobeam
PHC_GROUP_COUNT = 11
BEAM_SPACING = 6.5e3
EDGE_OFFSET = 5e3
CORNER_BEND_RAD = 3e3
CORNER_BEND_PTS = 41

# Parameters for the linker
LINKER_EDGE_OFFSET = 0
LINKER_WIDTH = 6e3
LINKER_NOTCH_SIZE = 500
LINKER_BEAM_SIZE = 500
LINKER_X_NOTCHES = 0

SAFETY = 60 # Space between membranes
SPACING = LENGTH + SAFETY + SUPPORT_WIDTH * 2
SPACING_Y = HEIGHT + SAFETY + SUPPORT_WIDTH * 2

PAIRS = 1

# photonic_components                      layout_management
# from classes.random_constants import (  | from classes.random_constants import (
#     CELL_MAIN,                          |     # CELL_MAIN,                       | C#
#     GRATING_CELL,                       |     # GRATING_CELL,                    | G#
#     REFLECTOR_CELL,                     |     # REFLECTOR_CELL,                  | R#
#     GRATING_TAPER_CELL,                 |     # GRATING_TAPER_CELL,              | G#
#     REFLECTOR_TAPER_CELL,               |     # REFLECTOR_TAPER_CELL,            | R#
#     # DEVICE_POS,                       |     # DEVICE_POS,                      | ##
#     # SUPPORT_MASK,                     |     # SUPPORT_MASK,                    | ##
#     # PHC_FRAME,                        |     # PHC_FRAME,                       | ##
#     # EBPG_MARKS,                       |     # EBPG_MARKS,                      | ##
#     # HEID_ALIGNMENT_MARK,              |     # HEID_ALIGNMENT_MARK,             | ##
#     DEVS,                               |     # DEVS,                            | D#
#     DEVICE_ROW_COUNT,                   |     DEVICE_ROW_COUNT,                  | DD
#     DEVICE_COL_COUNT,                   |     DEVICE_COL_COUNT,                  | DD
#     WRITE_FIELD_X_SIZE,                 |     WRITE_FIELD_X_SIZE,                | WW
#     WRITE_FIELD_Y_SIZE,                 |     WRITE_FIELD_Y_SIZE,                | WW
#     SHOT_PITCH,                         |     # SHOT_PITCH,                      | S#
#     PHC_GROUP_COUNT,                    |     PHC_GROUP_COUNT,                   | PP
#     BEAM_SPACING,                       |     BEAM_SPACING,                      | BB
#     NUM_CIRC_POINTS,                    |     # NUM_CIRC_POINTS,                 | N#
#     HOLE_POS_OFFSET,                    |     # HOLE_POS_OFFSET,                 | H#
#     EDGE_OFFSET,                        |     EDGE_OFFSET,                       | EE
#     CORNER_BEND_RAD,                    |     CORNER_BEND_RAD,                   | CC
#     CORNER_BEND_PTS,                    |     CORNER_BEND_PTS,                   | CC
#     # TAPER_NEFF,                       |     # TAPER_NEFF,                      | ##
#     SPACER,                             |     # SPACER,                          | S#
#     UNDER_EDGE_SIZE,                    |     # UNDER_EDGE_SIZE,                 | U#
#     LINKER_EDGE_OFFSET,                 |     LINKER_EDGE_OFFSET,                | LL
#     LINKER_WIDTH,                       |     LINKER_WIDTH,                      | LL
#     LINKER_CONNECTOR_WIDTH,             |     # LINKER_CONNECTOR_WIDTH,          | L#
#     LINKER_NOTCH_SIZE,                  |     # LINKER_NOTCH_SIZE,               | L#
#     LINKER_BEAM_SIZE,                   |     # LINKER_BEAM_SIZE,                | L#
#     LINKER_X_NOTCHES,                   |     # LINKER_X_NOTCHES,                | L#
#     LINKER_Y_NOTCHES,                   |     # LINKER_Y_NOTCHES,                | L#
#     # SUPPORT_CONNECTOR_WIDTH,          |     SUPPORT_CONNECTOR_WIDTH,           | #S
#     # SUPPORT_NOTCH_SIZE,               |     SUPPORT_NOTCH_SIZE,                | #S
#     # SUPPORT_BEAM_SIZE,                |     SUPPORT_BEAM_SIZE,                 | #S
#     # NUM_X_SUPPORT,                    |     NUM_X_SUPPORT,                     | #N
#     # NUM_Y_SUPPORT,                    |     NUM_Y_SUPPORT,                     | #N
#     # PERIOD,                           |     # PERIOD,                          | ##
#     # DUTY_CYCLE,                       |     # DUTY_CYCLE,                      | ##
#     GRATING_LINE_WIDTH,                 |     # GRATING_LINE_WIDTH,              | G#
#     GRATING_SPACING,                    |     # GRATING_SPACING,                 | G#
#     # CIRC_GRATING_SUPPORT,             |     # CIRC_GRATING_SUPPORT,            | ##
#     NUM_CIRC_GRATING_POINTS,            |     # NUM_CIRC_GRATING_POINTS,         | N#
#     NUM_GRATING,                        |     # NUM_GRATING,                     | N#
#     CIRC_GRATING_BASE,                  |     # CIRC_GRATING_BASE,               | C#
#     GRATING_ANGLE,                      |     # GRATING_ANGLE,                   | G#
#     ODD_SUPPORT_ANGLE,                  |     # ODD_SUPPORT_ANGLE,               | O#
#     EVEN_SUPPORT_ANGLE,                 |     # EVEN_SUPPORT_ANGLE,              | E#
#     SUPPORT_ANGLE_WIDTH,                |     # SUPPORT_ANGLE_WIDTH,             | S#
#     # K,                                |     # K,                               | ##
#     GAP,                                |     # GAP,                             | G#
#     SUPPORT_WIDTH,                      |     # SUPPORT_WIDTH,                   | S#
#     SUPPORT_PINCH_WIDTH,                |     # SUPPORT_PINCH_WIDTH,             | S#
#     SUPPORT_PINCH_LENGTH,               |     # SUPPORT_PINCH_LENGTH,            | S#
#     SUPPORT_PITCH,                      |     # SUPPORT_PITCH,                   | S#
#     LENGTH,                             |     # LENGTH,                          | L#
#     HEIGHT,                             |     # HEIGHT,                          | H#
#     # SAFETY,                           |     # SAFETY,                          | ##
#     SPACING,                            |     SPACING,                           | SS
#     SPACING_Y,                          |     SPACING_Y,                         | SS
#     PAIRS,                              |     PAIRS,                             | PP
#     # COLUMN_SIZE,                      |     # COLUMN_SIZE,                     | ##
# )                                       | )

# photonic_components                      layout_management
# from classes.random_constants import (  | from classes.random_constants import (
#     # SUPPORT_CONNECTOR_WIDTH,          |     SUPPORT_CONNECTOR_WIDTH,           | #S
#     # SUPPORT_NOTCH_SIZE,               |     SUPPORT_NOTCH_SIZE,                | #S
#     # SUPPORT_BEAM_SIZE,                |     SUPPORT_BEAM_SIZE,                 | #S
#     # NUM_X_SUPPORT,                    |     NUM_X_SUPPORT,                     | #N
#     # NUM_Y_SUPPORT,                    |     NUM_Y_SUPPORT,                     | #N

#     CELL_MAIN,                          |     # CELL_MAIN,                       | C#
#     GRATING_CELL,                       |     # GRATING_CELL,                    | G#
#     REFLECTOR_CELL,                     |     # REFLECTOR_CELL,                  | R#
#     GRATING_TAPER_CELL,                 |     # GRATING_TAPER_CELL,              | G#
#     REFLECTOR_TAPER_CELL,               |     # REFLECTOR_TAPER_CELL,            | R#
#     DEVS,                               |     # DEVS,                            | D#
#     SHOT_PITCH,                         |     # SHOT_PITCH,                      | S#
#     NUM_CIRC_POINTS,                    |     # NUM_CIRC_POINTS,                 | N#
#     HOLE_POS_OFFSET,                    |     # HOLE_POS_OFFSET,                 | H#
#     SPACER,                             |     # SPACER,                          | S#
#     UNDER_EDGE_SIZE,                    |     # UNDER_EDGE_SIZE,                 | U#
#     LINKER_CONNECTOR_WIDTH,             |     # LINKER_CONNECTOR_WIDTH,          | L#
#     LINKER_NOTCH_SIZE,                  |     # LINKER_NOTCH_SIZE,               | L#
#     LINKER_BEAM_SIZE,                   |     # LINKER_BEAM_SIZE,                | L#
#     LINKER_X_NOTCHES,                   |     # LINKER_X_NOTCHES,                | L#
#     LINKER_Y_NOTCHES,                   |     # LINKER_Y_NOTCHES,                | L#
#     GRATING_LINE_WIDTH,                 |     # GRATING_LINE_WIDTH,              | G#
#     GRATING_SPACING,                    |     # GRATING_SPACING,                 | G#
#     NUM_CIRC_GRATING_POINTS,            |     # NUM_CIRC_GRATING_POINTS,         | N#
#     NUM_GRATING,                        |     # NUM_GRATING,                     | N#
#     CIRC_GRATING_BASE,                  |     # CIRC_GRATING_BASE,               | C#
#     GRATING_ANGLE,                      |     # GRATING_ANGLE,                   | G#
#     ODD_SUPPORT_ANGLE,                  |     # ODD_SUPPORT_ANGLE,               | O#
#     EVEN_SUPPORT_ANGLE,                 |     # EVEN_SUPPORT_ANGLE,              | E#
#     SUPPORT_ANGLE_WIDTH,                |     # SUPPORT_ANGLE_WIDTH,             | S#
#     GAP,                                |     # GAP,                             | G#
#     SUPPORT_WIDTH,                      |     # SUPPORT_WIDTH,                   | S#
#     SUPPORT_PINCH_WIDTH,                |     # SUPPORT_PINCH_WIDTH,             | S#
#     SUPPORT_PINCH_LENGTH,               |     # SUPPORT_PINCH_LENGTH,            | S#
#     SUPPORT_PITCH,                      |     # SUPPORT_PITCH,                   | S#
#     LENGTH,                             |     # LENGTH,                          | L#
#     HEIGHT,                             |     # HEIGHT,                          | H#

#     DEVICE_ROW_COUNT,                   |     DEVICE_ROW_COUNT,                  | DD
#     DEVICE_COL_COUNT,                   |     DEVICE_COL_COUNT,                  | DD
#     WRITE_FIELD_X_SIZE,                 |     WRITE_FIELD_X_SIZE,                | WW
#     WRITE_FIELD_Y_SIZE,                 |     WRITE_FIELD_Y_SIZE,                | WW
#     PHC_GROUP_COUNT,                    |     PHC_GROUP_COUNT,                   | PP
#     BEAM_SPACING,                       |     BEAM_SPACING,                      | BB
#     EDGE_OFFSET,                        |     EDGE_OFFSET,                       | EE
#     CORNER_BEND_RAD,                    |     CORNER_BEND_RAD,                   | CC
#     CORNER_BEND_PTS,                    |     CORNER_BEND_PTS,                   | CC
#     LINKER_EDGE_OFFSET,                 |     LINKER_EDGE_OFFSET,                | LL
#     LINKER_WIDTH,                       |     LINKER_WIDTH,                      | LL
#     SPACING,                            |     SPACING,                           | SS
#     SPACING_Y,                          |     SPACING_Y,                         | SS
#     PAIRS,                              |     PAIRS,                             | PP
# )                                       | )