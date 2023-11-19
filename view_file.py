import gdspy
import matplotlib.pyplot as plt


# gds_filename = r"C:\LocalCodingProjects\KargEtcher\2023-11-10-O1.gds"

# gds_cell = gdspy.GdsLibrary().read_gds(gds_filename)

# # Create a new figure
# fig, ax = plt.subplots()

# # Get the polygons from the cell
# polygons = gds_cell.get_polygons()

# # Plot each polygon
# for polygon in polygons:
#     x, y = polygon.T  # Transpose the polygon to get x and y coordinates
#     ax.fill(x, y, edgecolor='k', facecolor='none')

# # Set the aspect ratio to be equal
# ax.set_aspect('equal')

# # Show the plot or save it to an image file
# # plt.show()  # Uncomment this line to display the image
# plt.savefig('layout.png')  # Save the image to a file

# # Close the figure to release resources
# plt.close()



# gds_filename = r"C:\LocalCodingProjects\KargEtcher\2023-11-10-O1.gds"
# gds_cell = gdspy.GdsLibrary().read_gds(gds_filename)

# # Create a new figure
# fig, ax = plt.subplots()

# # Define a function to convert GDS elements to polygons
# def gds_element_to_polygons(element):
#     polygons = []
#     if isinstance(element, gdspy.Polygon):
#         polygons.append(element)
#     elif isinstance(element, gdspy.CellReference):
#         cell = element.get_bounding_box()
#         polygons.extend(gds_element_to_polygons(cell))
#     elif isinstance(element, gdspy.CellArray):
#         for cell_ref in element.references:
#             cell = cell_ref.get_bounding_box()
#             polygons.extend(gds_element_to_polygons(cell))
#     return polygons

# # Get the polygons from the cell
# polygons = []
# for element in gds_cell.cells:
#     polygons.extend(gds_element_to_polygons(element))

# # Plot each polygon
# for polygon in polygons:
#     x, y = polygon.points.T
#     ax.fill(x, y, edgecolor='k', facecolor='none')

# # Set the aspect ratio to be equal
# ax.set_aspect('equal')

# # Show the plot or save it to an image file
# # plt.show()  # Uncomment this line to display the image
# plt.savefig('layout.png')  # Save the image to a file

# # Close the figure to release resources
# plt.close()

import gdspy
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def get_mpl_cell_from_gds_cell(cell):
    cell_geometries = cell.get_polygons(by_spec=True)
    return Polygon(list(cell_geometries.keys()), fill=True)


def create_plot(mpl_polygons):
    _, ax = plt.subplots()
    for mpl_polygon in mpl_polygons:
        ax.add_patch(mpl_polygon)
    # ax.autoscale(enable=True) 
    ax.set_aspect('equal', adjustable='box')
    # plt.tight_layout()
    plt.show()


def main():
    gds_filename = r"C:\LocalCodingProjects\KargEtcher\2023-11-10-O1.gds"
    gds_file = gdspy.GdsLibrary().read_gds(gds_filename)

    matplotlib_cells = []
    for name, cell in gds_file.cells.items():
        print(f"{name.rjust(20)}: {cell}...", end="")
        try:
            matplotlib_cells.append(
                get_mpl_cell_from_gds_cell(cell)
            )
            print("good!")
        except Exception as e:
            print("bad!")
        # print("-" * 40)
    
    create_plot(matplotlib_cells)

if __name__ == "__main__":
    main()


"""
name                        polygons paths labels references passed
main_cell                   10       0     0      1071       good
grating                     0        0     0      0          bad
reflector                   0        0     0      0          bad
grating_taper               0        0     0      0          bad
reflector_taper             0        0     0      0          bad
device_pos                  0        0     0      0          bad
support_mask                0        0     0      0          bad
phc_frame                   0        0     0      0          bad
ebpg_marks                  9        0     0      0          good
heid_alignment mark         18       0     0      0          good
membrane_width_sweep        0        0     0      0          bad
membrane_width_sweep_nohole 0        0     0      9          good
support_structure           4        0     0      0          good
0                           13       0     0      0          good
membrane00                  2        0     0      36         good
support0                    1        0     0      0          good
1                           13       0     0      0          good
membrane10                  2        0     0      36         good
support1                    1        0     0      0          good
2                           13       0     0      0          good
membrane20                  2        0     0      36         good
support2                    1        0     0      0          good
3                           13       0     0      0          good
membrane01                  2        0     0      36         good
support3                    1        0     0      0          good
4                           13       0     0      0          good
membrane11                  2        0     0      36         good
support4                    1        0     0      0          good
5                           13       0     0      0          good
membrane21                  2        0     0      36         good
support5                    1        0     0      0          good
6                           13       0     0      0          good
membrane02                  2        0     0      36         good
support6                    1        0     0      0          good
7                           13       0     0      0          good
membrane12                  2        0     0      36         good
support7                    1        0     0      0          good
8                           13       0     0      0          good
membrane22                  2        0     0      36         good
support8                    1        0     0      0          good

"""


# import gdspy
# import matplotlib.pyplot as plt
# from matplotlib.patches import Polygon as MatplotlibPolygon

# # Create a Polygon with initial vertices
# initial_vertices = [(0, 0), (1, 0), (1, 1), (0, 1)]
# polygon = gdspy.Polygon(initial_vertices)

# # # Modify the shape of the Polygon using set_xy
# # new_vertices = [(0, 0), (2, 0), (2, 2), (0, 2)]
# # polygon.set_xy(new_vertices)

# # Create a Matplotlib figure and axis
# fig, ax = plt.subplots()

# # Convert the GDSPY Polygon to a Matplotlib Polygon and add it to the axis
# matplotlib_polygon = MatplotlibPolygon(initial_vertices, closed=True, edgecolor='black', facecolor='none')
# ax.add_patch(matplotlib_polygon)

# # Set aspect ratio and display the plot
# ax.set_aspect('equal', adjustable='box')
# # plt.axis('off')  # Optional: Turn off axis if not needed
# plt.show()