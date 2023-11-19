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

# import gdspy
# import numpy

# gds_filename = r"C:\LocalCodingProjects\KargEtcher\2023-11-10-O1.gds"

# # Open the GDSII file
# gds_file = gdspy.GdsLibrary().read_gds(gds_filename)

# # Specify the cell name containing the content you want to visualize
# cell_name = 'main_cell'

# # Retrieve the cell by its name
# cell = gds_file.cells[cell_name]

# # Extract the cell's geometries
# cell_geometries = cell.get_polygons(by_spec=True)

# import matplotlib.pyplot as plt
# from matplotlib.patches import Polygon

# def gdspy_to_matplotlib(polygons):
#     matplotlib_polygons = []
#     for coordinate, data in polygons.items():
#         # print(coordinate, numpy.asarray(polygon), numpy.asarray(polygon).shape)
#         print(coordinate)
#         # print(data)
#     matplotlib_polygons.append(Polygon(list(polygons.keys()), fill=True))
#     return matplotlib_polygons

# matplotlib_cells = gdspy_to_matplotlib(cell_geometries)

# fig, ax = plt.subplots()

# for polygon in matplotlib_cells:
#     ax.add_patch(polygon)

# ax.set_aspect('equal', adjustable='box')
# # plt.axis('off')  # Optional: Turn off axis if not needed
# plt.show()




import gdspy
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MatplotlibPolygon

# Create a Polygon with initial vertices
initial_vertices = [(0, 0), (1, 0), (1, 1), (0, 1)]
polygon = gdspy.Polygon(initial_vertices)

# # Modify the shape of the Polygon using set_xy
# new_vertices = [(0, 0), (2, 0), (2, 2), (0, 2)]
# polygon.set_xy(new_vertices)

# Create a Matplotlib figure and axis
fig, ax = plt.subplots()

# Convert the GDSPY Polygon to a Matplotlib Polygon and add it to the axis
matplotlib_polygon = MatplotlibPolygon(initial_vertices, closed=True, edgecolor='black', facecolor='none')
ax.add_patch(matplotlib_polygon)

# Set aspect ratio and display the plot
ax.set_aspect('equal', adjustable='box')
# plt.axis('off')  # Optional: Turn off axis if not needed
plt.show()