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



gds_filename = r"C:\LocalCodingProjects\KargEtcher\2023-11-10-O1.gds"
gds_cell = gdspy.GdsLibrary().read_gds(gds_filename)

# Create a new figure
fig, ax = plt.subplots()

# Define a function to convert GDS elements to polygons
def gds_element_to_polygons(element):
    polygons = []
    if isinstance(element, gdspy.Polygon):
        polygons.append(element)
    elif isinstance(element, gdspy.CellReference):
        cell = element.get_bounding_box()
        polygons.extend(gds_element_to_polygons(cell))
    elif isinstance(element, gdspy.CellArray):
        for cell_ref in element.references:
            cell = cell_ref.get_bounding_box()
            polygons.extend(gds_element_to_polygons(cell))
    return polygons

# Get the polygons from the cell
polygons = []
for element in gds_cell.cells:
    polygons.extend(gds_element_to_polygons(element))

# Plot each polygon
for polygon in polygons:
    x, y = polygon.points.T
    ax.fill(x, y, edgecolor='k', facecolor='none')

# Set the aspect ratio to be equal
ax.set_aspect('equal')

# Show the plot or save it to an image file
# plt.show()  # Uncomment this line to display the image
plt.savefig('layout.png')  # Save the image to a file

# Close the figure to release resources
plt.close()
