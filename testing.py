import gdspy

um = 1000
chip_size = 6000 * um  # The size of our chip - assumed to be 6mm x 6mm - reduce slightly in case of variations
length = 220
safety = 200  # Space between membranes
pairs = 1
support_width = 10  # The width of the support for the suspended layer
spacing = length + safety + support_width * 2
devs = 3
num_cols = devs
num_rows = devs  # Number of times to write each row of membranes
column_size = num_rows * spacing


outer_frame = gdspy.Rectangle((-chip_size / 2, chip_size / 2), (chip_size / 2, -chip_size / 2), layer=4)
inner_frame = gdspy.Rectangle((-spacing * um * num_rows / 2, spacing * um * num_cols * pairs / 2),
                              (spacing * um * num_rows / 2, -spacing * um * num_cols * pairs / 2), layer=4)
frame = gdspy.boolean(outer_frame, inner_frame, "not", layer=4)





import gdspy
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MatplotlibPolygon


def create_mpl_polygon_from_polygon_set(polygon_set):
    vertices = polygon_set.polygons[0]
    return MatplotlibPolygon(vertices, closed=True, edgecolor='black', facecolor='none')


def create_plot(mpl_polygons):
    _, ax = plt.subplots()
    for mpl_polygon in mpl_polygons:
        ax.add_patch(mpl_polygon)
    ax.autoscale(enable=True) 
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.show()


def main():
    matplotlib_polygons = [
        create_mpl_polygon_from_polygon_set(frame)
    ]
    create_plot(matplotlib_polygons)


if __name__ == "__main__":
    main()
