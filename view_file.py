import gdspy
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def create_str_table(headers, values, delimiter="|", margin=" "):
    lengths = [max([len(str(d[i])) for d in [headers] + values]) for i in range(len(headers))]
    get_row = lambda data: delimiter + margin + (margin + delimiter + margin).join(
        [str(d).rjust(lengths[i]) for i, d in enumerate(data)]) + margin + delimiter
    ret_str = get_row(headers)
    ret_str += "\n" + delimiter + delimiter.join([(len(margin) * 2 + length) * "-" for length in lengths]) \
        + delimiter + "\n"
    ret_str += "\n".join([get_row(v) for v in values])
    return ret_str


def get_cell_dict(cell):
    return [cell.name, len(cell.polygons), len(cell.paths), len(cell.labels), len(cell.references)]


def generate_cell_report(gds_file):
    """ Creates a table of cell makeups. """
    return create_str_table(
        headers=["name", "polygons", "paths", "labels", "references"],
        values=[get_cell_dict(cell) for cell in gds_file.cells.values()]
    )


def get_mpl_polygons_from_gds_polygons(gds_polygons):
    return [
        Polygon(coords, closed=True, edgecolor='black', facecolor='none') 
        for coords in gds_polygons.polygons
    ]


def get_mpl_polygons_from_gds_cell(cell):  # make thi better. TODO
    cell_geometries = cell.get_polygons(by_spec=True)
    polygons = []
    for coords, nested_coords_list in cell_geometries.items():
        vertices = nested_coords_list[0]
        polygons.append(
            Polygon(vertices, closed=True, edgecolor='black', facecolor='none')
        )
    return polygons


def get_mpl_polygons_from_gds_cells(gds_file):
    max_len = max([len(name) for name in gds_file.cells.keys()])
    matplotlib_cells = []
    for i, (name, cell) in enumerate(gds_file.cells.items()):
        print(f"Graphing cell {name}...{' ' * (max_len - len(name))} ", end="")
        try:
            matplotlib_cells += get_mpl_polygons_from_gds_cell(cell)
            print("converted!")
        except Exception as e:
            print("failed!")
    return matplotlib_cells


def create_plot(mpl_polygons):
    _, ax = plt.subplots()
    for mpl_polygon in mpl_polygons:
        ax.add_patch(mpl_polygon)
    ax.autoscale(enable=True) 
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.show()


def main():
    gds_filename = r"C:\LocalCodingProjects\KargEtcher\2023-11-18-M12.gds"
    gds_file = gdspy.GdsLibrary().read_gds(gds_filename)

    print(generate_cell_report(gds_file))
    mpl_polygons = get_mpl_polygons_from_gds_cells(gds_file)
    create_plot(mpl_polygons)


if __name__ == "__main__":

    # import pya

    # # Create a new layout view
    # layout = pya.Layout()
    # view = pya.LayoutView(layout)

    # # Open the GDS file
    # gds_file_path = "path/to/your/file.gds"
    # gds_file_path = r"C:\LocalCodingProjects\KargEtcher\2023-11-18-M12.gds"
    # layout.read(gds_file_path)

    # # Display the layout
    # view.show()

    # import pya

    # # Open the GDS file
    # gds_file_path = "path/to/your/file.gds"
    # layout = pya.InputGDSFile(gds_file_path).read()

    # # Create a viewer and display the layout
    # viewer = pya.Viewer()
    # viewer.create_layout(1).add_cell(layout.cell(layout.top_cell()))
    # viewer.zoom_fit()
    # viewer.show()

    main()
