from dataclasses import dataclass, asdict
from datetime import datetime
from functools import cached_property
import os
import gdspy
from classes.generic import logger, Component, config_dataclass
import numpy

from classes.generic import Unit, logger
from classes.photonic_components import Unit, Hole, NanoBeam, GratingCoupler, Phc
from classes.layout_management import WorkArea


class PhcGenerator:

    @config_dataclass
    class HoleVariationConfig:
        hole_width_variations             : list  # previously hx_list
        hole_height_variations            : list  # previously hy_list
        hole_spacing_variations           : list  # previously acav
        hole_group_padding_variations     : list  # previously amir
        hole_modification_func_variations : any   # previously 
        hole_count_variations             : list  # previously mirror_list
        hole_count_with_acav_spacing      : list  # number of holes with acav spacing (previously cavity_list)
        taper_hole_count                  : list 
        taper_hole_modification_func      : any 

    @config_dataclass
    class BeamVariationConfig:
        width_variations  : list  # previously wy_list
        height_variations : list  # previously

    # @config_dataclass
    # class TaperedHoleVariationConfig:
    #     hole_width_variations             : list  # previously hx_list
    #     hole_height_variations            : list  # previously hy_list
    #     hole_spacing_variations           : list  # previously acav
    #     hole_group_padding_variations     : list  # previously amir
    #     hole_modification_func_variations : any   # previously 
    #     hole_count_variations             : list  # previously mirror_list
    #     hole_count_with_acav_spacing      : list  # number of holes with acav spacing (previously cavity_list)
    #     taper_hole_count                  : list 
    #     taper_hole_modification_func      : any 
    
    # @config_dataclass
    # class EdgeHoleVariationConfig:
    #     hole_width_variations             : list  # previously hx_list
    #     hole_height_variations            : list  # previously hy_list
    #     hole_spacing_variations           : list  # previously acav
    #     hole_group_padding_variations     : list  # previously amir
    #     hole_modification_func_variations : any   # previously 
    #     hole_count_variations             : list  # previously mirror_list
    #     hole_count_with_acav_spacing      : list  # number of holes with acav spacing (previously cavity_list)
    #     taper_hole_count                  : list 
    #     taper_hole_modification_func      : any 
    
    # @config_dataclass
    # class MiddleHoleVariationConfig:
    #     hole_width_variations             : list  # previously hx_list
    #     hole_height_variations            : list  # previously hy_list
    #     hole_spacing_variations           : list  # previously acav
    #     hole_group_padding_variations     : list  # previously amir
    #     hole_modification_func_variations : any   # previously 
    #     hole_count_variations             : list  # previously mirror_list
    #     hole_count_with_acav_spacing      : list  # number of holes with acav spacing (previously cavity_list)
    #     taper_hole_count                  : list 
    #     taper_hole_modification_func      : any 

    variation_config_to_config_classes = {
        HoleVariationConfig: Hole.Config,
        BeamVariationConfig: NanoBeam.Config,
    }

    @classmethod
    def get_phc_configs(cls, phc_generator_hole_variation_config, phc_generator_beam_variation_config):

        num_rows = 6  # Number of times to write each row of membranes
        num_cols = 6

        dev_list = numpy.arange(0, num_rows * num_cols, 1)

        device_name_to_params = zip(
            dev_list,
            numpy.array(
                numpy.meshgrid(
                    *phc_generator_hole_variation_config.list_type_attrib_values \
                        + phc_generator_beam_variation_config.list_type_attrib_values
                )
            ).T.reshape(
                -1, len(phc_generator_hole_variation_config.list_type_attrib_names \
                        + phc_generator_beam_variation_config.list_type_attrib_names)
            ).astype(int)
        )

        def get_config_dict(variation_config, component_config_class, params):
            data = {attrib_name: value for attrib_name, value in zip(variation_config.list_type_attrib_names, params)}
            data = {
                attrib_name: data[attrib_name] 
                if attrib_name in data else 
                getattr(variation_config, attrib_name) 
                for attrib_name in component_config_class.attribute_names
            }
            return component_config_class(**data)

        for device_id, params in device_name_to_params:
            print(get_config_dict(phc_generator_hole_variation_config, Hole.Config, params[:len(phc_generator_hole_variation_config.list_type_attrib_names)]))
            print(get_config_dict(phc_generator_beam_variation_config, NanoBeam.Config, params[len(phc_generator_beam_variation_config.list_type_attrib_names):]))
            # hole_config = Hole.Config
            # print(Hole.Config(**{name: params[i] for i, name in enumerate(Hole.Config.attribute_names)}))
            # print(Beam.Config(**{name: params[i + len(Hole.Config.attribute_names) - 1] for i, name in enumerate(Beam.Config.attribute_names)}))

            # print(device_id)
            # print(params)
        # exit()

    @classmethod
    def create_phc_variations(cls, phc_generator_hole_variation_config, phc_generator_beam_variation_config):
        cls.get_phc_configs(phc_generator_hole_variation_config, phc_generator_beam_variation_config)

        DEVICE_ROW_COUNT = 6
        DEVICE_COL_COUNT = 6

        device_id = 0
        phcs = []
        for col_i in range(DEVICE_COL_COUNT):
            if col_i == 2:
                return phcs
            for row_i in range(DEVICE_ROW_COUNT):
                phcs.append(
                    Phc(
                        device_id=device_id,
                        row_i=row_i, 
                        col_i=col_i, 
                        beam=NanoBeam(
                            width=60 * Unit.um.value, 
                            height= 30 * Unit.um.value
                        ),
                        grating_coupler=None,  # GratingCoupler()
                    )
                )
                device_id += 1
        return phcs


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

    def __init__(self, name, work_area, phc_generator_hole_variation_config, phc_generator_beam_variation_config) -> None:
        self.name = name
        self.work_area = work_area
        self.phc_generator_hole_variation_config = phc_generator_hole_variation_config
        self.phc_generator_beam_variation_config = phc_generator_beam_variation_config

        self.phc_variants = PhcGenerator.create_phc_variations(self.phc_generator_hole_variation_config, self.phc_generator_beam_variation_config)

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
    creator = Creator(
        name="test",
        work_area=WorkArea(
            x_length=6 * Unit.mm.value,  # The size of our chip - assumed to be 6mm x 6mm
            y_length=6 * Unit.mm.value,  #  - reduce slightly in case of variations
            layer=4,
            add_position_ids=True,
            add_litho_marks=True,
        ),
        phc_generator_hole_variation_config=PhcGenerator.HoleVariationConfig(
            hole_width_variations             = [63, 66, 69, 73, 76, 79],
            hole_height_variations            = [180, 184, 188, 192, 196, 200],
            hole_spacing_variations           = [155.5],
            hole_group_padding_variations     = [172.1],
            hole_modification_func_variations = lambda x: x ** 1,
            hole_count_variations             = [10],
            hole_count_with_acav_spacing      = [16],
            taper_hole_count                  = [0],
            taper_hole_modification_func      = None,
        ),
        phc_generator_beam_variation_config=PhcGenerator.BeamVariationConfig(
            width_variations             = [470],
            height_variations            = [239048],
        )
    )
    creator.save_file()


if __name__ == "__main__":
    main()
