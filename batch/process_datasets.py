import argparse
import os

import yaml

from apps_processing import AppsProcessing as apps
from magiqdataclasses import DatFile


def file_path(f):
    if not os.path.exists(f) and os.path.isfile(f):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist or is not a file".format(f))
    return f


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Process datasets")

    p.add_argument("--config", dest="config", type=file_path,
                   default=f"{os.path.dirname(os.path.abspath(__file__))}{os.sep}config.yml")

    command_args = p.parse_args()
    with open(command_args.config, 'r') as file:
        config = yaml.safe_load(file)

    study_directories = [d for d in os.listdir(config['output_directory']) if os.path.isdir(os.path.join(config['output_directory'], d))]
    for study in study_directories:
        print(f"---- Processing {study}")
        study_directory = os.path.join(config['output_directory'], study)

        # conversion
        conversion_config=config['process']['conversion']
        suppressed_directory = os.path.join(study_directory, 'sup')
        converted_directory = os.path.join(suppressed_directory, 'converted')
        os.makedirs(converted_directory, exist_ok=True)
        plot_path = os.path.join(
                study_directory,
                conversion_config['save_plot']
            ) if conversion_config['save_plot'] else None
        if plot_path:
            os.makedirs(os.path.dirname(plot_path), exist_ok=True)
        sup_file, uns_file = apps.run_bruker_conversion(
            out_name=converted_directory + os.sep,
            sup_file_path=suppressed_directory,
            ref_file_path=os.path.join(study_directory, 'uns'),
            baseline_correction=conversion_config['baseline_correction'],
            post_processing=conversion_config['post_processing'],
            quality_points=conversion_config['quality_points'],
            save_plot=plot_path,
        )

        # water removal
        wr_config = config['process']['water_removal']
        wr_input = os.path.join(
            converted_directory, f"corr_sup{'_bc' if conversion_config['baseline_correction'] else ''}_"
                                 f"{conversion_config['post_processing']}{conversion_config['quality_points']}.dat")

        input_dat = DatFile(wr_input)
        plot_path = os.path.join(
                study_directory,
                wr_config['save_plot']
            ) if wr_config['save_plot'] else None
        dat_hsvd, dat_wr = apps.run_water_removal(
            input_dat=input_dat,
            hsvd_points=wr_config['hsvd_points'],
            hsvd_ratio=wr_config['hsvd_ratio'],
            hsvd_components=wr_config['hsvd_components'],
            frequency_range_xmin=wr_config['frequency_range_xmin'],
            frequency_range_xmax=wr_config['frequency_range_xmax'],
            save_plot=plot_path,
        )

        apps.save_water_removal(input_dat, dat_wr)

        print(f"---- Finished processing {study}")

