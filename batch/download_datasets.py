import argparse
import os
import shutil
import subprocess
from urllib.parse import urlparse

from keycloak import KeycloakOpenID
from dicomweb_client.api import DICOMwebClient

import pydicom

import zipfile
import io
import yaml
import dirsync

def deflate_dataset(data, output_directory):
    with zipfile.ZipFile(io.BytesIO(data), "r") as zf:
        for zi in zf.infolist():
            try:
                deflate_dataset(zf.read(zi), output_directory)
            except zipfile.BadZipFile:
                zf.extract(zi, output_directory)


def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


def file_path(f):
    if not os.path.exists(f) and os.path.isfile(f):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist or is not a file".format(f))
    return f


def valid_url(arg):
    url = urlparse(arg)
    if all((url.scheme, url.netloc)):
        return arg
    raise argparse.ArgumentTypeError("Invalid URL")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Download DICOM datasets")

    p.add_argument("--config", dest="config", type=file_path,
                   default=f"{os.path.dirname(os.path.abspath(__file__))}{os.sep}config.yml")

    command_args = p.parse_args()
    with open(command_args.config, 'r') as file:
        config = yaml.safe_load(file)
        dl_config = config['download']

    keycloak_openid = KeycloakOpenID(
        server_url=dl_config['keycloak_server'],
        client_id=dl_config['keycloak_client'],
        realm_name=dl_config['keycloak_realm'],
    )

    password = os.environ["USER_PASSWORD"]
    auth = keycloak_openid.token(username=dl_config['user'], password=password)

    client = DICOMwebClient(
        url="https://dicom.cfmm.uwo.ca/dcm4chee-arc/aets/CFMMARC/rs",
        headers={"Authorization": f"Bearer {auth['access_token']}"},
    )

    studies = [
        pydicom.dataset.Dataset.from_json(s)
        for s in client.search_for_studies(
            # "20220201-20230901"
            search_filters=dl_config['study_filters']
        )
    ]

    studies.sort(key=lambda x: x.StudyDate)

    print("Downloading and pre-processing DICOM datasets.")
    for study in studies:
        spectroscopy_series = [
            pydicom.dataset.Dataset.from_json(s)
            for s in client.search_for_series(
                study.StudyInstanceUID, search_filters=dl_config['spectroscopy_filters']
            )
        ]
        structural_series = [
            pydicom.dataset.Dataset.from_json(s)
            for s in client.search_for_series(
                study.StudyInstanceUID, search_filters=dl_config['structural_filters']
            )
        ]
        if len(spectroscopy_series) == 1:

            spectroscopy_instances = [
                pydicom.dataset.Dataset.from_json(i)
                for i in client.search_for_instances(
                    study_instance_uid=study.StudyInstanceUID,
                    series_instance_uid=spectroscopy_series[0].SeriesInstanceUID,
                )
            ]
            spectroscopy_data = client.retrieve_bulkdata(
                url=f"{spectroscopy_instances[0].RetrieveURL}/bulkdata"
            )
            study_id = f"{study.StudyDate}_{study.PatientName}"
            print(f'--- Processing {study_id}')
            study_directory = f"{config['output_directory']}{os.sep}{study_id}"
            structural_directory = f"{study_directory}{os.sep}structural"
            os.makedirs(structural_directory, exist_ok=True)

            # extract spectroscopy data
            deflate_dataset(spectroscopy_data[0], study_directory)

            # download structural DICOM
            structural_data = client.retrieve_series(
                study_instance_uid=study.StudyInstanceUID,
                series_instance_uid=structural_series[0].SeriesInstanceUID,
            )

            for dataset in structural_data:
                dataset.save_as(
                    f"{structural_directory}{os.sep}{dataset.SOPInstanceUID}.dcm"
                )

            print(
                f"{study.StudyInstanceUID}: {study.StudyDate}, {len(spectroscopy_series)}, {len(structural_series)}"
            )

            if config['preprocess']['convert_to_nifti']:
                if not os.path.exists(os.path.join(study_directory, study_id + '.nii.gz')):
                    cmd = [
                        'dcm2niix', '-z', 'y', '-o', study_directory,
                        '-f', study_id, structural_directory
                    ]
                    print("Converting DICOM to Nifti")
                    print(subprocess.check_output(cmd).decode())
                else:
                    print('Nifti file exists, skipping')
            if config['preprocess']['rename']:
                print("Renaming")
                for source, destination in config['preprocess']['rename'].items():
                    try:
                        dirsync.sync(
                            f"{study_directory}{os.sep}{source}",
                            f"{study_directory}{os.sep}{destination}",
                            action="sync",
                            create=True,
                        )
                        shutil.rmtree(f"{study_directory}{os.sep}{source}")
                    except FileNotFoundError:
                        print(f"ERROR: did not find {study_directory}{os.sep}{source}")
        else:
            assert (
                    len(spectroscopy_series) == 0
            ), "Found study with >1 spectroscopy series"
    print("Finished downloading and pre-processing DICOM datasets.")