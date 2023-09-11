# Running MAGIQ in Docker

1. Download and install [Docker Desktop](https://www.docker.com/)


# Configuration

By default, the root directory of MAGIQ is volume-mapped to `/app` inside the container, e.g. this readme can be found in `/app/docker/README.md`.

To map your data directory, you can add another volume in the `docker-compose.yml` file, in the `volumes` section of the `magiq` service. A commented-out example is provided in `docker-compose.yml`.

To use a different MATLAB network licence, change the value of `MLM_LICENSE_FILE` environment variable in `docker-compose.yml`. See [Matlab Docker image documentation](https://hub.docker.com/r/mathworks/matlab) for more info.

## Running the GUI 

1. Navigate to the location of this directory and create the docker compose environment.
   ```bash
   docker-compose up -d
   ```
1. Navigate to http://localhost:8080/vnc.html in your browser. Click the 'Connect' button. 

1. To stop the containers and remove them:
   ```bash
   docker-compose down
   ```

## Running the batch scripts

1. The `batch` directory contains a `config.yml` file. This is used to configure the downloading and preprocessing.
 - set the `output_directory` to a directory that has been volume-mapped into the container (see `Configuration` section)
- set the `download.user` to a user allowed to query & retrieve the data from the DICOM server
 - Modify `download.study_filters` to select the correct studies to retrieve
 - Modify `download.spectroscopy_filters` and `download.structural_filters` to uniquely select the correct series.

1. If you haven't run the MAGIQ GUI using `docker-compose`, you'll need to build the Docker image before using it below.
   ```bash
   cd /path/to/MAGIQ/docker # substitute root directory of this repo for /path/to/MAGIQ
   docker build .
   ```

1. To download and preprocess the data, run the `download_datasets.py` script:
   ```bash
   cd /path/to/MAGIQ/docker # substitute root directory of this repo for /path/to/MAGIQ
   read -s -p "Password of DICOM user: " USER_PASSWORD && printf "\n" && docker run --rm -v $(dirname "$(pwd)"):/app -e USER_PASSWORD=${USER_PASSWORD} docker-magiq:latest python -u /app/batch/download_datasets.py
   ```

You can confirm that the datasets have been downloaded by looking at the location of the volume-mapped `output_directory` on the host machine.

1. To process the data (conversion, water removal), run the `process_datasets.py` script:

   ```bash
   docker run --rm -v $(dirname "$(pwd)"):/app docker-magiq:latest python -u batch/process_datasets.py
   ```

# Debugging issues

- [Docker Dashboard](https://docs.docker.com/desktop/use-desktop/) provides a way to examine the stdout/stderr printouts of the MAGIQ tools. This is where errors will show up.
- You can also [run the terminal](https://www.docker.com/blog/integrated-terminal-for-running-containers-extended-integration-with-containerd-and-more-in-docker-desktop-4-12/) inside the container.
