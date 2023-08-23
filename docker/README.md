# Running MAGIQ in Docker

1. Download and install [Docker Desktop](https://www.docker.com/)
1. Navigate to the location of this directory and create the docker compose environment.
   ```bash
   docker-compose up -d
   ```
1. Navigate to http://localhost:8080/vnc.html in your browser. Click the 'Connect' button. 
1.To stop the containers and remove them:
   ```bash
   docker-compose down
   ```

# Debugging issues

- [Docker Dashboard](https://docs.docker.com/desktop/use-desktop/) provides a way to examine the stdout/stderr printouts of the MAGIQ tools. This is where errors will show up.
- You can also [run the terminal](https://www.docker.com/blog/integrated-terminal-for-running-containers-extended-integration-with-containerd-and-more-in-docker-desktop-4-12/) inside the container.


# Configuration

By default, the root directory of MAGIQ is volume-mapped to `/app` inside the container, e.g. this readme can be found in `/app/docker/README.md`.

To map your data directory, you can add another volume in the `docker-compose.yml` file, in the `volumes` section of the `magiq` service. A commented-out example is provided.

To use a different MATLAB network licence, change the value of `MLM_LICENSE_FILE` environment variable in `docker-compose.yml`. See [Matlab Docker image documentation](https://hub.docker.com/r/mathworks/matlab) for more info.