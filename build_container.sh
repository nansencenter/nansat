#!/bin/bash
# build container only with Python libraries in conda
docker build . -t nansat

# compile Nansat in current host directory
docker run --rm -it -v `pwd`:/src nansat pip install -e /src

# remove container (if it exists)
docker rm nansat 2> /dev/null
# build container (mount the current directory and nansat) for development
docker create -it --name=nansat \
    -v `pwd`:/src \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    --security-opt seccomp=unconfined \
    --env=DISPLAY \
    nansat
