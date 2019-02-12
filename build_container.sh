#!/bin/bash

# get UID and GID from input (for linux users only)
DUID=${1-"1000"}
DGID=${2-"1000"}

# build container only with Python libraries in conda
docker build . -t nansat:conda --target conda

# build container with Python libraries and gcc
docker build . -t nansat:dev --target dev

# compile Nansat in current host directory
docker run --rm -it -v `pwd`:/src nansat:dev python setup.py build_ext --inplace

# remove container (if it exists)
docker rm nansat 2> /dev/null
# build container (mount the current directory and nansat) for development
docker create -it --name=nansat \
    -v `pwd`:/src \
    -v `pwd`/nansat:/opt/nansat \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    --security-opt seccomp=unconfined \
    --env=DISPLAY \
    nansat:dev

# build image for production (should be done by image maintainer only)
#docker build . -t nansat:latest --target latest
