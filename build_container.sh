#!/bin/bash

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
    nansat:dev

# build image for production
#docker build . -t nansat:latest --target latest
