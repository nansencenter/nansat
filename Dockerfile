FROM continuumio/miniconda3 as conda
LABEL maintainer="Anton Korosov <anton.korosov@nersc.no>"
LABEL purpose="Python libs for developing and running Nansat"
RUN conda config --add channels conda-forge \
&&  conda update -y conda \
&&  conda install -y \
    ipython \
    ipdb \
    gdal \
    matplotlib \
    mock \
    netcdf4 \
    nose \
    numpy \
    pillow \
    python-dateutil \
    scipy \
    urllib3 \
&&  conda remove qt pyqt --force \
&&  conda clean -a -y \
&&  rm /opt/conda/pkgs/* -rf \
&&  pip install pythesint \
&&  python -c 'import pythesint; pythesint.update_all_vocabularies()'

ENV PYTHONUNBUFFERED=1


FROM nansat:conda as dev
LABEL purpose="Python libs + gcc for developing Nansat"
RUN apt-get update \
&&  apt-get install -y --no-install-recommends \
    build-essential \
    gcc
# compile Nansat inplace and install (using symbolic link)
COPY nansat /tmp/nansat
COPY setup.py /tmp/
WORKDIR /tmp
RUN python setup.py build_ext --inplace \
&&  cp -r nansat /opt/ \
&&  ln -s /opt/nansat /opt/conda/lib/python3.7/site-packages/
WORKDIR /src

# NB! Relevant for Linux users only
# Default values of user and group IDs.
# That should be changed to the IDs of the actual user and group on the host system in the script
# for building the container. For example:
# docker build -t nansat:dev --build-arg DUID=12345 --build-arg DGID=54321 .
# where 12345 and 54321 are user and group ID on a host system
ARG DUID=1000
ARG DGID=1000
RUN groupadd user -g $DGID && useradd -ms /bin/bash -N -u $DUID -g $DGID user
USER user


# create image for production without GCC
FROM nansat:conda as latest
LABEL purpose="Python libs + Nansat"
COPY --from=dev /opt/nansat /opt/conda/lib/python3.7/site-packages/nansat
WORKDIR /src
