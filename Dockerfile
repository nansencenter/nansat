FROM continuumio/miniconda3

LABEL maintainer="Anton Korosov <anton.korosov@nersc.no>"
LABEL purpose="Python libs for developing and running Nansat"

ENV PYTHONUNBUFFERED=1 \
    PYTHONPATH=/src

RUN apt-get update \
&&  apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
&&  conda config --add channels conda-forge \
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

COPY utilities /tmp/utilities
COPY nansat /tmp/nansat
COPY setup.py /tmp/
WORKDIR /tmp
RUN python setup.py install

WORKDIR /src

