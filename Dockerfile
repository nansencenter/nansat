FROM continuumio/miniconda3 as base

LABEL maintainer="Anton Korosov <anton.korosov@nersc.no>"
LABEL purpose="Nansat and other geospatial Python libs"

ENV PYTHONUNBUFFERED=1 \
    PYTHONPATH=/src \
    MOD44WPATH=/usr/share/MOD44W/

RUN apt-get update \
&&  apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
&&  apt clean \
&&  rm -rf /var/lib/apt/lists/*

RUN conda install setuptools \
&&  conda update conda \
&&  conda config --add channels conda-forge  \
&&  conda install -y \
    ipython \
    ipdb \
    gdal=2.4.2 \
    matplotlib \
    mock \
    netcdf4 \
    nose \
    numpy \
    pillow \
    python-dateutil \
    scipy \
    urllib3 \
    coverage \
    coveralls \
&&  conda remove qt pyqt --force \
&&  conda clean -a -y \
&&  rm /opt/conda/pkgs/* -rf \
&&  pip install pythesint \
&&  python -c 'import pythesint; pythesint.update_all_vocabularies()' \
&&  wget -nc -P /usr/share/MOD44W ftp://ftp.nersc.no/nansat/test_data/MOD44W.tgz \
&&  tar -xzf /usr/share/MOD44W/MOD44W.tgz -C /usr/share/MOD44W/ \
&&  rm -f /usr/share/MOD44W/MOD44W.tgz

WORKDIR /src

FROM base

COPY utilities /tmp/utilities
COPY nansat /tmp/nansat
COPY setup.py /tmp/
WORKDIR /tmp
RUN python setup.py install \
&&  rm -rf /tmp/{utilities,nansat,setup.py}

WORKDIR /src
