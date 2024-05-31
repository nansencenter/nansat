ARG BASE_IMAGE=nansencenter/nansat_base
FROM ${BASE_IMAGE}
# Necessary to access the BASE_IMAGE variable during the build
ARG BASE_IMAGE

ARG NANSAT_RELEASE

COPY . /tmp/nansat/
WORKDIR /tmp/nansat
RUN apt update \
&&  apt install -y --no-install-recommends g++ \
&&  pip install . \
&&  rm -rf /tmp/nansat \
&&  apt autoremove -y \
&&  apt clean \
&&  rm -rf /var/lib/apt/lists/*

WORKDIR /src
