ARG BASE_IMAGE
FROM ${BASE_IMAGE}
# Necessary to access the BASE_IMAGE variable during the build
ARG BASE_IMAGE

ARG NANSAT_RELEASE

COPY utilities /tmp/utilities
COPY nansat /tmp/nansat
COPY setup.py /tmp/
WORKDIR /tmp
RUN apt update \
&&  apt install -y --no-install-recommends g++ \
&&  python setup.py install \
&&  rm -rf /tmp/{utilities,nansat,setup.py} \
&&  if [ -n "`echo $BASE_IMAGE | grep slim`" ];then apt remove -y gcc;fi \
&&  apt autoremove -y \
&&  apt clean \
&&  rm -rf /var/lib/apt/lists/*

WORKDIR /src
