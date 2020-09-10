ARG BASE_IMAGE
FROM ${BASE_IMAGE}

COPY utilities /tmp/utilities
COPY nansat /tmp/nansat
COPY setup.py /tmp/
WORKDIR /tmp
RUN apt update \
&&  apt install -y --no-install-recommends g++ \
&&  python setup.py install \
&&  rm -rf /tmp/{utilities,nansat,setup.py} \
&&  apt remove -y gcc \
&&  apt autoremove -y \
&&  apt clean \
&&  rm -rf /var/lib/apt/lists/*


WORKDIR /src