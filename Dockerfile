FROM nansencenter/nansat_base

COPY utilities /tmp/utilities
COPY nansat /tmp/nansat
COPY setup.py /tmp/
WORKDIR /tmp
RUN python setup.py install \
    &&  rm -rf /tmp/{utilities,nansat,setup.py}

WORKDIR /src