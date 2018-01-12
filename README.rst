|Build Status| |Coverage Status| |DOI|

.. image:: https://www.nersc.no/sites/www.nersc.no/files/images/nansat_logo_transp.png
   :align: right
   :width: 250px
   :target: https://github.com/nansencenter/nansat

About
-----

**Nansat** is a scientist friendly Python toolbox for processing 2D
satellite earth observation data.

The main **goal** of Nansat is to facilitate:

-  easy development and testing of scientific algorithms,
-  easy analysis of geospatial data, and
-  efficient operational processing.

We appreciate acknowledments of Nansat. Please add "The image analysis
was performed with the open-source Nansat
(https://github.com/nansencenter/nansat) python package" (or equivalent)
if you use Nansat in scientific publications.

Installation
------------

The easiest way to install Nansat on a Linux machine is to use
`anaconda <http://docs.continuum.io/anaconda/index>`__

::

    # download the latest version of miniconda
    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh

    # make it executable
    chmod +x miniconda.sh

    # install miniconda virtual environment
    ./miniconda.sh -b -f -p $HOME/miniconda

    # activate the miniconda environment
    export PATH=$HOME/miniconda/bin/:$PATH

    # update miniconda packages and metadata
    conda update -q --yes conda

    # install all requirements from conda-forge channel
    conda install -q --yes -c conda-forge gdal numpy pillow netcdf4 cfunits python-dateutil pythesint nose

    # finally install Nansat
    pip install https://github.com/nansencenter/nansat/archive/master.tar.gz

    # run nansat.tests
    nosetests nansat

    # Run all tests including nansat_integration_tests with coverage
    cd <nansat_repository_folder>
    nosetests -w . --with-coverage --cover-package=nansat

Fore more information see
`Install-Nansat <https://github.com/nansencenter/nansat/wiki/Install-Nansat>`__
section or use pre-configure virtual machines as explained on
`Nansat-lectures <https://github.com/nansencenter/nansat-lectures>`__

Usage
-----

.. code:: python

    # download a test file
    !wget https://github.com/nansencenter/nansat/raw/develop/nansat/tests/data/stere.tif

    # import main file opener
    from nansat import Nansat

    # open a test file
    n = Nansat('stere.tif')

    # see file content
    print n

    # view file footpring
    n.write_map('stere.footpring.png')

    # create RGB with auto-stretched histogram
    n.write_figure('stere_rgb.png', [1,2,3], clim='hist')

Fore more information see
`Tutorial <https://github.com/nansencenter/nansat/wiki/Tutorial>`__ or
notebooks for `Nansat
lectures <https://github.com/nansencenter/nansat-lectures/tree/master/notebooks>`__

License
-------

The project is licensed under the GNU general public license version 3.

Acknowledgments
----------------

Development is supported by the Research Council of Norway as a part of
`NORMAP <https://normap.nersc.no/>`__ project (grant no. 195397/V30).

.. |Build Status| image:: https://travis-ci.org/nansencenter/nansat.svg?branch=master
   :target: https://travis-ci.org/nansencenter/nansat
.. |Coverage Status| image:: https://coveralls.io/repos/nansencenter/nansat/badge.svg?branch=master
   :target: https://coveralls.io/r/nansencenter/nansat
.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.59998.svg
   :target: https://doi.org/10.5281/zenodo.59998
