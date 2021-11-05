+---------+------------------------+---------------------------+------------+
| Branch  | Travis CI              | Code Coverage             | Zenodo DOI |
+---------+------------------------+---------------------------+------------+
| Master  | |Build Status Master|  | |Coverage Status Master|  | |DOI|      |
+---------+------------------------+---------------------------+------------+
| Develop | |Build Status Develop| | |Coverage Status Develop| |            |
+---------+------------------------+---------------------------+------------+

.. NOTE: include statements doesn't work with github README.rst - the first section here is repeated
.. in docs/source/about.rst as well...

.. BEGIN REPETITION ===============================

.. image:: https://www.nersc.no/sites/www.nersc.no/files/images/nansat_logo_transp.png
   :align: right
   :width: 250px
   :target: https://github.com/nansencenter/nansat

**Nansat** is a scientist friendly Python toolbox for processing 2D
satellite earth observation data.

The main **goal** of Nansat is to facilitate:

-  easy development and testing of scientific algorithms,
-  easy analysis of geospatial data, and
-  efficient operational processing.


You can find a detailed description of Nansat in our `paper
<https://openresearchsoftware.metajnl.com/articles/10.5334/jors.120/>`_ published in `Journal of
Open Research Software <https://openresearchsoftware.metajnl.com/>`_ in 2016.

... and you can join the
`mailing list <https://groups.google.com/forum/#!forum/nansat-dev>`_.

We appreciate acknowledgments of Nansat. Please add a reference to the following paper
if you use Nansat in scientific publications:

Korosov A.A., Hansen M.W., Dagestad K.-F., Yamakawa A., Vines A., Riechert M., (2016). Nansat: a
Scientist-Orientated Python Package for Geospatial Data Processing. Journal of Open Research
Software. 4(1), p.e39. DOI: http://doi.org/10.5334/jors.120

.. END REPETITION =================================

Documentation
-------------

You will find complete documentation for Nansat at `Read the Docs`_.

.. _Read the Docs: http://nansat.readthedocs.io/

Contributing
------------

You will find information about contributing to Nansat at `Read the Docs`_.

.. _Read the Docs: http://nansat.readthedocs.io/

Installation
------------

An easy way to install Nansat requirements on any platform is to use Anaconda_ (`download installer <https://conda.io/miniconda.html>`_).

.. _Anaconda: http://docs.continuum.io/anaconda/index

::

    # create environment with key requirements
    conda create -y -n py3nansat gdal numpy pillow netcdf4 scipy
    # activate environment
    conda activate py3nansat
    # install nansat
    pip install nansat
    # launch python
    python

Another option is to use Docker containers (`read about Docker <https://docs.docker.com/>`_).:

::

    # download image with everything pre-installed and launch ipython
    docker run --rm -it -v /path/to/data:/data akorosov/nansat ipython


Example
-------

::

    # download a test file
    wget https://github.com/nansencenter/nansat/raw/develop/nansat/tests/data/stere.tif

.. code:: python

    # import main file opener
    from nansat import Nansat

    # open a test file
    n = Nansat('stere.tif')

    # see file content
    print(n)

    # create RGB with auto-stretched histogram
    n.write_figure('stere_rgb.png', [1,2,3], clim='hist')


Tests
-----

Nansat is outfitted with unittests, which you can use to ensure that all functionality works on your platform.

::

    # install testing packages from conda-forge
    conda install -c conda-forge nose mock

    # run nansat.tests
    nosetests nansat

    # Run all tests including nansat_integration_tests with coverage
    cd <nansat_repository_folder>
    nosetests -w . --with-coverage --cover-package=nansat

Nansat works on both Python 2 and Python 3 but automatic testing on TravisCI is done for Python 3.7 only.
Fore more information see `Read the Docs`_ or notebooks for `Nansat
lectures <https://github.com/nansencenter/nansat-lectures/tree/master/notebooks>`__

.. _Read the Docs: http://nansat.readthedocs.io/

License
-------

The project is licensed under the GNU general public license version 3.

Acknowledgments
----------------

Development is supported by the Research Council of Norway as a part of
`NORMAP <https://normap.nersc.no/>`__ project (grant no. 195397/V30).

.. |Build Status Master| image:: https://github.com/nansencenter/nansat/actions/workflows/tests_build.yml/badge.svg
   :target: https://github.com/nansencenter/nansat/actions/workflows/tests_build.yml
.. |Coverage Status Master| image:: https://coveralls.io/repos/nansencenter/nansat/badge.svg?branch=master&service=github
   :target: https://coveralls.io/github/nansencenter/nansat?branch=master
.. |Build Status Develop| image:: https://github.com/nansencenter/nansat/actions/workflows/tests_build.yml/badge.svg?branch=develop
   :target: https://github.com/nansencenter/nansat/actions/workflows/tests_build.yml
.. |Coverage Status Develop| image:: https://coveralls.io/repos/nansencenter/nansat/badge.svg?branch=develop&service=github
   :target: https://coveralls.io/github/nansencenter/nansat?branch=develop
.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.59998.svg
   :target: https://doi.org/10.5281/zenodo.59998
