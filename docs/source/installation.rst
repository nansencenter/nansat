Installation
============

Install Nansat from source
--------------------------

Get code and dependencies

* Install the `Requirements`_ 
* `Download the code <https://github.com/nansencenter/nansat/releases>`_

Requirements
^^^^^^^^^^^^

TODO: update!

Nansat depends on the following packages:
* Python 2.6 or higher
* `Numpy <http://www.numpy.org/>`_
* `Scipy <http://scipy.org/SciPy>`_
* `Matplotlib <http://matplotlib.org/>`_
* `Basemap <http://matplotlib.org/basemap/>`_
* `Python Imaging Library (PIL) <http://www.pythonware.com/products/pil/>`_
* `GDAL <http://www.gdal.org>`_
* `Pillow <https://python-pillow.github.io/>`_
* `py-thesaurus-interface <https://github.com/nansencenter/nersc-metadata>`_

The most tricky is to install GDAL and Basemap. One can find pre-built binaries
available for different platforms. However, we recommend to download the source and
build GDAL on your computer for reading, e.g., NetCDF or HDF5 files. Some hints are
given on the `GDAL web-site <http://trac.osgeo.org/gdal/wiki/BuildHints>`_.

Another option is to use a virtual machine managed by Virtualbox and provisioned
using Vagrant and Ansible. We provide 
`configurations for virtual machines <https://github.com/nansencenter/geo-spaas-vagrant>`_ 
for learning and development of Nansat.

Here are **instructions** how to initialize a virtual machine and materials for
`Nansat lectures <https://github.com/nansencenter/nansat-lectures>`_.

For regular users
^^^^^^^^^^^^^^^^^

::

  python setup.py install

Nansat will then be added to your site-packages and can be used like any regular Python package.

For developers (working directly with the source)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  python setup.py build_ext --inplace

The pixel functions C module is then compiled but no code is copied to site-packages and no linking
is performed.

Install from PyPi
-----------------

TODO: Add instructions when 

Install from Anaconda
---------------------

