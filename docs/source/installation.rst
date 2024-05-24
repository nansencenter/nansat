Installation
============

Quickstart
----------
The fastest way to install nansat:

* Install `Miniconda <https://conda.io/miniconda.html>`_ on your platform of choice.

.. code-block:: bash

    # create environment with key requirements
    conda create -y -n py3nansat gdal=2.4.2 numpy pillow netcdf4 scipy
    # activate environment
    source activate py3nansat
    # install nansat
    pip install nansat

Nansat is now installed.
For more details and other methods of installing Nansat, see below.

Requirements
------------

Nansat requires the following packages:

* Python 3.7 or higher
* `Numpy <http://www.numpy.org/>`_ >=1.11.3
* `GDAL <http://www.gdal.org>`_ >=2.2.3
* `Pillow <https://python-pillow.github.io/>`_ >=4.0.0
* `netCDF4 <https://github.com/Unidata/netcdf4-python>`_ >=1.3.1
* `py-thesaurus-interface <https://github.com/nansencenter/py-thesaurus-interface>`_

The following packages are optional:

* `Scipy <http://scipy.org/>`_ 0.18.1

  * Some mappers will not work without scipy. E.g. *sentinel1_l1*

* `Matplotlib <http://matplotlib.org/>`_ >=2.1.1,<3.9

  * matplotlib is required for Nansat methods *digitize_points()* and *crop_interactive()*

* `Basemap <http://matplotlib.org/basemap/>`_ >=1.0.8

  * basemap is required in *write_domain_map()*

The most tricky to compile yourself is GDAL and Basemap. But one can find pre-built binaries
available for different platforms. We recommend to install all dependencies with Conda, from the
conda-forge channel. See instructions on this below.

Installing Requirements
-----------------------

You have three main options on how to install the requirements. These are described in the
following three sections.


Install dependencies from Anaconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the recommended approach for installing dependencies.

* Download `Miniconda <https://conda.io/miniconda.html>`_ for your platform of choice.

* Install Miniconda

  * When you install Miniconda on Windows, you will get a new app called "Anaconda Prompt".
    Run this to access the conda installation.

  * On Linux/OS X use a regular terminal and make sure PATH is set to contain the installation
    directory as explained by the installer.

* Run the following three commands:

  * *conda create -n nansat Python=3.6*

  * *source activate nansat*

    * On windows you would omit 'source' and just run *'activate nansat'*

  * *conda install --yes -c conda-forge pythesint numpy scipy=0.18.1 matplotlib basemap netcdf4 gdal
    pillow urllib3*

Install Pre-built Binaries
^^^^^^^^^^^^^^^^^^^^^^^^^^
One can find pre-built binaries available for different platforms. We do not have an overview over
all the possible repositories where you can find binaries. But if you e.g. are on Ubuntu, the
following procedure can be used to install dependencies with *apt* and *pip*.

.. code-block:: bash

   sudo apt install -y virtualenv libgdal-dev gdal-bin python3
   cd
   virtualenv nansat_env
   source ~/nansat_env/bin/activate
   pip install pythesint pillow netcdf4 urllib3 requests "gdal==$(gdal-config --version)" numpy \
    scipy 'matplotlib<3.9' setuptools basemap


Compile and Build Yourself
^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have the technical expertise to build all dependencies, and need to do it yourself, feel
free to do so. If you need some aid, we would recommend you to look at how the corresponding
`conda-forge feedstocks <https://github.com/conda-forge/>`_ have been built.

Installing Nansat
-----------------

Install with pip (preferred method)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the following command:

::

  pip install nansat

Nansat will then be added to your site-packages and can be used like any regular Python package.


Install Nansat from source
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to install Nansat from source, you first need to install all requirements.
Then proceed with one of the following methods


Install from git repository
"""""""""""""""""""""""""""

git clone the master (most stable) or develop (cutting edge) branch, and install:

.. code-block:: bash

   git clone https://github.com/nansencenter/nansat.git
   cd nansat
   git checkout master (or a specific tag or branch)
   pip install .

Nansat will then be added to your site-packages and can be used like any regular Python package.

..
  Install from Anaconda
  ^^^^^^^^^^^^^^^^^^^^^
  TODO: Add instructions about installing from Anaconda when conda-forge has accepted the feedstock
  request. Basicall copy what's in Install dependencies from Anaconda but install only nansat.
  Also update the link to "simplest way to install Nansat" in basic info.

Special install for Nansat Developers
"""""""""""""""""""""""""""""""""""""
If you are working directly on the Nansat source, you need to install Nansat in the following way.

Git clone the branch you are working on, and run the following command from the nansat directory:

::

  pip install -e .

The pixel functions C module is then compiled but no code is copied to site-packages and no linking
is performed. Make sure to follow the `Nansat conventions <conventions.html>`_ if you want to
contribute to Nansat.

In addition to the regular dependencies, developers also need to install nose and mock. This can
easily be done with

::

  pip install mock


Use a self-provisioned Virtual Machine
--------------------------------------

Another option to install Nansat in a controlled environment is to use a virtual machine. Configuration
for `Vagrant <https://www.vagrantup.com/>`_ and `Ansible <https://www.ansible.com/>`_ that brings up and
provision a `VirtualBox <https://www.virtualbox.org/>`_ machine is provided in Nansat repository. To start
the machine you need to install Vagrant and VirtualBox on your computer; clone or download the nansat
source code; and start the machine:

::

  # download nansat source code
  git clone https://github.com/nansencenter/nansat.git
  cd nansat

  #start virtual machine
  vagrant up

That's it! The virtual machine will be started and all software will be installed automatically. To start using Nansat you need to log in to the virtual machine and start Python from the conda environment:

::

  vagrant ssh
  source activate py3nansat
  python


Use Docker
----------
Docker is a platform for developers and sysadmins to develop, deploy, and run applications with
containers (`Get started with Docker <https://docs.docker.com/get-started/>`_). We have developed
an image that containes compiled Nansat an a number of Python libraries needed for development
and running of Nansat. A user can start using the production version of Nansat Docker image:

::

    docker run --rm -it -v /path/to/data:/data nansencenter/nansat ipython

This will mound directory /path/to/data on your host to the directory /data in the container
and launch IPython where Nansat is available.

For developing Nansat you needs access to the code both from the container (to run it)
and from the host (to edit it). For this purpose you should clone Nansat repository and
do the following steps:
1. Build pixelfunctions inplace

::

    docker run --rm -it -v `pwd`:/src nansencenter/nansat_base pip install -e .

2. Run container with mounting of the current directory into /src. In this case Python
will use Nansat from /src/nansat (the directory shared between host and container):

::

    # launch Python with Nansat in container
    docker run --rm -it -v `pwd`:/src nansencenter/nansat_base python

    # ...or run tests
    docker run --rm -it -v `pwd`:/src nansencenter/nansat_base python -m unittest discover nansat.tests

Alternatively you can run the script *build_container.sh*. The script will build the image with
Python libraries from Anaconda, compile the Nansat code inplace and create a
container for running Nansat. You can then start container:

::

    docker start -i nansat
    # and run tests:
    (base) root@d1625f2ce873:~# python -m unittest discover nansat.tests
