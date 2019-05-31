Releasing Nansat
==================

General release procedure
-------------------------

In setup.py, Make sure the version is correct, by checking MAJOR, MINOR and MICRO, and that 'dev'
is removed from VERSION.

Releases should be made from the master branch. **Make sure that Nansat works on all operating
systems (Windows, Linux and OS X), with all combinations of python (py27, py35, py36) before
releasing.**

When you are ready to release a new version, create a tag for the release and push it.
The following creates a tag for version 0.6.17 and pushes it.

::

  git tag "v0.6.17"
  git push origin --tags

Go to https://github.com/nansencenter/nansat/releases and click "Draft a new release".

Select the tag version you just created, create a release title, for consistency, use the format of
Nansat-0.6.17

Write some release notes that describes the changes done in this release.


Releasing on PiPy
-----------------

First, wait until Nansat passes all tests on Travi-CI, Appveyor and Coverals. Then execute:

::

   conda create -n release_nansat -c conda-forge -y python pythesint scipy basemap netcdf4 gdal pillow mock nose urllib3 twine
   source activate release_nansat
   python setup.py sdist
   # Check the dist file that was just created
   ls dist
   # Should be a file on this format 'nansat-1.0.20.tar.gz'
   twine upload dist/nansat-1.0.20.tar.gz

Packaging documentation is found at `PyPA Packaging and Distributing Projects
<https://packaging.python.org/tutorials/distributing-packages/>`_

To avoid having to enter password when uploading, you can set $HOME/.pypirc as described in the
above link.

Releasing on Anaconda
---------------------

We are releasing Nansat through the conda-forge channel on Anaconda. First, wait until Nansat passes
all tests on Travi-CI, Appveyor and Coverals. Then execute:

::

   # install (or update) conda-smithy
   conda install -n root -c conda-forge conda-smithy
   git clone git@github.com:conda-forge/nansat-feedstock.git
   cd nansat-feedstock
   conda smithy rerender -c auto
   git push

Information how to use Conda-Smithy can be found at `The tool for managing conda-forge feedstocks
<https://github.com/conda-forge/conda-smithy>`_


Releasing on Docker
-------------------
To build an image with stable version of Nansat you need to build the image and push it to Docker Hub:

::
    # build the base image with conda
    docker build . -t nansat:conda --target conda
    # build the image for compiling Nansat
    docker build . -t nansat:dev --target dev
    # build the image for distributing Nansat
    docker build . -t nansat:latest --target latest
    # find the ide of the nansat:latest image (e.g. bb38976d03cf)
    docker images
    # tag the image (bb38976d03cf is the id of the nansat:latest image)
    docker tag bb38976d03cf akorosov/nansat:latest
    # authenticate on Docker Hub
    docker login --username=yourhubusername --email=youremail@company.com
    # push the image to Docker Hub
    docker push akorosov/nansat:latest
