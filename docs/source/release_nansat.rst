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

..
  TODO: add this when nansat is released on conda-forge
  Save draft (wait with publishing this release until the tag builds correctly on conda-forge)

Releasing on PiPy
-----------------

..
  TODO: add a note on waiting with releasing it on PyPi until the tag is built on conda-forge.

Packaging documentation is found at `PyPA Packaging and Distributing Projects 
<https://packaging.python.org/tutorials/distributing-packages/>`_

To avoid having to enter password when uploading, you can set $HOME/.pypirc as described in the
above link.

.. code-block:: bash

   conda create -n release_nansat -y python=3.6
   source activate release_nansat
   conda install -y -c conda-forge pythesint scipy=0.18.1 basemap netcdf4 gdal pillow
   pip install mock nose urllib3 twine
   python setup.py sdist
   # Check the dist file that was just created
   ls dist
   # Should be a file on this format 'nansat-0.6.17.tar.gz'
   twine upload dist/nansat-0.6.17.tar.gz

..
  Releasing on Anaconda
  ---------------------
  
  We will release Nansat through the conda-forge channel on Anaconda. 
  TODO: Update these instructions when the following pull-request has been merged and a feedstock
  for Nansat is created.
  https://github.com/conda-forge/staged-recipes/pull/4818