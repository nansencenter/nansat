nansat
======

Nansat is a scientist friendly Python toolbox for processing of 2D satellite earth observation data. This is a mirror of the subversion repository found at https://svn.nersc.no/nansat

Requirements
============

See https://svn.nersc.no/nansat/wiki/Nansat/Install/Requirements

Installation
============

See https://svn.nersc.no/nansat/wiki/Nansat/Install/Nansat

or:

## Ubuntu

<pre><code>
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get install gdal-bin
<<<<<<< HEAD
=======
sudo apt-get install libgdal-dev
>>>>>>> develop
sudo apt-get install libgeos-dev
sudo apt-get install libgeos-3.3.8
sudo apt-get install build-essential (?)
sudo apt-get install python-dev
sudo apt-get install gfortran
sudo apt-get install liblapack-dev libatlas-dev
sudo apt-get install freetype*
sudo apt-get install libpng12-dev
sudo apt-get install libgrib-api-dev
sudo apt-get install libgrib-api-tools
sudo apt-get install libgrib2c-dev
sudo apt-get install libhdf5-serial-dev (also takes care of libjpeg-dev, libjpeg-turbo8-dev, libjpeg8-dev)
sudo apt-get install libjasper-dev
sudo apt-get install libnetcdf-dev
sudo apt-get install libproj-dev
sudo apt-get install python-gdal
</code></pre>

Nansat requires basemap, which must be compiled. Check the installation
directory of libgeos (dpkg -L libgeos-3.3.8). Download basemap from
http://sourceforge.net/projects/matplotlib/files/matplotlib-toolkits/ and
compile (in the virtualenv):

<pre><code>
export GEOS_DIR=/usr/lib
python setup.py build 
python setup.py install
</code></pre>

The build may complain that there is no lgeos. In that case make a symbolic
link to libgeos-3.3.8.so called libgeos.so. If you work in a virtualenv, make
symbolic links to gdal and  osgeo in <virtualenv>/lib/python2.7/site-packages.

Finally:

<code><pre>
pip install requirements.txt
</code></pre>
