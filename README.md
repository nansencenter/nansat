nansat
======

Nansat is a scientist friendly Python toolbox for processing 2D satellite earth observation data.

The main **goal** of Nansat is to facilitate:

* easy development and testing of scientific algorithms,
* easy analysis of geospatial data, and
* efficient operational processing.

Downlodad
=========

https://github.com/nansencenter/nansat/wiki/Download

Requirements
============

See https://github.com/nansencenter/nansat/wiki/Required-libs

or:

## Ubuntu

<pre><code>
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get install gdal-bin
sudo apt-get install libgdal-dev
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

<code><pre>
pip install requirements.txt
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


