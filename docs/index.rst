.. Nansat documentation master file, created by
   sphinx-quickstart on Thu Jan 11 09:17:40 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Nansat's documentation!
==================================

.. include:: source/about.rst

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. About Nansat

.. ============

.. * Background scientifics needs

.. * Data formats readable with Nansat

.. * Projects employing Nansat

.. * How to get support

.. * Contributing to Nansat

User documentation
==================

.. toctree::
   :maxdepth: 2

   source/installation.rst

Feature documentation
=====================

Masking land/water in figures
-----------------------------

Using MODIS water-mask product
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To add simple land- or water-masks to your figures, you can use the watermask() method in the main
Nansat class. Download the prepared MODIS 250M water-mask product from our server and add the path
to the directory with this data to an environment variable named MOD44WPATH (e.g.
```MOD44WPATH=/Data/sat/auxdata/mod44w```). The water-mask can be downloaded from here:
`<ftp://ftp.nersc.no/pub/nansat/MOD44W.tgz>`_

Using a digital elevation model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The GMTED2010 datasets are provided by the `U.S. Geological Survey
<https://topotools.cr.usgs.gov/gmted_viewer/>`_. We have prepared a GDAL vrt file and a mapper that
can be used to open the 30 arcseconds Digital Elevation Model (DEM) with Nansat.

Global 30 Arc-Second Elevation (GTOPO30)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


Developer documentation
=======================

.. toctree::
   :maxdepth: 2

   source/conventions.rst
   source/create_mapper.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
