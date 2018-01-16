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

TODO: update the following text....

To add land- or water-masks to your figures, you will need to use the watermask() method in the main
Nansat class. Download the prepared MODIS 250M water-mask product from our server and add the path
to the directory with this data to an environment variable named MOD44WPATH (e.g.
```MOD44WPATH=/Data/sat/auxdata/mod44w```). The water-mask can be downloaded from here:
`<ftp://ftp.nersc.no/pub/nansat/MOD44W.tgz>`_

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
