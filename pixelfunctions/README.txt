:Author:    Antonio Valentino
:Contact:   a_valentino@users.sf.net
:Date:      2011-05-13
:Copyright: This document has been placed in the public domain.

GDAL pixfun plugin
==================

This package provides:

* the implementation of a set of GDALDerivedPixelFunc(s) to be used with
  source raster band of virtual GDAL datasets
* a fake GDAL driver to register pixel functions

.. note:: using the plugin mechanism is a hack aimed to enable python users
          to use pixel functions without C++ coding

To use the plugin just build the gdal_PIXFUN.so shared object::

    make

and set the GDAL_DRIVER_PATH accordingly::

    export GDAL_DRIVER_PATH=<PATH_TO_GDAL_PIXFUN_SO>:$GDAL_DRIVER_PATH

.. seealso:: http://lists.osgeo.org/pipermail/gdal-dev/2011-May/028737.html
