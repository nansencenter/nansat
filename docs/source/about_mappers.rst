About Nansat mappers
============================

Concept 
-------

GDAL can read most satellite EO and NetCDF-CF compliant raster data relevant for Earth sciences.
However, GDAL does not attach meaning to the contained information, i.e., it does not
specify what kind of data is contained in a given band, e.g., if it is *water-leaving-radiance*
used for monitoring of water quality. Nansat provides a mapping between
geophysical variables of known meaning and the raster bands provided in the "Datasets" returned by GDAL.

The modules within the mappers package each contain a class defining the mapping between the bands
returned from GDAL and metadata vocabularies provided via the `py-thesaurus-interface
<https://github.com/nansencenter/py-thesaurus-interface>`_ package. For example, the simplest mapper for Meris level-1 data
explicitly states that the first 15 bands are *upwelling_spectral_radiance* at 15 wavelengths and
that the 16th band contains quality flags. The description of these bands require some compulsory
metadata attributes to be defined in the mapper. These attributes follow certain given conventions:

* `NetCDF-CF <http://cfconventions.org/>`_
* `NASA GCMD keywords <https://earthdata.nasa.gov/about/gcmd/global-change-master-directory-gcmd-keywords>`_
* `nersc-vocabularies <https://github.com/nansencenter/nersc-vocabularies>`_

By using these conventions, the mappers thus attach unambiguous geophysical meaning to variables
following the given standards. This allows the user to open and use a geo-refenced raster dataset
with Nansat without depending on detailed apriori knowledge about the origin or type of the data.

Workflow
--------

When we open a file with Nansat:

.. code-block:: python

   #!/usr/bin/env python
   n = Nansat(filename)

these steps follow:

* The Nansat constructor calls gdal.Open(filename) to open the file with GDAL, and returns a GDAL Dataset with a list of available raster bands
* The Nansat constructor loops through available mappers and parses the Dataset to the mapper

  * Each mapper checks if the input Dataset is appropriate for the mapper, i.e., if the format, the metadata and the set of bands in the Dataset corresponds to what is expected in the mapper

    * If the Dataset is not valid, the mapper silently fails and the next mapper is tested
    * If the Dataset fits the mapper:

        * the mapper creates a  `GDAL VRT file <http://www.gdal.org/gdal_vrttut.html>`_ with georeference and raster bands corresponding to the "well known variables" in `nersc-vocabularies <https://github.com/nansencenter/nersc-vocabularies>`_ and adds respective metadata to each band (standard_name, units, etc).

* The mapper object, which is an instance of the VRT-class, is then returned to the Nansat instance as an attribute named ``vrt`` (``Nansat.vrt``)


The VRT has the following properties:

* we can use any available GDAL API functions, e.g., warping or exporting
* it contains georeferencing recognised by GDAL
* we can add PixelFunctions for, e.g., calculation of *speed* given two vector components of wind or current
* still it contains only Raster Bands with metadata which correspond to any of the NANSAT "Well Known Variables"

The Dataset may, e.g., be subsetted, reprojected, merged,
etc., by simply modifying the VRT-file, either automatically by the GDAL high level
applications/functions, or with NANSAT-specific Python logic. An important benefit of this approach
is that we employ the *lazy processing concept* in GDAL.

.. note::

   No processing or file reading/writing is performed before it is needed. 
   
The VRT file defines a set of operations in xml format. When information is needed, data is
extracted as numpy arrays for further processing or plotting. As such, we basically use the GDAL
Datamodel and do not need to design our own.

Technical details
-----------------

* The VRT-file is stored in memory using GDAL VSI-approach
* The VRT-class is a wrapper around the VRT-file. It has methods for generating, modifying, copying and other operations with VRT-files. VRT-class uses both GDAL methods and direct writing for modifying the VRT-file.
* Each mapper inherits the VRT-class.

Where to put new mappers?
-------------------------

If you have created a new mapper, you can either submit a pull request for inclusion in the nansat
mappers package, or you can make a namespace package to let nansat discover your mapper
automatically. This is done the following way:

1. Create a directory called *nansat_mappers* within a directory on your $PYTHONPATH 
2. Inside *nansat_mappers*, create the file *__init__.py* with the following lines:

.. code-block:: python

   # __init__.py
   from pkgutil import extend_path
   __path__ = extend_path(__path__, __name__)

3. Add your mapper module (the filename should start with *mapper_* and end with *.py*) to the *nansat_mappers* folder 
4. Reload your shell and start Python
5. Nansat should now find you mapper. 
   
Note that user defined mappers have higher priority than standard mappers.

Required metadata added in the mappers
--------------------------------------

* TODO: add list of required metadata

Adding mapper tests
-------------------

TODO: add documentation about how to write mapper tests


