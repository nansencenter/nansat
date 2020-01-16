Differentiating between land and water
---------------------------------------

To add simple land- or water-masks to your figures, you can use the watermask() method in the main
Nansat class. Download the prepared `MODIS 250M water-mask product
<ftp://ftp.nersc.no/nansat/MOD44W.tgz>`_ from our server and add the path to the directory with this
data to an environment variable named MOD44WPATH (e.g.  ``MOD44WPATH=/Data/sat/auxdata/mod44w``).

Distance to the Nearest coast
------------------------------

To get inforamtion about distace to the nearest coastline within the domain of interest, you can use 
``distance2coast`` method from the ``nanasat.toolbox``. Download NASA's Ocean 
Biology Processing Group `0.1x0.1 degree Distance to the Nearest Coast product 
<https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/GMT_intermediate_coast_distance_01d.zip`_ and add the 
path to the directory with this data to an environment variable named DIST2COAST
(e.g. ``DIST2COAST=/path/to/file``). For more inforamtion about the product wisit `NASA's Ocean 
Biology Processing Group <https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/>`.

Digital Elevation Models (DEMs)
--------------------------------

Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The GMTED2010 datasets are provided by the `U.S. Geological Survey
<https://topotools.cr.usgs.gov/gmted_viewer/>`_. We have prepared a GDAL vrt file that can be used
together with `mapper_topography.py <nansat.mappers.html#module-nansat.mappers.mapper_topography>`__
to open the 30 arcseconds Digital Elevation Model (DEM) with Nansat. To use it, the vrt file must be
downloaded from `<ftp://ftp.nersc.no/nansat/dem/gmted2010_30.vrt>`_ and stored in the same folder as
the tif files of *mean elevation* available at
`<https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Global_tiles_GMTED/300darcsec/mea/>`_.

In case you want to use a different DEM, the procedure is as follows:

#. Download the relevant GDAL readable files to a local folder
#. Generate a vrt file using *gdalbuildvrt*, e.g.:

   .. code-block:: bash

      gdalbuildvrt gmted2010_30.vrt *.tif

#. Update *mapper_topography.py* to accept the new kind of file(s)



Global 30 Arc-Second Elevation (GTOPO30)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have also created a vrt-file for the GTOPO30 dataset. This is available as
`<ftp://ftp.nersc.no/nansat/dem/gtopo30.vrt>`_. The vrt-file should be placed in the same folder as
the .DEM files available at `<https://dds.cr.usgs.gov/srtm/version2_1/SRTM30/>`_.


