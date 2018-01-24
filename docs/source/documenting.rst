Documenting Nansat
=====================

Documentation should follow the `conventions
<conventions.html#example-function-with-complete-docstring>`_. 

.. note::

   Documentation for classes should be given after the class definition, not within the
   ``__init__``-method.

To build documentation locally, the best is to create a virtual environment with only the sphinx
tools installed. This is done as follows:

.. code-block:: bash

   mkvirtualenv nansatdocs
   workon nansatdocs
   pip install sphinx sphinx-autobuild sphinx_rtd_theme recommonmark

Then, the following commands should build the documentation:

.. code-block:: bash

   cd docs
   make clean
   sphinx-apidoc -fo source/ ../nansat
   make html

Some documentation remains to be written. This is marked by ``TODO`` in the rst source files. Find
open tasks by:

.. code-block:: bash

   cd docs/source
   grep TODO *
