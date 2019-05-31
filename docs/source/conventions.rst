Nansat conventions
==================

Git branching and merging
-------------------------

We adopt the following system for branching and merging:

1. **master** branch: numbered releases of the code. Never edited. Merged from *develop* and *hot fix* branches (see notes on workflow below). Long living.
2. **develop** branch: rather stable version of code under development. Never edited. Merged from topic specific issue branches. Long living.
3. **issue<NNN>_<short-heading>** - issue specific branches (NNN = issue number). Main working area. Short living. Branched from, and merged back into develop.
4. **hotfix<NNN>_<short-heading>** - branches that are specific to a hotfix issue. Hotfixes are bugfixes on master that can not wait until next release.

.. note::

   1. Never edit code in the master or develop branch. Always make a new branch for your edits.
   2. A new branch should be very specific to only one problem. It should be short living.
   3. Commit often.
   4. Branch often.
   5. Branch only from master or from develop.
   6. Create pull requests for your branches and **always** assign a reviewer to merge, delete the branch, and close the issue (this is easy in github)

How to report and handle new issues
-----------------------------------

If you discover a bug in Nansat or you would like to suggest improvements or new features to Nansat, the following procedure should be followed:

1. Check that no one else has reported the same issue at https://github.com/nansencenter/nansat/issues
2. If not, add a new issue at https://github.com/nansencenter/nansat/issues
3. If Nansat respository is not accessible for writing (non team members), fork Nansat and clone your own fork locally
4. Create an issue branch on your local system named **issue<NNN>_<short-heading>** where NNN is the issue number from GitHub. This will be the main (short living) working area. The issue branch should originate from develop.
5. Add tests to reproduce the bug or test the new functionality
6. Write the necessary code
7. Push the new issue branch to your own fork at GitHub
8. Create a pull request of the issue branch on your fork at GitHub. A member of our team will the review the code, merge, delete the branch, and close the issue (this is easy in github)

How to report and handle a hotfix
---------------------------------

If a bug is relatively quick and easy to handle, we call it a hotfix. A hotfix is branched from master (by team members), and the following workflow applies:

1. Branch from master into the hotfix specific branch (**hotfix<NNN>_<short-heading>**)

  a) Update the tests
  b) Fix the bug
  c) Increment micro version in ``setup.py``
  d) Commit to the hotfix branch

2. If needed rebase the hotfix branch on top of master:

  a) Checkout the hotfix branch
  b) Use ``git rebase master``
  c) Fix conflicts if any
  d) Test the code
  e) Push your hotfix branch to GitHub. Note: You have to use 'git push -f' in order to rewrite the history on github. Rebase will not cause problems if only one person works with the hotfix branch.

3. Go to `<https://github.com/nansencenter/nansat>`_ and add a pull request for the newly pushed
   hotfix branch and assign a reviewer

4. Let the reviewer do the following:

  a) Wait for tests on Travis CI to pass
  b) Check the code
  c) Request changes or merge the pull request into master using 'Rebase and Merge' button in the online tool.
  d) Delete the branch
  e) Merge master into develop (``git checkout develop``; ``git merge master``; test the code; push to GitHUB; check Travis CI status)
  f) Close the related issue

.. note::

    If you work on a project using Nansat, this project should use the master version of Nansat. Different
    situations may then occur:

    1. You discover a **bug** that is easy to fix. This should be solved in a **hotfix** branch originating from master.
    2. You discover a **bug** that is difficult to fix. This should be solved in an **issue** branch originating from develop.
    3. You want to add **functionality** to Nansat. This should be solved in an **issue** branch originating from develop. When this is completed and merged back into develop, we may have a new release.
    4. You want to add a mapper. This can be done by adding a package ``nansat_mappers`` to your
       project. When the mapper is completed, create an issue and an issue-branch in Nansat, and
       finally submit a pull request with suggestion of a reviewer. You can still use the mapper via
       your ``nansat_mappers`` package until the new mapper is implemented in a release version of
       Nansat.

General conventions
-------------------

* Nansat coding style generally follows `PEP-8 (General style guide)
  <http://www.python.org/dev/peps/pep-0008/>`_ and `PEP-257 (Docstrings)
  <http://www.python.org/dev/peps/pep-0257/>`_
* Max line length is set to 100 chars
* Every unit of code must be properly tested (see unit-test) and documented
* All class/function/method/variable names have to be explicit and should contain no more than 3 words
* Single quotes should be used consistently instead of double quotes (except for cases where quotes
  are required, and for docstrings)
* GNU v3 licence should be inserted in all files. Mappers should have a standard header like this:

.. code-block:: python

   # Name:         mapper_asar.py
   # Purpose:      Mapper for Envisat/ASAR data
   # Authors:      Asuka Yamakava, Anton Korosov
   # Licence:      This file is part of NANSAT. You can redistribute it or modify
   #               under the terms of GNU General Public License, v.3
   #               http://www.gnu.org/licenses/gpl-3.0.html
   #
   # Additional mapper/format specific links and information

* Docstrings should follow the `Numpy style
  <https://numpydoc.readthedocs.io/en/latest/>`_
  , with an extensive example file here `example.py
  <https://numpydoc.readthedocs.io/en/latest/example.html#example>`_
* Valid headers for functions and classes are, in the following order, 'Parameters', 'Attributes' (only class),
  'Returns' (only function), 'Yields', 'Other parameters', 'Raises', 'Warns', 'Warnings',
  'See also', 'Notes', 'References' and 'Examples'.

Example function with complete Docstring
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   def some_function(start = 0, stop, step = 1):
       """ Return evenly spaced values within a given interval.

       Values are generated within the half-open interval ''[start, stop)''
       (in other words, the interval including 'start' but excluding 'stop').

       For integer arguments the function is equivalent to the Python built-in
       'range '_ function, but returns a ndarray rather than a list.

       Parameters
       ----------
       start : number, optional
           Start of interval.  The interval includes this value.  The default start value is 0.
       stop : number
           End of interval.  The interval does not include this value.
       step : number, optional
           Spacing between values.  For any output 'out', this is the distance between two adjacent values, ''out[i+1] - out[i]''. The default step size is 1. If 'step' is specified, 'start' must also be given.
       dtype : dtype
           The type of the output array. If 'dtype' is not given, infer the data type from the other input arguments.

       Returns
       -------
       out : ndarray
           Array of evenly spaced values.

           For floating point arguments, the length of the result is ''ceil((stop - start)/step)''. Because of floating point overflow, this rule may result in the last element of 'out' being greater than 'stop'.

       See Also
       --------
       linspace : Evenly spaced numbers with careful handling of endpoints
       ogrid: Arrays of evenly spaced numbers in N-dimensions
       mgrid: Grid-shaped arrays of evenly spaced numbers in N-dimensions

       Examples
       --------
       >>> np.arange(3)
       array([0, 1, 2])

       """

Naming conventions
-------------------------

* when a variable points to the GDALDataset, GDALDriver, etc. its name must always contain word "dataset", "driver", etc. representatively (raw_dataset, src_dataset, example_driver)
* when a variable points to a string with name it should contain 'name' (band_name)
* when longitude and latitude are input to (or output from) a function, they should be given in this order: (lon, lat). These variables should always be named 'lon' and 'lat' (i.e. never 'long').
* source and destination are prefixed as 'src' and 'dst' (src_dataset,  dst_raster_xsize)
* band numbers should be called ‘band_number’
* GDAL bands should be called 'band' or, e.g., ‘dst_band’ when prefixed (GDAL is actually in-consistent here: gdal.Dataset.!GetRasterBand returns a 'Band'-object; hence 'Band' is the name of the class and the Python datatype)
* We use ‘filename’ (as in Python standard library)

Style checking
--------------

In your IDE/editor, it is highly recommended to activate/install a plugin for/script a save hook for
doing automatic style checks and/or corrections, eg autopep8, pylint, pyflakes.

Tests
------------

In general:

* Every function must be accompanied with a test suite
* Tests should be both positive (testing that the function work as intended with valid data) and negative (testing that the function behaves as expected with invalid data e.g. that correct exceptions are thrown)
* If a function has optional arguments, separate tests for all options should be created

Testing core Nansat functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Tests for Nansat, Domain, etc should be added to nansat/tests/test_<module_name>.py file;
* These tests should be added as functions of classes inheriting from unittest.TestCase (e.g. DomainTest);
* Tests sharing similar set-up may inherit from the same class which has a setUp function;
* The core tests are run at `Travis CI <https://travis-ci.org/nansencenter/nansat>`_ (continuous integration) which integrates with `Coveralls <https://coveralls.io/r/nansencenter/nansat>`_ for providing test coverage

Integration testing
^^^^^^^^^^^^^^^^^^^^

Products read by Nansat mappers are tested in modules within the nansat_integration_tests folder in
the repository root. These tests should have access to all the kinds of data read by nansat. Since
this is a very large amount of data, and since we cannot share every data product openly, these
tests are not presently executed at Travis CI. Every developer should add new end-to-end tests and
execute them when new mappers or workflows are added. Unavailable test data will lead to fewer tests
being executed, i.e. they won't fail because of missing data. If possible, datasets used in new
tests should be made available to the Nansen Center such that we can run the full test suite.


Testing mappers
^^^^^^^^^^^^^^^

General tests checking that the mappers don't violate the functionality of nansat and checks that
some specific metadata is added, are collected in the nansat_integration_tests.mapper_tests module.

Also, we aim to create proper unit tests that use mock object for all the mappers. This will help to
significantly increase the test coverage.

Testing specific data products or workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In typical scientific workflows, a data product is opened with Nansat and some operations are
performed, e.g., adding new derived bands and exporting the results to a netcdf, or creating figures.
To make sure that new versions of nansat do not harm these workflows with bugs or sudden interface
changes, we collect tests for typical workflows in separate modules within the
nansat_integration_tests package, e.g. test_sar, test_radarsat2, etc. We encourage users and
developers to add such tests to avoid such potential problems

Doctests
^^^^^^^^^^^^

TODO: add information about how to use doctests
