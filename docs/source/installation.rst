Installation
============

Pre-requisites
---------------

In order to use this library, it is necessary to install these additional packages:

Voro++
```````

`Voro++`_ is a library for carrying out three-dimensional computations of the Voronoi tessellation, needed by several *darkmix* functions.
To install it, follow the instructions in their website.


FFTW (recommended)
```````````````````

`FFTW`_ is a standard library to compute discrete Fourier transforms in a fast way.
It is necessary to compile the ``kernel_fourier_residuals`` programme.
If FFTW is not installed, the code will estimate the residuals instead using the ``kernel_absolute_residuals``, which uses a much slower algorithm.

The FFTW library can be downloaded and installed from its website, and it is also available as a package in the standard repositories of most operating systems.
If installed using the later method (recommended), the *development* version of the library must be installed (e.g., the package ``libfftw3-dev`` in Ubuntu).

Additional R packages
``````````````````````

A few additional *R* packages are needed to run the *darkmix* code.
These can be installed from within *R* using:

.. code-block:: r

  packages <- c("spatstat", "spatstat.utils", "LaplacesDemon", "expint", "plotrix", "misc3d", "moments")
  install.packages(packages)




Compiling the code
--------------------

The *darkmix* consists of a series of functions contained in ``darkmix.R`` that can be directly loaded and used locally.
These functions rely on a series of *C++* programmes that should be previously compiled.

After downloading and unpacking the code, access the ``Code`` folder:

.. code-block:: bash

   cd PATH/darkmix/Code

A ``Makefile`` has been created for convenience.
Modify it if any change is needed for your local configuration, and compile the code using:

.. code-block:: bash

  make all

This should create the three needed executables: ``kernel_lambda``, ``kernel_absolute_residuals`` and ``kernel_fourier_residuals``.


R setup
---------

Once the code is compiled, to use the *darkmix* functions from within *R*, we need to set the working directory and load the ``darkmix.R`` functions:

.. code-block:: r

   setwd("~/PATH/darkmix/")
   source("Code/darkmix.R")





.. _`Voro++`: http://math.lbl.gov/voro++/
.. _`FFTW`: http://www.fftw.org/
