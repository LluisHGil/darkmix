Structure of the repository
===========================

The `darkmix repository`_ contains the code of the project, the sample data used in the :ref:`publication`, as well as the output files obtained in that analysis.


Code
-----

The folder ``Code`` contains the following files:

* ``darkmix.R``: contains the *R* functions used in the library. It must be loaded every time the library is used.
* ``kernel_lambda.cpp``: C++ code for the kernel density calculations. Not meant to be used directly by the user.
* ``kernel_absolute_residuals.cpp``: C++ code for the residuals kernel density estimation. Not meant to be used directly by the user.
* ``kernel_fourier_residuals.cpp``: fast C++ code for the residuals kernel density estimation. Not meant to be used directly by the user. Requires the installation of the `FFTW library <http://www.fftw.org/>`_ (see :doc:`installation`).

* ``darkmix_steps.Rmd``: notbook containing the execution steps as used in the publication. Uses the data stored in folder ``Data``. You can see this example of execution in :doc:`example1` (**CHECK!!**).


Data
-----

The folder ``Data`` contains the following files:

* ``datacat.txt``: data set containing the coordinates of 2081 dark matter particles from the `Bolshoi simulation`_ . Used in the ``darkmix_steps.Rmd`` example code.
* ``halocat.txt``: catalog of halos found by the `BDM algorithm`_ in our region of interest.
* ``bdm_halos.txt``: selection of 10 halos from the ``halocat.txt`` catalog coincident with the ones found in the ``datacat.txt`` sample following ``darkmix_steps.Rmd``.

Output
-------

The folder ``Output`` contains the following files:

* ``parameters.txt``: best fit paramters of our 11-component mixture model. The last row corresponds to the background component and column *N* is the number of estimated particles per component.
* ``membership.txt``: probability for each of the 2081 particles in ``datacat.txt`` of belonging to each of the 11 components. Column *class* is the component to which the particle has been assigned in our example.


.. _`darkmix repository`: https://github.com/LluisHGil/darkmix
.. _`Bolshoi simulation`: https://www.cosmosim.org/cms/documentation/projects/multidark-bolshoi-project/
.. _`BDM algorithm`: https://www.cosmosim.org/cms/documentation/database-structure/tables/bdm/
