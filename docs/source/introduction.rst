Introduction
============
*darkmix*: Mixture Models for detection and characterization of Dark Matter halos


This `darkmix repository`_ contains the code included in the *darkmix* project.
This library uses Finite Mixture Models to estimate the distribution of dark matter particle samples.
Dark matter halos present a variety of shapes and sizes that needs to be correctly classified and characterized.
Traditional non parametric methods do not provide enough information or misclassify individuals in merging situations.
Using the Nelder-Mead algorithm and MCMC routines, we can correctly estimate a Finite Mixture Model, and obtain quantities of physical interest.
The library assumes the Einasto profile for dark matter halos.



.. _publication:

Publication
-----------

The motivation for this project, as well as a detailed description of the algorithm and the code can be found in:

- Ll. Hurtado-Gil, M. A. Kuhn, P. Arnalte-Mur, E. D. Feigelson, V. J. Martínez:
  *DarkMix: Mixture Models for Detection and Characterization of Dark Matter Halos* (submitted to ApJ)


Authors
-------

The *darkmix* code was written by Lluís Hurtado-Gil (eDreams ODIGEO / Universitat de València) and Michael A. Kuhn (Caltech), with contributions from Pablo Arnalte-Mur (Universitat de València).


Support and Contributions
-------------------------

This code is maintained by Lluís Hurtado-Gil (Lluis.Hurtado@uv.es) and Pablo Arnalte-Mur (Pablo.Arnalte@uv.es).
To get help when using the code, and for any bug, feature requests or general ideas, please `raise an issue`_ at the GitHub repository.

If you want to contribute code to the project, please simply fork the project on GitHub and then raise a pull request.
Pull requests will be reviewed and we will comment, suggest changes or alternatives and/or accept them via GitHub.

License
-------

This code is open source, according to the terms of the `MIT License`_ (see ``LICENSE.txt``)



.. _`darkmix repository`: https://github.com/LluisHGil/darkmix

.. _`MIT License`: https://choosealicense.com/licenses/mit/

.. _`raise an issue`: https://github.com/LluisHGil/darkmix/issues
