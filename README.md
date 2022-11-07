[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6854940.svg)](https://doi.org/10.5281/zenodo.6854940)

# darkmix
Mixture Models for detection and characterization of Dark Matter halos


This repository contains the code included in the `darkmix` project.
This library uses Finite Mixture Models to estimate the distribution of dark matter particle samples.
Dark matter halos present a variety of shapes and sizes that needs to be correctly classified and characterized.
Traditional non parametric methods do not provide enough information or misclassify individuals in merging situations.
Using the Nelder-Mead algorithm and MCMC routines, we can correctly estimate a Finite Mixture Model, and obtain quantities of physical interest.
The library assumes the Einasto profile for dark matter halos.


## Documentation

You can find the online documentation for `darkmix` at: <https://darkmix.readthedocs.io>


## Publication

The motivation for this project, as well as a detailed description of the algorithm and the code can be found in:

- Ll. Hurtado-Gil, M. A. Kuhn, P. Arnalte-Mur, E. D. Feigelson, V. J. Martínez (2022):
  *DarkMix: Mixture Models for Detection and Characterization of Dark Matter Halos*,
  The Astrophysical Journal, 939:34; [DOI:10.3847/1538-4357/ac88d4](https://doi.org/10.3847/1538-4357/ac88d4);
  [arXiv:2208.04194](https://arxiv.org/abs/2208.04194)


## License

This code is open source, according to the terms of the [MIT License](https://choosealicense.com/licenses/mit/) (see `LICENSE.txt`)
