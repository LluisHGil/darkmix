# darkmix
Mixture Models for detection and characterization of Dark Matter halos

This repository contains the codes included in the `darkmix` project (Hurtado-Gil et al., in prep.). This library uses Finite Mixture Models to estimate the distribution of dark matter particle samples. Dark matter halos present a variety of shapes and sizes that needs to be correctly classified and characterized. Traditional non parametric methods do not provide enough information or misclassify individuals in merging situations. Using the Nelder-Mead algortihm and MCMC routines, we can correctly estimate a Finite Mixture Model, and obtain quantities of phyisical interest. The library assumes the Einasto profile for dark matter halos.

In folder Code one can found the following files:

- `darkmix.R`: contains the `R` functions used in the library. It must be loaded everytime the library is used.
- `darkmix_steps.R`: execution steps as used in the publication. Uses the data stored in folder `Data`.
- `kernel_lambda.cpp`: C++ code for the kernel density calculations. Not meant to be used directly by the user.
- `kernel_absolute_residuals.cpp`: C++ code for the residuals kernel density estimation. Not meant to be used directly by the user.
- `kernel_fourier_residuals.cpp`: fast C++ code for the residuals kernel density estimation. Not meant to be used directly by the user. Requires the installation of the [FFTW library](http://www.fftw.org/) (see below). 

In folder Data, one can found the following files:

- `datacat.txt`: data set containing the coordinates of 2081 dark matter particles from [Bolshoi](https://www.cosmosim.org/cms/documentation/projects/multidark-bolshoi-project/) simulation (Spanish MultiDark Consolider project). Used in the `darkmix_steps.R` example code.
- `halocat.txt`: catalog of halos found by the [BDM algorithm](https://www.cosmosim.org/cms/documentation/database-structure/tables/bdm/) in our region of interest.
- `bdm_halos.txt`: selection of 10 halos from the `halocat.txt` catalog coincident with our found halos in `datacat.txt` sample.

On folder Output one can found the following files:

- `parameters.txt`: best fit paramters of out 11 components mixtire model. The last row correspond to the background component and column *N* is the number of estimated particles per component.
- `membership.txt`: probability of each of the 2081 particles of belonging to each of the 11 components. Column *class* is the component to which the particle has been assigned in our example.

Installation:

This library is not prepared to be installed but only as a list of functions that can be loaded and used locally. To do so, the C++ files must be compiled and all executions have to be made from the parent folder. After downloading the repository, access to folder

`../darkmix/Code`

and compile the C++ code with

```
g++ -o kernel_lambda kernel_lambda.cpp
g++ -o kernel_absolute_residuals kernel_absolute_residuals.cpp
```

If the used decides to use the Fast Fourier transformations for the residuals kernel density estimation (recommended), the [FFTW](http://www.fftw.org/) library should be installed first. Once the installations is completed, acces the `Code` folder again and compile the C++ function with

```
g++ -o kernel_fourier_residuals kernel_fourier_residuals.cpp -lfftw3 -lm
```

Another necessary software is the [Voro++](http://math.lbl.gov/voro++/) library, used to calculate the volume of 3d voronoi tessellations. Download and install before using the `darkmix.R` functions.

Once these steps are completed, one can open `R` and install the following necessary `R` packages (as explained in `darkmix_steps.R`):

```
install.packages(c("spatstat", "spatstat.utils", "LaplacesDemon", "expint", "plotrix", "misc3d", "moments"))
```

Now we only have to set the working directory and load the `darkmix.R` functions:

```	
setwd("~/PATH/darkmix/")
source("Code/darkmix.R")
```
