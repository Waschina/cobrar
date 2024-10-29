
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cobrar <img src="man/figures/logo.svg" align="right" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/Waschina/cobrar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Waschina/cobrar/actions/workflows/R-CMD-check.yaml)
[![Anaconda-Server
Badge](https://anaconda.org/bioconda/r-cobrar/badges/version.svg?branch=master&kill_cache=1)](https://anaconda.org/bioconda/r-cobrar)
[![Anaconda-Server
Badge](https://anaconda.org/bioconda/r-cobrar/badges/downloads.svg?branch=master&kill_cache=1)](https://anaconda.org/bioconda/r-cobrar)
<!-- badges: end -->

The R package *cobrar* provides structures and functions for
constraint-based metabolic network analysis, e.g. the prediction of
metabolic fluxes using fluxes using Flux Balance Analysis (FBA).
*cobrar* is inspired by the former CRAN R package *sybil*[(1)](#R1).

## Installation

Please note that *cobrar* requires the two system libraries *libSBML*
and *glpk*. The following installation instructions for different
operating systems install the dependencies first, then *cobrar*. If you
already have *libSBML* and *glpk* installed, you can skip to the last
part of the instructions.

*cobrar* is under development and the installation instructions
described will install the latest development version.

#### Ubuntu/Debian/Mint

Install *libSBML* and *glpk*:

``` sh
sudo apt install libsbml-dev libglpk-dev
```

Install cobrar (in *R*):

``` r
# install.packages("remotes")
remotes::install_github("Waschina/cobrar", build_vignettes = TRUE)
```

#### Centos/Fedora/RHEL

Install *libSBML* and *glpk*:

``` sh
sudo yum install libsbml-devel glpk-devel
```

Install cobrar (in *R*):

``` r
# install.packages("remotes")
remotes::install_github("Waschina/cobrar")
```

#### MacOS

*libSBML* and *glpk* can be installed using homebrew:

``` sh
brew install glpk brewsci/bio/libsbml
```

Install cobrar (in *R*):

``` r
# install.packages("remotes")
remotes::install_github("Waschina/cobrar")
```

#### Windows

*cobrar* is currently not available for Windows.

#### Conda

Releases of *cobrar* for linux-64 and osx-64 systems are available via
[bioconda](https://bioconda.github.io/recipes/r-cobrar/README.html#package-r-cobrar).

``` sh
conda install bioconda::r-cobrar
```

## Usage

The full documentation including illustrative examples is available
[here](https://waschina.github.io/cobrar/).

The vignettes of the package can also be accessed:

``` r
library(cobrar)
vignette(package = "cobrar")
vignette(topic = "cobrar")
```

A simple Flux Balance Analysis (FBA) for the core metabolism of
*Escherichia coli*:

``` r
library(cobrar)
#> Loading required package: Matrix
#> cobrar uses...
#>  - libSBML (v. 5.19.0)
#>  - glpk (v. 5.0)

fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)

fba(mod)
#> Algorithm:              FBA 
#> Solver status:          solution is optimal 
#> Optimization status:    optimization process was successful 
#> Objective fct. value:   0.8739215 
#> Secondary objective:    NA
```

## Key differences to sybil

- cobrar is fully functional from reading SBML files until optimisation
  of linear programs, without the need of additional packages such as
  *sybilSBML* or *glpkAPI*.
- cobrar links to libsbml via libsbml’s C++ API, not the C API as the
  *sybilSBML* package.
- The GNU glpk library is a system requirement. There are no
  links/references to IBM’s CPLEX in the *cobrar* package, which were
  one reason why *sybil* was discontinued on CRAN.
- [roxygen2](https://roxygen2.r-lib.org/) is used for documenting
  functions and classes.
- Feature trim: A range of functions in *sybil* are not part of
  *cobrar*.
- In *cobrar*, R’s garbage collector handles memory management,
  including memory associated to C++-objects and the pointers to these.
- Simplifications in class and function architecture
  - No more Class “SysBiolAlg” nor sub-classes. Different algorithms
    have their own function and detailed documentation of their return
    values
- In *sybil*, columns named “annotation” were actually concatenated
  CVTerms
  (<https://synonym.caltech.edu/software/libsbml/5.20.0/cpp-api/class_c_v_term.html>).
  To avoid confusion with other levels of annotation, the columns
  (e.g. in `react_attr` or `met_attr`) are named ‘CVTerms’ in *cobrar*.
- *cobrar* allows to assign SBOTerms to reactions, metabolites, genes.
- Performance (i.e., computation time) improvements in certain
  procedures:
  - identification of dead-end metabolites
  - reading/exporting SBML files,
  - pFBA algorithm (a.k.a MTF ‘Minimization of Total Flux’).
  - FVA; also now allows relaxed constraints on optimal growth
    (e.g. flux variability with 90-100% optimal growth)

## Notes

- *cobrar* exports SBML files level 3 version 2 with `fbc` version 3 and
  `groups` version 1.
- Group assignments are only supported for reactions.
- GLPK is the default solver and is required to build the package. A
  plugin for IBM’s ILOG CPLEX is available [here
  (cobrarCPLEX)](https://github.com/Waschina/cobrarCPLEX).
- *Multiple objectives*. The SBML standard with its `fbc` extension
  allows to specify more than one objective (class `ListOfObjectives`).
  However, *cobrar* can only handle one current objective function per
  model, which is defined as an objective coefficient vector in slot
  `obj_coef` of an object of class `modelorg`. Note that when reading
  SBML models, *cobrar* will only use the objective, which is defined as
  `activeObjective` in the SBML file, or the first objective if no
  active objective is defined.

## References

1.  <a id="R1"></a>Gelius-Dietrich, G., Desouki, A.A., Fritzemeier,
    C.J., Lercher, M.J. sybil – Efficient constraint-based modelling
    in R. BMC Syst Biol 7, 125 (2013).
    <https://doi.org/10.1186/1752-0509-7-125>
