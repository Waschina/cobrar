
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cobrar

<!-- badges: start -->

[![R-CMD-check](https://github.com/Waschina/cobrar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Waschina/cobrar/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The R-package *cobrar* provides structures and functions to perform
constraint-based metabolic network analysis, e.g., the prediction of
metabolic fluxes using Flux Balance Analysis (FBA). *cobrar* is inspired
by the former CRAN R-package *sybil*.

## Installation

### Prerequisites

*cobrar* requires the two system libraries libSBML and glpk. These
libraries are available from most OS package manager.

#### Ubuntu/Debian/Mint

``` sh
sudo apt install libsbml-dev libglpk-dev
```

#### Centos/Fedora/RHEL

``` sh
sudo yum install glpk-devel
```

#### MacOS

TODO.

### Installation via CRAN

*cobrar* is not yet available on CRAN, but we a working on it.

### Installation via GitHub

You can install the development version of *cobrar* from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Waschina/cobrar")
```

## Usage

The full documentation including illustrative examples is available
[here](https://waschina.github.io/cobrar/).

The following (very brief) example illustrates how to use `cobrar` to
perform flux balance analysis for a core metabolism model of
*Escherichia coli*. In a first simulation, fluxes and growth are
predicted for the default constraints of the model, which represents a
minimal medium with glucose as the sole carbon source and the presence
of oxygen. The second simulation performs the same flux balance analysis
but without oxygen to simulate an anoxic growth environment.

``` r
# First, load the cobrar package and the metabolic model:
library(cobrar)
#> Loading required package: Matrix
#> cobrar uses...
#>  - libSBML (v. 5.18.0)
#>  - glpk (v. 4.65)
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)

# Next, perform simulation 1: Aerobic growth of E. coli
res_aero <- fba(mod)

# Finally, simulation 2: Anaerobic growth of E. coli
mod <- changeBounds(mod, react = "EX_o2_e", lb = 0) # before: -1000
res_anaero <- fba(mod)

# Report simulation results (i.e., growth and acetate production)
cat("[Aerobic growth]\n",
    " Growth rate:        ", res_aero@obj,"\n",
    " Acetate production: ", res_aero@fluxes[mod@react_id == "EX_ac_e"],"\n\n",
    "[Anaerobic growth]\n",
    " Growth rate:        ", res_anaero@obj,"\n",
    " Acetate prodcution: ", res_anaero@fluxes[mod@react_id == "EX_ac_e"],"\n",
    sep = "")
#> [Aerobic growth]
#>  Growth rate:        0.8739215
#>  Acetate production: 0
#> 
#> [Anaerobic growth]
#>  Growth rate:        0.2116629
#>  Acetate prodcution: 8.503585
```

## Key differences to sybil

- cobrar is fully functional from reading SBML files until optimisation
  of linear programs, without the need of additional packages such as
  *sybilSBML* or *glpkAPI*.
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
- currently only glpk is supported as solver. A plugin for IBM’s ILOG
  CPLEX is planned.
- *Multiple objectives*. The SBML standard with its fbc extension allows
  to specify more than one objective (class `ListOfObjectives`).
  However, *cobrar* can only handle one current objective function per
  model, which is defined as an objective coefficient vector in slot
  `obj_coef` of an object of class `modelorg`. Note that when reading
  SBML models, *cobrar* will only use the first objective defined in the
  SBML document.
