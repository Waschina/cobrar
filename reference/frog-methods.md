# FROG: Reproducible fitness and robustness diagnostics for constraint-based models

Run a standardized FROG analysis on a model including (1) Objective
Function Values, (2) Flux Variability Analysis (FVA), (3) Gene Deletion
Fluxes, (4) Reaction Deletion Fluxes.

## Usage

``` r
frog(model, ...)

# S4 method for class 'ModelOrg'
frog(model, outdir = "frog_report")

# S4 method for class 'character'
frog(model, outdir = "frog_report")
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)
  or path to SBML file.

- ...:

  Additional parameters passed on to the specific method instance.

- outdir:

  Path to output directory, where report files will be saved.

## Value

TRUE if FROG report export was successful.

## References

K. Raman et al., “FROG Analysis Ensures the Reproducibility of Genome
Scale Metabolic Models.,” bioRxiv, Sept. 2024. doi:
10.1101/2024.09.24.614797.

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)

output_dir <- tempdir()
frog(mod, outdir = output_dir)
#> [1] TRUE

# list files in output directory
dir(output_dir)
#> [1] "01_objective.tsv"                      
#> [2] "02_fva.tsv"                            
#> [3] "03_gene_deletion.tsv"                  
#> [4] "04_reaction_deletion.tsv"              
#> [5] "bslib-b4cb493c4e74a343d7a14e9c9b491bd8"
#> [6] "downlit"                               
#> [7] "e_coli_core.xml"                       
#> [8] "metadata.json"                         
```
