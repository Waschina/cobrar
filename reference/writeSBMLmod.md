# Exports a Metabolic Network in SBML Format

Export a constraint-based metabolic network model from a S4 object of
class [ModelOrg](ModelOrg-class.md) to a SBML file.

## Usage

``` r
writeSBMLmod(model, file_path = NULL)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- file_path:

  SBML file name for exporting the model. Default is the model's ID with
  ".xml" suffix.

## Value

TRUE if file export was successful.

## Details

Exported SBML files are of level 3, version 2. FBC-package version 2.

What content from the data.frames `react_attr`, `met_attr`, and
`mod_attr` is exported to SBML files? Currently only the columns named
"CVTerms".
