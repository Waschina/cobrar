# Structure of ModelOrg Class

This class represents a model organism with various attributes related
to central model structures, metabolites, reactions, and genes.

## Slots

- `mod_id`:

  A character vector representing the model identifier.

- `mod_desc`:

  A character vector describing the model.

- `mod_name`:

  A character vector containing the model name.

- `mod_compart`:

  A character vector indicating the model compartment.

- `mod_compart_name`:

  A character vector with the name of the model compartment.

- `mod_attr`:

  A data frame with additional model attributes.

- `mod_notes`:

  A character string that can contain an XML block with additional
  information about the model.

- `S`:

  A sparse numeric matrix of
  [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
  representing the Stoichiometric matrix.

- `obj_coef`:

  A numeric vector containing coefficients for the objective function.

- `obj_dir`:

  Character specifying the objective direction. Either "maximize" or
  "minimize".

- `subSys`:

  A sparse Boolean matrix of
  [lgCMatrix-class](https://rdrr.io/pkg/Matrix/man/lsparseMatrix-class.html)
  defining subsystems.

- `subSys_id`:

  A character vector representing subsystem identifiers.

- `subSys_name`:

  A character vector containing the subsystem names.

- `constraints`:

  An object of class
  [Constraints](https://waschina.github.io/cobrar/reference/Constraints-class.md)
  which specifies constraints in a model in addition to stationarity and
  individual flux bounds.

- `met_id`:

  A character vector representing metabolite identifiers.

- `met_name`:

  A character vector with metabolite names.

- `met_comp`:

  A character vector indicating metabolite compartments.

- `met_attr`:

  A data.frame that enables the storage of additional data for
  metabolites. Only specific columns are exported to SBML files. See
  [writeSBMLmod](https://waschina.github.io/cobrar/reference/writeSBMLmod.md)
  for details.

- `react_id`:

  A character vector representing reaction identifiers.

- `react_name`:

  A character vector with reaction names.

- `react_comp`:

  A character vector indicating reaction compartments.

- `lowbnd`:

  A character vector containing lower bounds for reactions.

- `uppbnd`:

  A character vector containing upper bounds for reactions.

- `react_attr`:

  A data.frame that enables the storage of additional data for
  reactions. Only specific columns are exported to SBML files. See
  [writeSBMLmod](https://waschina.github.io/cobrar/reference/writeSBMLmod.md)
  for details.

- `gprRules`:

  A character vector with Gene-Protein-Reaction association rules (with
  gene product indices corresponding to the order in slot 'genes').

- `genes`:

  A list of character vectors. Each vector contains the IDs of gene
  products associated to the respective reaction.

- `allGenes`:

  A character vector with all gene identifiers.

- `genes_attr`:

  A data.frame that enables the storage of additional data (e.g., name
  and CVTerms) for genes/gene products. Only specific columns are
  exported to SBML files. See
  [writeSBMLmod](https://waschina.github.io/cobrar/reference/writeSBMLmod.md)
  for details.

- `metadata`:

  A list, which can be used to include any arbitrary content describing
  the metabolic model. Important: No content of 'metadata' is written to
  SBML exports of the model.
