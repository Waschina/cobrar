# Helper functions for constructing minimal models used in tests

make_minimal_flux_model <- function(import_bound = 10) {
  react_ids <- c("EX_a", "BIOMASS")
  met_ids <- c("a", "biomass")

  S <- Matrix::Matrix(
#    c( 1, 0,
#      -1, 1
#    ),
    c( 1, -1,
       0,  1
    ),
    nrow = length(met_ids),
    ncol = length(react_ids),
    byrow = FALSE,
    sparse = TRUE,
    dimnames = list(met_ids, react_ids)
  )
  S <- as(S, "generalMatrix")

  constraints <- new(
    "Constraints",
    coeff = Matrix::Matrix(0, nrow = 0, ncol = length(react_ids), sparse = TRUE),
    lb = numeric(),
    ub = numeric(),
    rtype = character()
  )

  new(
    "ModelOrg",
    mod_id = "toy",
    mod_desc = "toy",
    mod_name = "Toy model",
    mod_compart = "c",
    mod_compart_name = "cytosol",
    mod_attr = data.frame(
      CVTerms = "",
      SBOTerm = "SBO:0000624",
      stringsAsFactors = FALSE
    ),
    mod_notes = NA_character_,
    S = S,
    obj_coef = c(0, 1),
    obj_dir = "maximize",
    subSys = Matrix::Matrix(
      FALSE,
      nrow = 0,
      ncol = length(react_ids),
      sparse = TRUE
    ),
    subSys_id = character(),
    subSys_name = character(),
    constraints = constraints,
    met_id = met_ids,
    met_name = met_ids,
    met_comp = rep("c", length(met_ids)),
    met_attr = data.frame(
      chemicalFormula = rep("", length(met_ids)),
      charge = rep(NA_real_, length(met_ids)),
      CVTerms = rep("", length(met_ids)),
      SBOTerm = rep("", length(met_ids)),
      stringsAsFactors = FALSE
    ),
    react_id = react_ids,
    react_name = react_ids,
    react_comp = rep("c", length(react_ids)),
    lowbnd = c(0, 0),
    uppbnd = c(import_bound, 1000),
    react_attr = data.frame(
      CVTerms = rep("", length(react_ids)),
      SBOTerm = rep("", length(react_ids)),
      stringsAsFactors = FALSE
    ),
    gprRules = rep("", length(react_ids)),
    genes = replicate(length(react_ids), character(), simplify = FALSE),
    allGenes = character(),
    genes_attr = data.frame(
      name = character(),
      CVTerms = character(),
      SBOTerm = character(),
      stringsAsFactors = FALSE
    ),
    metadata = list()
  )
}

make_dead_end_model <- function() {
  react_ids <- c("EX_a", "a_to_b", "b_to_a_rev", "a_to_c")
  met_ids <- c("a", "b", "c")

  S <- Matrix::Matrix(
    c(
      1, -1,  1, -1,
      0,  1, -1,  0,
      0,  0,  0,  1
    ),
    nrow = length(met_ids),
    ncol = length(react_ids),
    byrow = TRUE,
    sparse = TRUE,
    dimnames = list(met_ids, react_ids)
  )
  S <- as(S, "generalMatrix")

  constraints <- new(
    "Constraints",
    coeff = Matrix::Matrix(0, nrow = 0, ncol = length(react_ids), sparse = TRUE),
    lb = numeric(),
    ub = numeric(),
    rtype = character()
  )

  new(
    "ModelOrg",
    mod_id = "dead_end",
    mod_desc = "dead_end",
    mod_name = "Dead-end toy",
    mod_compart = "c",
    mod_compart_name = "cytosol",
    mod_attr = data.frame(
      CVTerms = "",
      SBOTerm = "SBO:0000624",
      stringsAsFactors = FALSE
    ),
    mod_notes = NA_character_,
    S = S,
    obj_coef = c(0, 0, 0, 0),
    obj_dir = "maximize",
    subSys = Matrix::Matrix(
      FALSE,
      nrow = 0,
      ncol = length(react_ids),
      sparse = TRUE
    ),
    subSys_id = character(),
    subSys_name = character(),
    constraints = constraints,
    met_id = met_ids,
    met_name = met_ids,
    met_comp = rep("c", length(met_ids)),
    met_attr = data.frame(
      chemicalFormula = rep("", length(met_ids)),
      charge = rep(NA_real_, length(met_ids)),
      CVTerms = rep("", length(met_ids)),
      SBOTerm = rep("", length(met_ids)),
      stringsAsFactors = FALSE
    ),
    react_id = react_ids,
    react_name = react_ids,
    react_comp = rep("c", length(react_ids)),
    lowbnd = c( 0,    0, -10, 0),
    uppbnd = c(10, 1000,   0, 5),
    react_attr = data.frame(
      CVTerms = rep("", length(react_ids)),
      SBOTerm = rep("", length(react_ids)),
      stringsAsFactors = FALSE
    ),
    gprRules = rep("", length(react_ids)),
    genes = replicate(length(react_ids), character(), simplify = FALSE),
    allGenes = character(),
    genes_attr = data.frame(
      name = character(),
      CVTerms = character(),
      SBOTerm = character(),
      stringsAsFactors = FALSE
    ),
    metadata = list()
  )
}
