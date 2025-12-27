# Helper functions for constructing minimal models used in tests

make_minimal_flux_model <- function(import_bound = 10) {
  met_ids <- c("a", "biomass")

  model <- new(
    "ModelOrg",
    mod_id = "toy",
    mod_desc = "toy",
    mod_name = "Toy model"
  )
  
  model <- addCompartment(model, "c", "cytosol")
  
  model <- addMetabolite(model, met_ids)
  
  model <- addReact(model, id = "EX_a",
                    met = "a", Scoef = 1, lb = 0, ub = import_bound)
  model <- addReact(model, id = "a_to_biomass",
                    met = c("a","biomass"), Scoef = c(-1,1), lb = 0, ub = 1000)
  model <- addReact(model, id = "EX_biomass",
                    met = "biomass", Scoef = -1, lb = 0, ub = 1000, obj = 1)
  
  return(model)
}

make_dead_end_model <- function() {
  react_ids <- c("EX_a", "a_to_b", "EX_b", "a_to_c")
  met_ids <- c("a", "b", "c")

  model <- new(
    "ModelOrg",
    mod_id = "dead_end",
    mod_desc = "dead_end",
    mod_name = "Dead-end toy"
  )
  
  model <- addCompartment(model, "c", "cytosol")
  model <- addMetabolite(model, met_ids)
  
  model <- addReact(model, react_ids[1],
                    met = "a", Scoef = -1, lb = -10, ub = 0)
  model <- addReact(model, react_ids[2],
                    met = c("a","b"), Scoef = c(-1,1), lb = 0, ub = 1000)
  model <- addReact(model, react_ids[3],
                    met = "b", Scoef = -1, lb = 0, ub = 10)
  model <- addReact(model, react_ids[4],
                    met = c("a","c"), Scoef = c(-1,1), lb = 0, ub = 5)
  
  return(model)
}
