# S4-Object for metabolic model

setClass("Modelorg",

         slots = c(
           # central model structures
           mom_id = "character",
           mod_desc = "character",
           mod_name = "character",
           mod_compart = "character",
           mod_attr = "data.frame",
           S = "Matrix",
           obj_coef = "numeric",
           subSys = "Matrix",

           # metabolites,
           met_id = "character",
           met_name = "character",
           met_comp = "character",
           met_attr = "character",

           # reactions
           react_id = "character",
           react_name = "character",
           lowbnd = "character",
           uppbnd = "character",

           # genes
           gprRules = "character",
           genes = "list",
           allGenes = "character"


         ))
