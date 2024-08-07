% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_modifications.R
\name{addReact}
\alias{addReact}
\title{Add or modify a reaction}
\usage{
addReact(
  model,
  id,
  met,
  Scoef,
  reversible = FALSE,
  lb = 0,
  ub = COBRAR_SETTINGS("MAXIMUM"),
  obj = 0,
  subsystem = NA,
  subsystemName = NA,
  gprAssoc = NA,
  reactName = NA,
  metName = NA,
  metComp = NA,
  metCharge = NA,
  metChemicalFormula = NA,
  CVTerms = NA,
  SBOTerm = "guess"
)
}
\arguments{
\item{model}{Model of class \link{ModelOrg}}

\item{id}{Character for the reaction ID}

\item{met}{Character vector providing the IDs of metabolites that participate
in the reaction}

\item{Scoef}{Numeric vector (same length as \code{met}) of stoichiometric
coefficients for the metabolites in \code{met}. The value in \code{Scoef[i]} is the
stoichiometric coefficient of the metabolite in \code{met[i]}.}

\item{reversible}{This option has now effect and is only here for legacy
reasons. Whether a reaction is reversible or not is inferred by cobrar based
on the lower and upper bounds.}

\item{lb, ub}{Single numeric values that define the lower and upper flux
limits, respectively.}

\item{obj}{Single numeric value for the coefficient of the reaction in the
objective function.}

\item{subsystem}{A vector of character strings containing the sub system IDs
to which the reaction belongs.}

\item{subsystemName}{A character vector (same length as \code{subsystem}) for the
names of the subsystems. If the subsystem is already part of the model and
you do not want to change its name, just use NA the corresponding entry.}

\item{gprAssoc}{A single character string giving the Gene-Product-Reaction
(GPR) association for the reaction. If NA, no GRP association is created.}

\item{reactName}{A single character string giving the name for the reaction.
If NA, the value of argument \code{id} is used.}

\item{metName}{A vector of character strings of the same length as \code{met}
containing the metabolites names for the metabolites given in argument \code{met}.}

\item{metComp}{A vector of character strings of the same length as \code{met}
specifying the compartment IDs for the metabolites given in argument \code{met}.}

\item{metCharge}{A numeric vector of the same length as \code{met} defining the
charges for the metabolites given in argument \code{met}.}

\item{metChemicalFormula}{A character vector of the same length as \code{met}
defining the chemical formulas for the metabolites given in argument \code{met}.}

\item{CVTerms}{Cross-references to other resources.}

\item{SBOTerm}{A term ID from the Systems Biology Ontology. If 'guess', the
SBO will be inferred from the reaction details.}
}
\description{
The function can be used to add or modify a reaction in an existing model.
}
\details{
If you want to use the function to update data of a pre-existing reaction but
not its stoichiometry, use NA for the parameters 'met' and 'Scoeff'.
If the reaction is already part of the model, any reaction value (e.g., lb,
ub, reactName), that is set to NA has the effect that the old value will be
used.
If the reaction is already part of the model, and values for the parameter
\code{subsystem} are provided, all previous set subsystem associations of the
reaction will be removed.
If metabolites or subsystems are not part of the model yet, they will be
added.
}
\examples{
# This example adds the 4-aminobutyrate degradation pathway to the E. coli
# core metabolic model
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)

fba(mod)

# 4abut transport: 4abut_e + h_e <=> 4abut_c + h_c
mod <- addReact(mod, id = "ABUTt", Scoef = c(-1,-1,1,1),
                met = c("4abut_e","h_e","4abut_c","h_c"), reversible = TRUE,
                lb = -1000, ub = 1000,
                reactName = "4-aminobutyrate transport in via proton symport",
                metName = c("4-aminobutyrate",NA, "4-aminobutyrate",NA),
                metComp = c("e","e","c","c"), metCharge = c(0,NA,0,NA),
                metChemicalFormula = c("C4H9NO2",NA,"C4H9NO2",NA))

# exchange reaction for 4abut (with 1.5 mmol/gDW/hr availability)
mod <- addReact(mod, id = "EX_4abut_e", Scoef = c(-1), met = "4abut_e",
                lb = -1.5, ub = 1000, reactName = "4-aminobutyrate exchange")

# 4abut amninotransferase (EC 2.6.1.19)
mod <- addReact(mod, id = "ABTA", Scoef = c(-1,-1,1,1),
                met = c("4abut_c","akg_c","glu__L_c","sucsal_c"),
                lb = 0,
                reactName = "4-aminobutyrate transaminase",
                metName = c(NA,NA,NA,"Succinic semialdehyde"),
                metComp = c(NA,NA,NA,"c"), metCharge = c(NA,NA,NA,-1),
                metChemicalFormula = c(NA,NA,NA,"C4H5O3"),
                subsystem = "GABAdegr", subsystemName = "4-aminobutyrate degradation",
                CVTerms = "bqbiol_is;http://identifiers.org/ec-code/2.6.1.19",
                gprAssoc = "b2662 | b1302")

# Succinate-semialdehyde dehydrogenase (NAD) (EC 1.2.1.24)
mod <- addReact(mod, id = "SSALx", Scoef = c(-1,-1,-1,2,1,1),
                met = c("h2o_c","nad_c","sucsal_c","h_c","nadh_c","succ_c"),
                lb = 0,
                reactName = "Succinate-semialdehyde dehydrogenase (NAD)",
                subsystem = "GABAdegr",
                CVTerms = "bqbiol_is;http://identifiers.org/ec-code/1.2.1.24",
                gprAssoc = "b1525")

printReaction(mod, "SSALx")

fba(mod)

}
