% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_modifications.R
\name{rmCompartment}
\alias{rmCompartment}
\title{Remove compartments from a model}
\usage{
rmCompartment(model, comp)
}
\arguments{
\item{model}{Model of class \link{ModelOrg}}

\item{comp}{A character vector stating the compartment IDs in a model or a
numeric vector providing the compartment indices.}
}
\value{
An updated model of class \link{ModelOrg}
}
\description{
This function removes specified compartments from a model.
}
\note{
If at least one of the provided compartments still has metabolites associated
with it, the function stops with an error message.
}
\seealso{
Other Model manipulation tools: 
\code{\link{addCompartment}()},
\code{\link{addConstraint}()},
\code{\link{addGene}()},
\code{\link{addMetabolite}()},
\code{\link{addReact}()},
\code{\link{addSubsystem}()},
\code{\link{changeBounds}()},
\code{\link{rmConstraint}()},
\code{\link{rmGene}()},
\code{\link{rmMetabolite}()},
\code{\link{rmReact}()},
\code{\link{rmSubsystem}()}
}
\concept{Model manipulation tools}