% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_modifications.R
\name{changeBounds}
\alias{changeBounds}
\title{Change flux bounds}
\usage{
changeBounds(model, react, lb = NULL, ub = NULL)
}
\arguments{
\item{model}{Model of class \link{ModelOrg}}

\item{react}{A character vector stating the reaction IDs in a model or a
numeric vector providing the reaction indices.}

\item{lb}{A numeric vector giving the new lower flux bounds for reactions
\code{react}. If \code{lb} is of length 1, the same value will be used for
all reactions.}

\item{ub}{A numeric vector giving the new upper flux bounds for reactions
\code{react}. If \code{ub} is of length 1, the same value will be used for
all reactions.}
}
\value{
An updated model of class \link{ModelOrg}
}
\description{
The function changes either upper bounds, lower bounds, or both for specific
reactions.
}
