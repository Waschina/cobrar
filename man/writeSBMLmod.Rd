% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{writeSBMLmod}
\alias{writeSBMLmod}
\title{Exports a Metabolic Network in SBML Format}
\usage{
writeSBMLmod(model, file_path = NULL)
}
\arguments{
\item{model}{Model of class \link{ModelOrg}}

\item{file_path}{SBML file name for exporting the model. Default is the
model's ID with ".xml" suffix.}
}
\value{
TRUE if file export was successful.
}
\description{
Export a constraint-based metabolic network model from a S4 object of class
\link{ModelOrg} to a SBML file.
}
\details{
Exported SBML files are of level 3, version 2. FBC-package version 2.

What content from the data.frames \code{react_attr}, \code{met_attr}, and \code{mod_attr} is
exported to SBML files? Currently only the columns named "CVTerms".
}
