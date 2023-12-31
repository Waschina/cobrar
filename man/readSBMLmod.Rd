% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbml_io.R
\name{readSBMLmod}
\alias{readSBMLmod}
\title{Reads an SBML file and constructs an object of class 'modelorg'}
\usage{
readSBMLmod(file_path)
}
\arguments{
\item{file_path}{Path to SBML file.}
}
\value{
A \link{modelorg-class} object.
}
\description{
Reads an SBML file and constructs an object of class 'modelorg'
}
\examples{
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)
mod

}
