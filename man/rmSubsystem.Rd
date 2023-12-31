% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_modifications.R
\name{rmSubsystem}
\alias{rmSubsystem}
\title{Remove subsystems from a model}
\usage{
rmSubsystem(model, subsystem)
}
\arguments{
\item{model}{Model of class \link{modelorg}}

\item{subsystem}{A character vector stating the subsystem IDs in a model or a
numeric vector providing the subsystem indices.}
}
\value{
An updated model  of class \link{modelorg}
}
\description{
This function removes specified subsystems from a model.
}
