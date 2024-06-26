% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LPproblemClass.R, R/LPproblem_glpkClass.R
\docType{methods}
\name{setObjDirection}
\alias{setObjDirection}
\alias{setObjDirection,LPproblem_glpk-method}
\alias{setObjDirection,LPproblem_glpk}
\title{Set objective direction}
\usage{
setObjDirection(lp, ...)

\S4method{setObjDirection}{LPproblem_glpk}(lp, lpdir)
}
\arguments{
\item{lp}{Object of class \link{LPproblem}}

\item{...}{Additional parameters passed on to the specific method instance.}

\item{lpdir}{Objective direction ("max" or "min")}
}
\description{
Set the objective function direction in an \link{LPproblem}.
}
