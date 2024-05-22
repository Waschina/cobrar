#' Structure of ModelComm Class
#'
#' This class represents a community of model organisms.
#'
#' @slot community A data.frame
#'
#' @aliases ModelComm
#'
#' @include ModelOrgClass.R
#'
#' @family Object classes
#' @exportClass ModelComm
setClass("ModelComm",
         slots = c(community = "data.frame"),
         prototype = list(
           community = data.frame(index = character(0L),
                                  id = character(0L),
                                  desc = character(0L),
                                  name = character(0L),
                                  abun = double(0L))
         ),
         contains = "ModelOrg")
