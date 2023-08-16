# some GLPK specific codes:
glpkPar <- list(
  GLP_FR = 1,
  GLP_LO = 2,
  GLP_UP = 3,
  GLP_DB = 4,
  GLP_FX = 5,

  GLP_MIN = 1,
  GLP_MAX = 2,

  GLP_CV = 1,
  GLP_IV = 2,
  GLP_BV = 3,

  GLP_EBADB = 1,
  GLP_ESING = 2,
  GLP_ECOND = 3,
  GLP_EBOUND = 4,
  GLP_EFAIL = 5,
  GLP_EOBJLL = 6,
  GLP_EOBJUL = 7,
  GLP_EITLIM = 8,
  GLP_ETMLIM = 9,
  GLP_ENOPFS = 10,
  GLP_ENODFS = 11,
  GLP_EROOT = 12,
  GLP_ESTOP = 13,
  GLP_EMIPGAP = 14,
  GLP_ENOFEAS = 15,
  GLP_ENOCVG = 16,
  GLP_EINSTAB = 17,
  GLP_EDATA = 18,
  GLP_ERANGE = 19,

  GLP_OPT = 5,
  GLP_UNDEF = 1,
  GLP_FEAS = 2,
  GLP_INFEAS = 3,
  GLP_NOFEAS = 4,
  GLP_UNBND = 6

)

setClass(Class = "LPproblem_glpk",
         contains = "LPproblem"
)

#------------------------------------------------------------------------------#
# Function instances for the default solver: GLPK                              #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature = "LPproblem_glpk",
          definition = function(.Object,
                                name,
                                method) {

            .Object@ptr <- initProb(name)
            .Object@solver = "glpk"
            .Object@method = method

            return(.Object)
          }
)

# function for cplex could look something like this:
# setMethod(f = "initialize",
#           signature = "LPproblem_cplex",
#           definition = function(.Object,
#                                 name,
#                                 dir = COBRAR_SETTINGS("OPT_DIRECTION")) {
#
#             .Object@ptr <- cobrarCPLEX::initProb(name, dir)
#             .Object@solver = "cplex"
#
#             return(.Object)
#           }
# )


setMethod("loadLPprob", signature(lp = "LPproblem_glpk"),

          function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype, lpdir,
                   rub = NULL, ctype = NULL) {

            # optimization direction
            lpdir <- switch(EXPR = lpdir,
                            "max" = glpkPar$GLP_MAX,
                            "min" = glpkPar$GLP_MIN)
            setObjDirection(lp, lpdir = lpdir)

            # problem dimensions
            addCols(lp, ncols = nCols)
            addRows(lp, nrows = nRows)


            # constraint matrix
            TMPmat <- as(mat, "TsparseMatrix")
            loadMatrix(lp,
                       ne = length(TMPmat@x),
                       ia = TMPmat@i + 1,
                       ja = TMPmat@j + 1,
                       ra = TMPmat@x)

            # column (variable) bounds and objective function
            setColsBndsObjCoefs(lp,
                                j = c(1:nCols),
                                lb = lb,
                                ub = ub,
                                obj_coef = obj)


            # variable type
            if (!is.null(ctype)) {
              cctype <- sapply(ctype,
                               function(x) switch(EXPR = x,
                                                  "C" = glpkPar$GLP_CV,
                                                  "I" = glpkPar$GLP_IV,
                                                  "B" = glpkPar$GLP_BV,
                                                        glpkPar$GLP_CV))

              setColsKind(lp, j = c(1:nCols), kind = cctype)
            }

            # right hand side
            # Note: This is the steady state condition: Production of internal
            # metabolites should equal the consumption. In over words:
            # row lower bound = row upper bounds = 0
            if (is.null(rub)) {
              # The values in rlb will be copied to rub. GLPK ignores rlb and rub,
              # depending on the constraint type (e.g. an upper bound, if the
              # constraint type says, it has a lower bound):
              # Constraint type "L": ignore rub
              # Constraint type "U": ignore rlb
              # Constraint type "E": ignore rub
              # Constraint type "F": ignore rlb and rub

              crub <- rlb
            }
            else {
              crub <- rub
            }
            stopifnot(length(rlb) == length(crub))
            setRowsBnds(lp,
                        i = c(1:nRows),
                        lb = rlb,
                        ub = crub,
                        type = rtype)


          }
)

setMethod("setObjDirection", signature(lp = "LPproblem_glpk"),
          function(lp, lpdir) {
            setObjDir(lp@ptr, lpdir)
          }
)

setMethod("addCols", signature(lp = "LPproblem_glpk"),
          function(lp, ncols) {
            addColsLP(lp@ptr, as.integer(ncols))
          }
)

setMethod("addRows", signature(lp = "LPproblem_glpk"),
          function(lp, nrows) {
            addRowsLP(lp@ptr, as.integer(nrows))
          }
)

setMethod("loadMatrix", signature(lp = "LPproblem_glpk"),
          function(lp, ne, ia, ja, ra) {
            loadMatrixLP(lp@ptr,
                         as.integer(ne),
                         as.integer(ia),
                         as.integer(ja),
                         as.numeric(ra))
          }
)

setMethod("setColsBndsObjCoefs", signature(lp = "LPproblem_glpk"),
          function(lp, j, lb, ub, obj_coef, type = NULL) {


            if (is.null(type)) {
              Ctype <- as.null(type)
            }
            else {
              Ctype <- as.integer(type)
            }

            setColsBndsObjCoefsLP(lp@ptr,
                                  as.integer(j),
                                  Ctype,
                                  as.numeric(lb),
                                  as.numeric(ub),
                                  as.numeric(obj_coef))
          }
)

setMethod("setColsKind", signature(lp = "LPproblem_glpk"),
          function(lp, j, kind) {
            setColsKindLP(lp@ptr,
                          as.integer(j),
                          as.integer(kind))
          }
)


setMethod("setRowsBnds", signature(lp = "LPproblem_glpk"),
          function(lp, i, lb, ub , type) {


            type <- sapply(type,
                           function(x) switch(EXPR = x,
                                              "F" = glpkPar$GLP_FR,
                                              "L" = glpkPar$GLP_LO,
                                              "U" = glpkPar$GLP_UP,
                                              "D" = glpkPar$GLP_DB,
                                              "E" = glpkPar$GLP_FX,
                                              glpkPar$GLP_FX))

            if (is.null(type)) {
              Ctype <- as.null(type)
            }
            else {
              Ctype <- as.integer(type)
            }

            setRowsBndsLP(lp@ptr,
                          as.integer(i),
                          Ctype,
                          as.numeric(lb),
                          as.numeric(ub))

          }
)

setMethod("solveLp", signature(lp = "LPproblem_glpk"),
          function(lp) {
            out          <- FALSE

            switch(EXPR = lp@method,
                   "interior" = {
                     out <- solveInterior(lp@ptr)
                   },
                   "exact" = {
                     out <- solveSimplexExact(lp@ptr)
                   },
                   "mip" = {
                     out <- solveMIP(lp@ptr)
                   },
                   {
                     out <- solveSimplex(lp@ptr)
                   }
            )

            # get optimization status term for code in "out"
            term <- switch(EXPR = out+1,
                           "optimization process was successful",
                           "invalid basis",
                           "singular matrix",
                           "ill-conditioned matrix",
                           "invalid bounds",
                           "solver failed",
                           "objective lower limit reached",
                           "objective upper limit reached",
                           "iteration limit exceeded",
                           "time limit exceeded",
                           "no primal feasible solution",
                           "no dual feasible solution",
                           "root LP optimum not provided",
                           "search terminated by application",
                           "relative mip gap tolerance reached",
                           "no primal/dual feasible solution",
                           "no convergence",
                           "numerical instability",
                           "invalid data",
                           "result out of range")
            if(is.null(out))
              term <- paste("Failed to obtain solution, unknown error code:", out)

            return(list(code= out,
                        term = term))
          }
)

setMethod("getObjValue", signature(lp = "LPproblem_glpk"),
          function(lp) {
            switch(EXPR = lp@method,
                   "interior" = {
                     out <- getObjValIpt(lp@ptr)
                   },
                   "mip" = {
                     out <- mipObjVal(lp@ptr)
                   },
                   {
                     out <- getObjVal(lp@ptr)
                   }
            )

          }
)

setMethod("getSolStat", signature(lp = "LPproblem_glpk"),
          function(lp) {

            out <- getSolStatLP(lp@ptr)

            # get term
            term <- switch(EXPR = out,
                           "solution is undefined",
                           "solution is feasible",
                           "olution is infeasible",
                           "problem has no feasible solution",
                           "solution is optimal",
                           "problem has unbounded solution")
            if(is.null(out))
              term <- paste("unknown status code:", out)

            return(list(code = out,
                        term = term))
          }
)

setMethod("getColsPrimal", signature(lp = "LPproblem_glpk"),
          function(lp) {

            out <- getColsPrimalLP(lp@ptr)

            return(out)
          }
)

setMethod("getRedCosts", signature(lp = "LPproblem_glpk"),
          function(lp) {

            if (lp@method == "interior") {
              out <- getColsDualIptLP(lp@ptr)
            }
            else {
              out <- getColsDualLP(lp@ptr)
            }

            return(out)
          }
)

setMethod("addSingleConstraint", signature(lp = "LPproblem_glpk"),
          function(lp, coeffs, lb, ub, type) {

            # add new row to constraint matrix
            addRows(lp, nrows = 1)
            i_newrow <- getNumRowsLP(lp@ptr)

            nz <- which(coeffs != 0)

            # add coeffs to new row
            setMatRowLP(lp@ptr,
                        as.integer(i_newrow),
                        as.integer(length(nz)),
                        as.integer(nz),
                        as.numeric(coeffs[nz]))

            # set bounds
            setRowsBnds(lp,
                        i = i_newrow,
                        lb = lb,
                        ub = ub,
                        type = type)


          }
)

setMethod("fvaJob", signature(lp = "LPproblem_glpk"),
          function(lp, ind) {

            fvares <- fvaLP(lp@ptr,
                            as.integer(ind))

            return(fvares)
          }
)
