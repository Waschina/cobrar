#include <iostream>
#include <RcppArmadillo.h>
#include <glpk.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/* structure for glpk parameters */
glp_smcp parmS;
glp_iptcp parmI;
glp_iocp parmM;

/*
 * Get glpk version number
 */
// [[Rcpp::export]]
Rcpp::CharacterVector getGLPKVersion() {
  return glp_version();
}

// Forward declaration of glp_prob
struct glp_prob;

// Delete GLP problem
void lpXPtrFinalizer(SEXP lp_ptr) {
  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(lp_ptr);
  glp_delete_prob(lp);
}

// [[Rcpp::export]]
SEXP initProb(const char* name, double tol_bnd)
{
  glp_prob *lp;

  glp_init_smcp(&parmS);
  glp_init_iptcp(&parmI);
  glp_init_iocp(&parmM);

  parmS.msg_lev = GLP_MSG_ERR;
  parmI.msg_lev = GLP_MSG_ERR;
  parmM.msg_lev = GLP_MSG_ERR;

  // set feasibility/bound tolerance for simplex
  parmS.tol_bnd = tol_bnd;

  lp = glp_create_prob();

  if (lp == nullptr) {
    Rcpp::stop("Error initializing LP (glpk).");
  }

  glp_set_prob_name(lp, name);

  // if (strcmp("max", dir) == 0) {
  //   glp_set_obj_dir(lp, GLP_MAX);
  // }
  // if (strcmp("min", dir) == 0) {
  //   glp_set_obj_dir(lp, GLP_MIN);
  // }

  SEXP xp = R_MakeExternalPtr(lp, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(xp, lpXPtrFinalizer);

  return xp;
}

// [[Rcpp::export]]
SEXP setObjDirLP(SEXP xp, int dir)
{
  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  SEXP out = R_NilValue;

  glp_set_obj_dir(lp, dir);

  return out;
}

// [[Rcpp::export]]
SEXP addColsLP(SEXP xp, SEXP ncols)
{
  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  SEXP out = R_NilValue;
  int fcol = 0;

  fcol = glp_add_cols(lp, Rf_asInteger(ncols));

  out = Rf_ScalarInteger(fcol);

  return out;
}

// [[Rcpp::export]]
SEXP addRowsLP(SEXP xp, SEXP nrows)
{
  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  SEXP out = R_NilValue;
  int frow = 0;

  frow = glp_add_rows(lp, Rf_asInteger(nrows));

  out = Rf_ScalarInteger(frow);

  return out;
}

// [[Rcpp::export]]
SEXP loadMatrixLP(SEXP xp,
                  SEXP ne, SEXP ia, SEXP ja, SEXP ra) {

  SEXP out = R_NilValue;

  const int *ria = INTEGER(ia);
  const int *rja = INTEGER(ja);
  const double *rra = REAL(ra);

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  glp_load_matrix(lp, Rf_asInteger(ne),
                  &(ria[-1]), &(rja[-1]), &(rra[-1]));

  return out;
}

// [[Rcpp::export]]
SEXP setColsBndsObjCoefsLP(SEXP xp, SEXP j, SEXP type,
                           SEXP lb, SEXP ub, SEXP obj_coef) {

  SEXP out = R_NilValue;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  int k, nj;
  int *rj = INTEGER(j);
  double *rlb = REAL(lb), *rub = REAL(ub), *robj_coef = REAL(obj_coef);
  const int *rtype;

  if (type == R_NilValue) {
    rtype = NULL;
  }
  else {
    rtype = INTEGER(type);
  }

  nj = Rf_length(j);

  if (rtype == NULL) {
    for (k = 0; k < nj; k++) {
      if (std::islessgreater(rlb[k],rub[k])) {
        glp_set_col_bnds(lp,
                         rj[k], GLP_DB, rlb[k], rub[k]
        );
      }
      else {
        glp_set_col_bnds(lp,
                         rj[k], GLP_FX, rlb[k], rub[k]
        );
      }
      glp_set_obj_coef(lp, rj[k], robj_coef[k]);
    }
  }
  else {
    for (k = 0; k < nj; k++) {
      glp_set_col_bnds(lp,
                       rj[k], rtype[k], rlb[k], rub[k]
      );
      glp_set_obj_coef(lp, rj[k], robj_coef[k]);
    }

  }

  return out;

}

// [[Rcpp::export]]
SEXP setColsKindLP(SEXP xp, SEXP j, SEXP kind) {

  SEXP out = R_NilValue;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  int k, nj;
  int *rj = INTEGER(j);
  int *rkind = INTEGER(kind);


  nj = Rf_length(j);
  for (k = 0; k < nj; k++) {
    glp_set_col_kind(lp, rj[k], rkind[k]);
  }

  return out;

}

// [[Rcpp::export]]
SEXP setRowsBndsLP(SEXP xp, SEXP i, SEXP type, SEXP lb, SEXP ub) {

  SEXP out = R_NilValue;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  int k, ni;
  int *ri = INTEGER(i);
  double *rlb = REAL(lb), *rub = REAL(ub);
  const int *rtype;

  if (type == R_NilValue) {
    rtype = NULL;
  }
  else {
    rtype = INTEGER(type);
  }

  ni = Rf_length(i);

  if (rtype == NULL) {
    for (k = 0; k < ni; k++) {
      if (std::islessgreater(rlb[k], rub[k])) {
        glp_set_row_bnds(lp,
                         ri[k], GLP_DB, rlb[k], rub[k]
        );
      }
      else {
        glp_set_row_bnds(lp,
                         ri[k], GLP_FX, rlb[k], rub[k]
        );
      }
    }
  }
  else {
    for (k = 0; k < ni; k++) {
      glp_set_row_bnds(lp,
                       ri[k], rtype[k], rlb[k], rub[k]
      );
    }
  }

  return out;

}

// [[Rcpp::export]]
SEXP getSolStatLP(SEXP xp) {

  SEXP out = R_NilValue;
  int stat = 0;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  stat = glp_get_status(lp);

  out = Rf_ScalarInteger(stat);

  return out;
}

// [[Rcpp::export]]
SEXP getColsPrimalLP(SEXP xp) {

  std::vector<double> prim;
  int stat = 0;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  unsigned int nc = glp_get_num_cols(lp);

  // in case of infeasible solution
  stat = glp_get_status(lp);
  if(stat != 2 && stat != 5) {
    return DoubleVector(nc, DoubleVector::get_na());
  }

  for(unsigned int i=1; i <= nc; i++) {
    prim.push_back(glp_get_col_prim(lp, i));
  }

  return Rcpp::wrap(prim);
}


SEXP setDefaultSmpParm() {

  SEXP parmext = R_NilValue;

  glp_init_smcp(&parmS);

  return parmext;
}

SEXP setDefaultIptParm() {

  SEXP parmext = R_NilValue;

  glp_init_iptcp(&parmI);

  return parmext;
}

SEXP setDefaultMIPParm() {

  SEXP parmext = R_NilValue;

  glp_init_iocp(&parmM);

  return parmext;
}


/*
 * Method-specific functions
 */

/* Simplex */
// [[Rcpp::export]]
SEXP solveSimplex(SEXP xp) {

  SEXP out = R_NilValue;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  int ret = 0;

  ret = glp_simplex(lp, &parmS);

  out = Rf_ScalarInteger(ret);

  return out;
}

// [[Rcpp::export]]
SEXP solveSimplexExact(SEXP xp) {

  SEXP out = R_NilValue;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  int ret = 0;

  ret = glp_exact(lp, &parmS);

  out = Rf_ScalarInteger(ret);

  return out;
}

static int hook(void *info, const char *s)
{
  return 1;
}

// [[Rcpp::export]]
SEXP scaleSimplexProb(SEXP xp) {

  SEXP out = R_NilValue;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  glp_term_hook(hook, NULL); // redirects output to somewhere NULL
  glp_scale_prob(lp, GLP_SF_SKIP);
  glp_term_hook(NULL, NULL); // get output back to console

  return out;
}

// [[Rcpp::export]]
SEXP getObjVal(SEXP xp) {

  SEXP out = R_NilValue;
  double obj_val = 0;
  int stat = 0;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  // in case of infeasible solution
  stat = glp_get_status(lp);
  if(stat != 2 && stat != 5) {
    return DoubleVector(1, DoubleVector::get_na());
  }

  //in case of feasible/optimal
  obj_val = glp_get_obj_val(lp);

  out = Rf_ScalarReal(obj_val);

  return out;
}

/* Interior */
// [[Rcpp::export]]
SEXP solveInterior(SEXP xp) {

  SEXP out = R_NilValue;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  int ret = 0;

  ret = glp_interior(lp, &parmI);

  out = Rf_ScalarInteger(ret);

  return out;
}

// [[Rcpp::export]]
SEXP getObjValIpt(SEXP xp) {

  SEXP out = R_NilValue;
  double obj_val = 0;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  obj_val = glp_ipt_obj_val(lp);

  out = Rf_ScalarReal(obj_val);

  return out;
}

/* MIP */
// [[Rcpp::export]]
SEXP solveMIP(SEXP xp) {

  SEXP out = R_NilValue;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  int ret = 0;

  ret = glp_intopt(lp, &parmM);

  out = Rf_ScalarInteger(ret);

  return out;
}

// [[Rcpp::export]]
SEXP mipObjVal(SEXP xp) {

  SEXP out = R_NilValue;
  double obj_val = 0;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  obj_val = glp_mip_obj_val(lp);

  out = Rf_ScalarReal(obj_val);

  return out;
}

/* retrieve column dual value (interior) for all columns (interior) */
// [[Rcpp::export]]
SEXP getColsDualIptLP(SEXP xp) {

  SEXP out = R_NilValue;
  double col_dual = 0;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  int num_cols, k;

  num_cols = glp_get_num_cols(lp);

  out = Rf_allocVector(REALSXP, num_cols);
  for (k = 1; k <= num_cols; k++) {
    col_dual = glp_ipt_col_dual(lp, k);
    REAL(out)[k-1] = col_dual;
  }

  return out;
}

/* retrieve column dual value for all columns */
// [[Rcpp::export]]
SEXP getColsDualLP(SEXP xp) {

  SEXP out = R_NilValue;
  double col_dual = 0;
  int stat = 0;

  int num_cols, k;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  num_cols = glp_get_num_cols(lp);

  // in case of infeasible solution
  stat = glp_get_status(lp);
  if(stat != 2 && stat != 5) {
    return DoubleVector(num_cols, DoubleVector::get_na());
  }

  // in case of feasible/optimal solution
  out = Rf_allocVector(REALSXP, num_cols);
  for (k = 1; k <= num_cols; k++) {
    col_dual = glp_get_col_dual(lp, k);
    REAL(out)[k-1] = col_dual;
  }

  return out;
}

/* set or replace row of constraint matrix */
// [[Rcpp::export]]
SEXP setMatRowLP(SEXP xp, SEXP i, SEXP len, SEXP ind, SEXP val) {
  SEXP out = R_NilValue;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  int *rind = INTEGER(ind);
  double *rval = REAL(val);

  glp_set_mat_row(lp, Rf_asInteger(i), Rf_asInteger(len), &(rind[-1]), &(rval[-1]));

  return out;
}

// /* set specific values in constraint matrix */
// // [[Rcpp::export]]
// SEXP setMatEntriesLP(SEXP xp, IntegerVector ci, IntegerVector cj, NumericVector val) {
//   SEXP out = R_NilValue;
//
//   glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);
//
//   int n = ci.size();
//
//
//   for(int i = 1; i <= n; i ++) {
//     int row_ind = ci[i-1];
//     int col_ind[] = {cj[i-1]};
//     double ival[] = {val[i-1]};
//     glp_set_mat_row(lp, row_ind, 1, col_ind, ival);
//   }
//
//   return out;
// }

/* get number of rows */
// [[Rcpp::export]]
SEXP getNumRowsLP(SEXP xp) {

  SEXP out = R_NilValue;
  int nrows = 0;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  nrows = glp_get_num_rows(lp);

  out = Rf_ScalarInteger(nrows);

  return out;
}


/* wrapper for FVA */
// [[Rcpp::export]]
Rcpp::DataFrame fvaLP(SEXP xp, SEXP ind) {

  // Get the indices as integers
  IntegerVector indices(ind);

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  // overwrite current objective function
  for(int i = 1; i <= glp_get_num_cols(lp); i++) {
    glp_set_obj_coef(lp, i, 0.0);
  }

  // Create vectors to store results
  NumericVector minVals(indices.size());
  NumericVector maxVals(indices.size());

  // MAX
  glp_set_obj_dir(lp, GLP_MAX);
  for(unsigned int i = 0; i < indices.size(); i++) {
    int columnIndex = indices[i];
    glp_set_obj_coef(lp, columnIndex, 1.0);

    // Maximize the variable
    glp_simplex(lp, &parmS);

    // Store the maximum value
    maxVals[i] = glp_get_obj_val(lp);

    // Reset the objective coefficient to 0
    glp_set_obj_coef(lp, columnIndex, 0.0);
  }

  // MIN
  glp_set_obj_dir(lp, GLP_MIN);
  for(unsigned int i = 0; i < indices.size(); i++) {
    int columnIndex = indices[i];
    glp_set_obj_coef(lp, columnIndex, 1.0);

    // Minimize the variable
    glp_simplex(lp, &parmS);

    // Store the minimum value
    minVals[i] = glp_get_obj_val(lp);

    // Reset the objective coefficient to 0
    glp_set_obj_coef(lp, columnIndex, 0.0);
  }

  // Create a DataFrame to store the results
  DataFrame result = DataFrame::create(_["min.flux"] = minVals, _["max.flux"] = maxVals);

  return result;
}


