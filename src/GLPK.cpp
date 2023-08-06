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
SEXP initProb(const char* name)
{
  glp_prob *lp;

  glp_init_smcp(&parmS);
  glp_init_iptcp(&parmI);
  glp_init_iocp(&parmM);

  parmS.msg_lev = GLP_MSG_ERR;
  parmI.msg_lev = GLP_MSG_ERR;
  parmM.msg_lev = GLP_MSG_ERR;

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
SEXP setObjDir(SEXP xp, int dir)
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

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

  unsigned int nc = glp_get_num_cols(lp);

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

// [[Rcpp::export]]
SEXP getObjVal(SEXP xp) {

  SEXP out = R_NilValue;
  double obj_val = 0;

  glp_prob* lp = (glp_prob*)R_ExternalPtrAddr(xp);

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
