// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getGLPKVersion
Rcpp::CharacterVector getGLPKVersion();
RcppExport SEXP _cobrar_getGLPKVersion() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getGLPKVersion());
    return rcpp_result_gen;
END_RCPP
}
// initProb
SEXP initProb(const char* name);
RcppExport SEXP _cobrar_initProb(SEXP nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type name(nameSEXP);
    rcpp_result_gen = Rcpp::wrap(initProb(name));
    return rcpp_result_gen;
END_RCPP
}
// setObjDirLP
SEXP setObjDirLP(SEXP xp, int dir);
RcppExport SEXP _cobrar_setObjDirLP(SEXP xpSEXP, SEXP dirSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< int >::type dir(dirSEXP);
    rcpp_result_gen = Rcpp::wrap(setObjDirLP(xp, dir));
    return rcpp_result_gen;
END_RCPP
}
// addColsLP
SEXP addColsLP(SEXP xp, SEXP ncols);
RcppExport SEXP _cobrar_addColsLP(SEXP xpSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(addColsLP(xp, ncols));
    return rcpp_result_gen;
END_RCPP
}
// addRowsLP
SEXP addRowsLP(SEXP xp, SEXP nrows);
RcppExport SEXP _cobrar_addRowsLP(SEXP xpSEXP, SEXP nrowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nrows(nrowsSEXP);
    rcpp_result_gen = Rcpp::wrap(addRowsLP(xp, nrows));
    return rcpp_result_gen;
END_RCPP
}
// loadMatrixLP
SEXP loadMatrixLP(SEXP xp, SEXP ne, SEXP ia, SEXP ja, SEXP ra);
RcppExport SEXP _cobrar_loadMatrixLP(SEXP xpSEXP, SEXP neSEXP, SEXP iaSEXP, SEXP jaSEXP, SEXP raSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ne(neSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ia(iaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ja(jaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ra(raSEXP);
    rcpp_result_gen = Rcpp::wrap(loadMatrixLP(xp, ne, ia, ja, ra));
    return rcpp_result_gen;
END_RCPP
}
// setColsBndsObjCoefsLP
SEXP setColsBndsObjCoefsLP(SEXP xp, SEXP j, SEXP type, SEXP lb, SEXP ub, SEXP obj_coef);
RcppExport SEXP _cobrar_setColsBndsObjCoefsLP(SEXP xpSEXP, SEXP jSEXP, SEXP typeSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP obj_coefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type j(jSEXP);
    Rcpp::traits::input_parameter< SEXP >::type type(typeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< SEXP >::type obj_coef(obj_coefSEXP);
    rcpp_result_gen = Rcpp::wrap(setColsBndsObjCoefsLP(xp, j, type, lb, ub, obj_coef));
    return rcpp_result_gen;
END_RCPP
}
// setColsKindLP
SEXP setColsKindLP(SEXP xp, SEXP j, SEXP kind);
RcppExport SEXP _cobrar_setColsKindLP(SEXP xpSEXP, SEXP jSEXP, SEXP kindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type j(jSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kind(kindSEXP);
    rcpp_result_gen = Rcpp::wrap(setColsKindLP(xp, j, kind));
    return rcpp_result_gen;
END_RCPP
}
// setRowsBndsLP
SEXP setRowsBndsLP(SEXP xp, SEXP i, SEXP type, SEXP lb, SEXP ub);
RcppExport SEXP _cobrar_setRowsBndsLP(SEXP xpSEXP, SEXP iSEXP, SEXP typeSEXP, SEXP lbSEXP, SEXP ubSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    Rcpp::traits::input_parameter< SEXP >::type type(typeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ub(ubSEXP);
    rcpp_result_gen = Rcpp::wrap(setRowsBndsLP(xp, i, type, lb, ub));
    return rcpp_result_gen;
END_RCPP
}
// getSolStatLP
SEXP getSolStatLP(SEXP xp);
RcppExport SEXP _cobrar_getSolStatLP(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(getSolStatLP(xp));
    return rcpp_result_gen;
END_RCPP
}
// getColsPrimalLP
SEXP getColsPrimalLP(SEXP xp);
RcppExport SEXP _cobrar_getColsPrimalLP(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(getColsPrimalLP(xp));
    return rcpp_result_gen;
END_RCPP
}
// solveSimplex
SEXP solveSimplex(SEXP xp);
RcppExport SEXP _cobrar_solveSimplex(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(solveSimplex(xp));
    return rcpp_result_gen;
END_RCPP
}
// solveSimplexExact
SEXP solveSimplexExact(SEXP xp);
RcppExport SEXP _cobrar_solveSimplexExact(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(solveSimplexExact(xp));
    return rcpp_result_gen;
END_RCPP
}
// getObjVal
SEXP getObjVal(SEXP xp);
RcppExport SEXP _cobrar_getObjVal(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(getObjVal(xp));
    return rcpp_result_gen;
END_RCPP
}
// solveInterior
SEXP solveInterior(SEXP xp);
RcppExport SEXP _cobrar_solveInterior(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(solveInterior(xp));
    return rcpp_result_gen;
END_RCPP
}
// getObjValIpt
SEXP getObjValIpt(SEXP xp);
RcppExport SEXP _cobrar_getObjValIpt(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(getObjValIpt(xp));
    return rcpp_result_gen;
END_RCPP
}
// solveMIP
SEXP solveMIP(SEXP xp);
RcppExport SEXP _cobrar_solveMIP(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(solveMIP(xp));
    return rcpp_result_gen;
END_RCPP
}
// mipObjVal
SEXP mipObjVal(SEXP xp);
RcppExport SEXP _cobrar_mipObjVal(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(mipObjVal(xp));
    return rcpp_result_gen;
END_RCPP
}
// getColsDualIptLP
SEXP getColsDualIptLP(SEXP xp);
RcppExport SEXP _cobrar_getColsDualIptLP(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(getColsDualIptLP(xp));
    return rcpp_result_gen;
END_RCPP
}
// getColsDualLP
SEXP getColsDualLP(SEXP xp);
RcppExport SEXP _cobrar_getColsDualLP(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(getColsDualLP(xp));
    return rcpp_result_gen;
END_RCPP
}
// setMatRowLP
SEXP setMatRowLP(SEXP xp, SEXP i, SEXP len, SEXP ind, SEXP val);
RcppExport SEXP _cobrar_setMatRowLP(SEXP xpSEXP, SEXP iSEXP, SEXP lenSEXP, SEXP indSEXP, SEXP valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    Rcpp::traits::input_parameter< SEXP >::type len(lenSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ind(indSEXP);
    Rcpp::traits::input_parameter< SEXP >::type val(valSEXP);
    rcpp_result_gen = Rcpp::wrap(setMatRowLP(xp, i, len, ind, val));
    return rcpp_result_gen;
END_RCPP
}
// getNumRowsLP
SEXP getNumRowsLP(SEXP xp);
RcppExport SEXP _cobrar_getNumRowsLP(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(getNumRowsLP(xp));
    return rcpp_result_gen;
END_RCPP
}
// fvaLP
Rcpp::DataFrame fvaLP(SEXP xp, SEXP ind);
RcppExport SEXP _cobrar_fvaLP(SEXP xpSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(fvaLP(xp, ind));
    return rcpp_result_gen;
END_RCPP
}
// getSBMLVersion
Rcpp::String getSBMLVersion();
RcppExport SEXP _cobrar_getSBMLVersion() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getSBMLVersion());
    return rcpp_result_gen;
END_RCPP
}
// readSBMLfile
SEXP readSBMLfile(std::string file_path);
RcppExport SEXP _cobrar_readSBMLfile(SEXP file_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file_path(file_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(readSBMLfile(file_path));
    return rcpp_result_gen;
END_RCPP
}
// getModelObj
SEXP getModelObj(SEXP sbml_document_ptr);
RcppExport SEXP _cobrar_getModelObj(SEXP sbml_document_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sbml_document_ptr(sbml_document_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getModelObj(sbml_document_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getModelId
Rcpp::String getModelId(SEXP model_ptr);
RcppExport SEXP _cobrar_getModelId(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getModelId(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getModelName
Rcpp::String getModelName(SEXP model_ptr);
RcppExport SEXP _cobrar_getModelName(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getModelName(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getModelCompartments
Rcpp::DataFrame getModelCompartments(SEXP model_ptr);
RcppExport SEXP _cobrar_getModelCompartments(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getModelCompartments(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getStoichiometricMatrix
arma::sp_mat getStoichiometricMatrix(SEXP model_ptr);
RcppExport SEXP _cobrar_getStoichiometricMatrix(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getStoichiometricMatrix(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getModelCVTerms
Rcpp::StringVector getModelCVTerms(SEXP model_ptr);
RcppExport SEXP _cobrar_getModelCVTerms(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getModelCVTerms(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getModelSBOTerm
Rcpp::String getModelSBOTerm(SEXP model_ptr);
RcppExport SEXP _cobrar_getModelSBOTerm(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getModelSBOTerm(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getModelNotes
Rcpp::String getModelNotes(SEXP model_ptr);
RcppExport SEXP _cobrar_getModelNotes(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getModelNotes(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getObjectiveFunction
Rcpp::NumericVector getObjectiveFunction(SEXP model_ptr);
RcppExport SEXP _cobrar_getObjectiveFunction(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getObjectiveFunction(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getSubsystems
Rcpp::List getSubsystems(SEXP model_ptr);
RcppExport SEXP _cobrar_getSubsystems(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getSubsystems(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getReactionIds
Rcpp::CharacterVector getReactionIds(SEXP model_ptr);
RcppExport SEXP _cobrar_getReactionIds(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getReactionIds(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getReactionNames
Rcpp::CharacterVector getReactionNames(SEXP model_ptr);
RcppExport SEXP _cobrar_getReactionNames(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getReactionNames(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getReactionCompartment
Rcpp::CharacterVector getReactionCompartment(SEXP model_ptr);
RcppExport SEXP _cobrar_getReactionCompartment(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getReactionCompartment(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getReactionFluxBounds
Rcpp::List getReactionFluxBounds(SEXP model_ptr);
RcppExport SEXP _cobrar_getReactionFluxBounds(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getReactionFluxBounds(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getReactionCVTerms
Rcpp::List getReactionCVTerms(SEXP model_ptr);
RcppExport SEXP _cobrar_getReactionCVTerms(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getReactionCVTerms(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getReactionSBOTerms
Rcpp::StringVector getReactionSBOTerms(SEXP model_ptr);
RcppExport SEXP _cobrar_getReactionSBOTerms(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getReactionSBOTerms(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getMetaboliteIds
Rcpp::CharacterVector getMetaboliteIds(SEXP model_ptr);
RcppExport SEXP _cobrar_getMetaboliteIds(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getMetaboliteIds(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getMetaboliteNames
Rcpp::CharacterVector getMetaboliteNames(SEXP model_ptr);
RcppExport SEXP _cobrar_getMetaboliteNames(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getMetaboliteNames(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getMetaboliteAnnotation
Rcpp::DataFrame getMetaboliteAnnotation(SEXP model_ptr);
RcppExport SEXP _cobrar_getMetaboliteAnnotation(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getMetaboliteAnnotation(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getMetaboliteCVTerms
Rcpp::List getMetaboliteCVTerms(SEXP model_ptr);
RcppExport SEXP _cobrar_getMetaboliteCVTerms(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getMetaboliteCVTerms(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getMetaboliteSBOTerms
Rcpp::StringVector getMetaboliteSBOTerms(SEXP model_ptr);
RcppExport SEXP _cobrar_getMetaboliteSBOTerms(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getMetaboliteSBOTerms(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getMetaboliteCompartments
Rcpp::CharacterVector getMetaboliteCompartments(SEXP model_ptr);
RcppExport SEXP _cobrar_getMetaboliteCompartments(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getMetaboliteCompartments(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getGeneProducts
Rcpp::DataFrame getGeneProducts(SEXP model_ptr);
RcppExport SEXP _cobrar_getGeneProducts(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getGeneProducts(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getGeneProductCVTerms
Rcpp::List getGeneProductCVTerms(SEXP model_ptr);
RcppExport SEXP _cobrar_getGeneProductCVTerms(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getGeneProductCVTerms(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getGeneProductSBOTerms
Rcpp::StringVector getGeneProductSBOTerms(SEXP model_ptr);
RcppExport SEXP _cobrar_getGeneProductSBOTerms(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getGeneProductSBOTerms(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// getGPRs
Rcpp::List getGPRs(SEXP model_ptr);
RcppExport SEXP _cobrar_getGPRs(SEXP model_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_ptr(model_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(getGPRs(model_ptr));
    return rcpp_result_gen;
END_RCPP
}
// writeSBML
bool writeSBML(String file_path, String mod_id, String mod_name, String mod_desc, StringVector mod_cvterms, String mod_notes, int mod_sbo, StringVector comp_id, StringVector comp_name, StringVector met_id, StringVector met_name, NumericVector met_charge, StringVector met_formula, StringVector met_comp, Rcpp::ListOf<StringVector> met_cvterms, IntegerVector met_sbo, Rcpp::ListOf<IntegerVector> subsys, StringVector subsys_id, StringVector subsys_name, StringVector gene_id, StringVector gene_name, Rcpp::ListOf<StringVector> gene_cvterms, IntegerVector gene_sbo, StringVector param_id, NumericVector param_val, IntegerVector param_sbo, StringVector react_id, StringVector react_name, Rcpp::ListOf<NumericVector> Scoeff, Rcpp::ListOf<StringVector> react_mets, StringVector react_lb, StringVector react_ub, LogicalVector react_rev, Rcpp::ListOf<StringVector> react_cvterms, IntegerVector react_sbo, StringVector gpr, NumericVector obj_coef);
RcppExport SEXP _cobrar_writeSBML(SEXP file_pathSEXP, SEXP mod_idSEXP, SEXP mod_nameSEXP, SEXP mod_descSEXP, SEXP mod_cvtermsSEXP, SEXP mod_notesSEXP, SEXP mod_sboSEXP, SEXP comp_idSEXP, SEXP comp_nameSEXP, SEXP met_idSEXP, SEXP met_nameSEXP, SEXP met_chargeSEXP, SEXP met_formulaSEXP, SEXP met_compSEXP, SEXP met_cvtermsSEXP, SEXP met_sboSEXP, SEXP subsysSEXP, SEXP subsys_idSEXP, SEXP subsys_nameSEXP, SEXP gene_idSEXP, SEXP gene_nameSEXP, SEXP gene_cvtermsSEXP, SEXP gene_sboSEXP, SEXP param_idSEXP, SEXP param_valSEXP, SEXP param_sboSEXP, SEXP react_idSEXP, SEXP react_nameSEXP, SEXP ScoeffSEXP, SEXP react_metsSEXP, SEXP react_lbSEXP, SEXP react_ubSEXP, SEXP react_revSEXP, SEXP react_cvtermsSEXP, SEXP react_sboSEXP, SEXP gprSEXP, SEXP obj_coefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type file_path(file_pathSEXP);
    Rcpp::traits::input_parameter< String >::type mod_id(mod_idSEXP);
    Rcpp::traits::input_parameter< String >::type mod_name(mod_nameSEXP);
    Rcpp::traits::input_parameter< String >::type mod_desc(mod_descSEXP);
    Rcpp::traits::input_parameter< StringVector >::type mod_cvterms(mod_cvtermsSEXP);
    Rcpp::traits::input_parameter< String >::type mod_notes(mod_notesSEXP);
    Rcpp::traits::input_parameter< int >::type mod_sbo(mod_sboSEXP);
    Rcpp::traits::input_parameter< StringVector >::type comp_id(comp_idSEXP);
    Rcpp::traits::input_parameter< StringVector >::type comp_name(comp_nameSEXP);
    Rcpp::traits::input_parameter< StringVector >::type met_id(met_idSEXP);
    Rcpp::traits::input_parameter< StringVector >::type met_name(met_nameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type met_charge(met_chargeSEXP);
    Rcpp::traits::input_parameter< StringVector >::type met_formula(met_formulaSEXP);
    Rcpp::traits::input_parameter< StringVector >::type met_comp(met_compSEXP);
    Rcpp::traits::input_parameter< Rcpp::ListOf<StringVector> >::type met_cvterms(met_cvtermsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type met_sbo(met_sboSEXP);
    Rcpp::traits::input_parameter< Rcpp::ListOf<IntegerVector> >::type subsys(subsysSEXP);
    Rcpp::traits::input_parameter< StringVector >::type subsys_id(subsys_idSEXP);
    Rcpp::traits::input_parameter< StringVector >::type subsys_name(subsys_nameSEXP);
    Rcpp::traits::input_parameter< StringVector >::type gene_id(gene_idSEXP);
    Rcpp::traits::input_parameter< StringVector >::type gene_name(gene_nameSEXP);
    Rcpp::traits::input_parameter< Rcpp::ListOf<StringVector> >::type gene_cvterms(gene_cvtermsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type gene_sbo(gene_sboSEXP);
    Rcpp::traits::input_parameter< StringVector >::type param_id(param_idSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param_val(param_valSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type param_sbo(param_sboSEXP);
    Rcpp::traits::input_parameter< StringVector >::type react_id(react_idSEXP);
    Rcpp::traits::input_parameter< StringVector >::type react_name(react_nameSEXP);
    Rcpp::traits::input_parameter< Rcpp::ListOf<NumericVector> >::type Scoeff(ScoeffSEXP);
    Rcpp::traits::input_parameter< Rcpp::ListOf<StringVector> >::type react_mets(react_metsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type react_lb(react_lbSEXP);
    Rcpp::traits::input_parameter< StringVector >::type react_ub(react_ubSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type react_rev(react_revSEXP);
    Rcpp::traits::input_parameter< Rcpp::ListOf<StringVector> >::type react_cvterms(react_cvtermsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type react_sbo(react_sboSEXP);
    Rcpp::traits::input_parameter< StringVector >::type gpr(gprSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type obj_coef(obj_coefSEXP);
    rcpp_result_gen = Rcpp::wrap(writeSBML(file_path, mod_id, mod_name, mod_desc, mod_cvterms, mod_notes, mod_sbo, comp_id, comp_name, met_id, met_name, met_charge, met_formula, met_comp, met_cvterms, met_sbo, subsys, subsys_id, subsys_name, gene_id, gene_name, gene_cvterms, gene_sbo, param_id, param_val, param_sbo, react_id, react_name, Scoeff, react_mets, react_lb, react_ub, react_rev, react_cvterms, react_sbo, gpr, obj_coef));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_sbml_module();

static const R_CallMethodDef CallEntries[] = {
    {"_cobrar_getGLPKVersion", (DL_FUNC) &_cobrar_getGLPKVersion, 0},
    {"_cobrar_initProb", (DL_FUNC) &_cobrar_initProb, 1},
    {"_cobrar_setObjDirLP", (DL_FUNC) &_cobrar_setObjDirLP, 2},
    {"_cobrar_addColsLP", (DL_FUNC) &_cobrar_addColsLP, 2},
    {"_cobrar_addRowsLP", (DL_FUNC) &_cobrar_addRowsLP, 2},
    {"_cobrar_loadMatrixLP", (DL_FUNC) &_cobrar_loadMatrixLP, 5},
    {"_cobrar_setColsBndsObjCoefsLP", (DL_FUNC) &_cobrar_setColsBndsObjCoefsLP, 6},
    {"_cobrar_setColsKindLP", (DL_FUNC) &_cobrar_setColsKindLP, 3},
    {"_cobrar_setRowsBndsLP", (DL_FUNC) &_cobrar_setRowsBndsLP, 5},
    {"_cobrar_getSolStatLP", (DL_FUNC) &_cobrar_getSolStatLP, 1},
    {"_cobrar_getColsPrimalLP", (DL_FUNC) &_cobrar_getColsPrimalLP, 1},
    {"_cobrar_solveSimplex", (DL_FUNC) &_cobrar_solveSimplex, 1},
    {"_cobrar_solveSimplexExact", (DL_FUNC) &_cobrar_solveSimplexExact, 1},
    {"_cobrar_getObjVal", (DL_FUNC) &_cobrar_getObjVal, 1},
    {"_cobrar_solveInterior", (DL_FUNC) &_cobrar_solveInterior, 1},
    {"_cobrar_getObjValIpt", (DL_FUNC) &_cobrar_getObjValIpt, 1},
    {"_cobrar_solveMIP", (DL_FUNC) &_cobrar_solveMIP, 1},
    {"_cobrar_mipObjVal", (DL_FUNC) &_cobrar_mipObjVal, 1},
    {"_cobrar_getColsDualIptLP", (DL_FUNC) &_cobrar_getColsDualIptLP, 1},
    {"_cobrar_getColsDualLP", (DL_FUNC) &_cobrar_getColsDualLP, 1},
    {"_cobrar_setMatRowLP", (DL_FUNC) &_cobrar_setMatRowLP, 5},
    {"_cobrar_getNumRowsLP", (DL_FUNC) &_cobrar_getNumRowsLP, 1},
    {"_cobrar_fvaLP", (DL_FUNC) &_cobrar_fvaLP, 2},
    {"_cobrar_getSBMLVersion", (DL_FUNC) &_cobrar_getSBMLVersion, 0},
    {"_cobrar_readSBMLfile", (DL_FUNC) &_cobrar_readSBMLfile, 1},
    {"_cobrar_getModelObj", (DL_FUNC) &_cobrar_getModelObj, 1},
    {"_cobrar_getModelId", (DL_FUNC) &_cobrar_getModelId, 1},
    {"_cobrar_getModelName", (DL_FUNC) &_cobrar_getModelName, 1},
    {"_cobrar_getModelCompartments", (DL_FUNC) &_cobrar_getModelCompartments, 1},
    {"_cobrar_getStoichiometricMatrix", (DL_FUNC) &_cobrar_getStoichiometricMatrix, 1},
    {"_cobrar_getModelCVTerms", (DL_FUNC) &_cobrar_getModelCVTerms, 1},
    {"_cobrar_getModelSBOTerm", (DL_FUNC) &_cobrar_getModelSBOTerm, 1},
    {"_cobrar_getModelNotes", (DL_FUNC) &_cobrar_getModelNotes, 1},
    {"_cobrar_getObjectiveFunction", (DL_FUNC) &_cobrar_getObjectiveFunction, 1},
    {"_cobrar_getSubsystems", (DL_FUNC) &_cobrar_getSubsystems, 1},
    {"_cobrar_getReactionIds", (DL_FUNC) &_cobrar_getReactionIds, 1},
    {"_cobrar_getReactionNames", (DL_FUNC) &_cobrar_getReactionNames, 1},
    {"_cobrar_getReactionCompartment", (DL_FUNC) &_cobrar_getReactionCompartment, 1},
    {"_cobrar_getReactionFluxBounds", (DL_FUNC) &_cobrar_getReactionFluxBounds, 1},
    {"_cobrar_getReactionCVTerms", (DL_FUNC) &_cobrar_getReactionCVTerms, 1},
    {"_cobrar_getReactionSBOTerms", (DL_FUNC) &_cobrar_getReactionSBOTerms, 1},
    {"_cobrar_getMetaboliteIds", (DL_FUNC) &_cobrar_getMetaboliteIds, 1},
    {"_cobrar_getMetaboliteNames", (DL_FUNC) &_cobrar_getMetaboliteNames, 1},
    {"_cobrar_getMetaboliteAnnotation", (DL_FUNC) &_cobrar_getMetaboliteAnnotation, 1},
    {"_cobrar_getMetaboliteCVTerms", (DL_FUNC) &_cobrar_getMetaboliteCVTerms, 1},
    {"_cobrar_getMetaboliteSBOTerms", (DL_FUNC) &_cobrar_getMetaboliteSBOTerms, 1},
    {"_cobrar_getMetaboliteCompartments", (DL_FUNC) &_cobrar_getMetaboliteCompartments, 1},
    {"_cobrar_getGeneProducts", (DL_FUNC) &_cobrar_getGeneProducts, 1},
    {"_cobrar_getGeneProductCVTerms", (DL_FUNC) &_cobrar_getGeneProductCVTerms, 1},
    {"_cobrar_getGeneProductSBOTerms", (DL_FUNC) &_cobrar_getGeneProductSBOTerms, 1},
    {"_cobrar_getGPRs", (DL_FUNC) &_cobrar_getGPRs, 1},
    {"_cobrar_writeSBML", (DL_FUNC) &_cobrar_writeSBML, 37},
    {"_rcpp_module_boot_sbml_module", (DL_FUNC) &_rcpp_module_boot_sbml_module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cobrar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}