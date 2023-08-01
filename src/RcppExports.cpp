// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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

RcppExport SEXP _rcpp_module_boot_sbml_module();

static const R_CallMethodDef CallEntries[] = {
    {"_cobrar_readSBMLfile", (DL_FUNC) &_cobrar_readSBMLfile, 1},
    {"_cobrar_getModelObj", (DL_FUNC) &_cobrar_getModelObj, 1},
    {"_cobrar_getStoichiometricMatrix", (DL_FUNC) &_cobrar_getStoichiometricMatrix, 1},
    {"_cobrar_getReactionIds", (DL_FUNC) &_cobrar_getReactionIds, 1},
    {"_cobrar_getReactionNames", (DL_FUNC) &_cobrar_getReactionNames, 1},
    {"_cobrar_getReactionFluxBounds", (DL_FUNC) &_cobrar_getReactionFluxBounds, 1},
    {"_rcpp_module_boot_sbml_module", (DL_FUNC) &_rcpp_module_boot_sbml_module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cobrar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
