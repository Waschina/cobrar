#include <Rcpp.h>
#include "sbml/SBMLTypes.h"
#include "sbml/SBMLDocument.h"

// Function to read an SBML file and return a pointer to the SBMLDocument
// [[Rcpp::export]]
Rcpp::XPtr<SBMLDocument> readSBMLFile(const std::string& filename) {
  SBMLReader reader;
  SBMLDocument* doc = reader.readSBMLFromFile(filename);

  if (doc == nullptr) {
    Rcpp::stop("Failed to read SBML file: " + filename);
  }

  Rcpp::XPtr<SBMLDocument> xptr(doc, true);
  return xptr;
}
