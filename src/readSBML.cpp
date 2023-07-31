#include <iostream>
#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>
#include <util.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
LIBSBML_CPP_NAMESPACE_USE
BEGIN_C_DECLS

// [[Rcpp::export]]
int readSBML (char filename)
{
  SBMLDocument* document;
  SBMLReader reader;
#ifdef __BORLANDC__
  unsigned long start, stop;
#else
  unsigned long long start, stop;
#endif
  start    = getCurrentMillis();
  document = reader.readSBML(filename);
  stop     = getCurrentMillis();
  unsigned int errors = document->getNumErrors();
  cout << endl;
  cout << "            filename: " << filename              << endl;
  cout << "           file size: " << getFileSize(filename) << endl;
  cout << "      read time (ms): " << stop - start          << endl;
  cout << " validation error(s): " << errors << endl;
  cout << endl;
  document->printErrors(cerr);
  delete document;
  return errors;
}
END_C_DECLS
