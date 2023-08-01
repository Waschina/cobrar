#include <iostream>
#include <RcppArmadillo.h>
#include <sbml/SBMLTypes.h>
#include <sbml/packages/fbc/common/FbcExtensionTypes.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
SEXP readSBMLfile(std::string file_path) {
  SBMLReader reader;
  SBMLDocument* document = reader.readSBML(file_path);

  if (document == nullptr) {
    Rcpp::stop("Error reading SBML document.");
  }

  return Rcpp::XPtr<SBMLDocument>(document, true);
}

// [[Rcpp::export]]
SEXP getModelObj(SEXP sbml_document_ptr) {
  SBMLDocument* document = Rcpp::XPtr<SBMLDocument>(sbml_document_ptr);

  if (document == nullptr) {
    Rcpp::stop("Invalid SBMLDocument pointer.");
  }

  Model* model = document->getModel();

  if (model == nullptr) {
    Rcpp::stop("No Model object found in the SBMLDocument.");
  }

  return Rcpp::XPtr<Model>(model, true);
}

// [[Rcpp::export]]
arma::sp_mat getStoichiometricMatrix(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  // Create the stoichiometric matrix (as a matrix or data frame) using the reaction and species information
  unsigned int num_reactions = model->getNumReactions();
  unsigned int num_species = model->getNumSpecies();
  arma::sp_mat S(num_species, num_reactions);

  for (unsigned int i = 0; i < num_reactions; i++) {
    Reaction* reaction = model->getReaction(i);
    for (unsigned int j = 0; j < num_species; j++) {
      Species* species = model->getSpecies(j);
      SpeciesReference* srefP = reaction->getProduct(species->getId());
      SpeciesReference* srefR = reaction->getReactant(species->getId());

      if(srefP != NULL) {
        S(j, i) = srefP->getStoichiometry();
      }
      else if(srefR != NULL) {
        S(j, i) = -srefR->getStoichiometry();
      }

    }
  }

  // Return the stoichiometric matrix
  return S; // Replace this with the actual stoichiometric matrix
}

// [[Rcpp::export]]
Rcpp::CharacterVector getReactionIds(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  unsigned int num_reactions = model->getNumReactions();
  Rcpp::CharacterVector rxn_id;

  for(unsigned int i = 0; i < num_reactions; i++) {
    Reaction* reaction = model->getReaction(i);
    rxn_id.push_back(reaction->getId());
  }

  return rxn_id;
}

// [[Rcpp::export]]
Rcpp::CharacterVector getReactionNames(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  unsigned int num_reactions = model->getNumReactions();
  Rcpp::CharacterVector rxn_name;

  for(unsigned int i = 0; i < num_reactions; i++) {
    Reaction* reaction = model->getReaction(i);
    rxn_name.push_back(reaction->getName());
  }

  return rxn_name;
}


// [[Rcpp::export]]
Rcpp::List getReactionFluxBounds(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  Rcpp::NumericVector lower_bounds;
  Rcpp::NumericVector upper_bounds;

  ListOfParameters* lpara = model->getListOfParameters();


  // Loop over the reactions in the model to extract their flux bounds
  for (unsigned int i = 0; i < model->getNumReactions(); i++) {
    Reaction* reaction = model->getReaction(i);

    FbcReactionPlugin* rplugin = static_cast<FbcReactionPlugin*>(reaction->getPlugin("fbc"));
    std::string lower_bnd = rplugin->getLowerFluxBound();
    std::string upper_bnd = rplugin->getUpperFluxBound();

    double lower_bndval = (lpara->get(lower_bnd))->getValue();
    double upper_bndval = (lpara->get(upper_bnd))->getValue();

    lower_bounds.push_back(lower_bndval);
    upper_bounds.push_back(upper_bndval);

  }

  // Create a named list to return reaction ids and flux bounds
  Rcpp::List reaction_flux_bounds;
  reaction_flux_bounds["lower_bound"] = lower_bounds;
  reaction_flux_bounds["upper_bound"] = upper_bounds;

  return reaction_flux_bounds;
}

RCPP_MODULE(sbml_module) {
  function("readSBMLfile", &readSBMLfile, "Read SBML document using libSBML");
  function("getModelObj", &getModelObj, "Get the Model object from the SBMLDocument");
  function("getStoichiometricMatrix", &getStoichiometricMatrix, "Gets the models stoichiomentric matrix S.");
  function("getReactionIds", &getReactionIds, "Gets reaction Ids.");
  function("getReactionNames", &getReactionNames, "Gets reaction Ids.");
  function("getReactionFluxBounds", &getReactionFluxBounds, "Gets reaction Ids.");
}
