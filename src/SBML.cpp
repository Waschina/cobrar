#include <iostream>
#include <vector>
#include <RcppArmadillo.h>
#include <sbml/SBMLTypes.h>
#include <sbml/packages/fbc/common/FbcExtensionTypes.h>
#include <sbml/packages/groups/common/GroupsExtensionTypes.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/*
 * Get libSBML version number
 */
// [[Rcpp::export]]
Rcpp::String getSBMLVersion() {
  return getLibSBMLDottedVersion();
}

/*
 * Create pointers visible to R
 */

// [[Rcpp::export]]
SEXP readSBMLfile(std::string file_path) {
  SBMLReader reader;
  SBMLDocument* document = reader.readSBML(file_path);

  if (document == nullptr) {
    Rcpp::stop("Error reading SBML document.");
  }

  return Rcpp::XPtr<SBMLDocument>(document, false);
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

  return Rcpp::XPtr<Model>(model, false);
}

/*
 * Model fields
 */

// [[Rcpp::export]]
Rcpp::String getModelId(SEXP model_ptr) {
  // Get the Model object
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  std::string modId = model->getId();
  if (modId.size() == 0)
    return NA_STRING;

  return modId;
}

// [[Rcpp::export]]
Rcpp::String getModelName(SEXP model_ptr) {
  // Get the Model object
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  std::string modName = model->getName();
  if (modName.size() == 0)
    return NA_STRING;

  return modName;
}

// [[Rcpp::export]]
Rcpp::DataFrame getModelCompartments(SEXP model_ptr) {
  // Get the Model object
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  ListOfCompartments* complist = model->getListOfCompartments();
  Rcpp::CharacterVector compIds;
  Rcpp::CharacterVector compNames;

  for (unsigned int i = 0; i < model->getNumCompartments(); i++) {
    String cid = complist->get(i)->getId();
    String cname = complist->get(i)->getName();
    if(cname == "")
      cname = NA_STRING;

    compIds.push_back(cid);
    compNames.push_back(cname);
  }

  Rcpp::List cols = Rcpp::List::create(Named("id") = compIds,
                                       Named("name") = compNames);
  Rcpp::DataFrame df(cols);

  return df;
}

// [[Rcpp::export]]
arma::sp_mat getStoichiometricMatrix(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  // Create the stoichiometric matrix
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
Rcpp::String getModelAnnotation(SEXP model_ptr) {
  // Get the Model object
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  std::string modAnno = model->getAnnotationString();
  if (modAnno.size() == 0)
    return NA_STRING;

  return modAnno;
}

// [[Rcpp::export]]
Rcpp::String getModelNotes(SEXP model_ptr) {
  // Get the Model object
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  std::string modNotes = model->getNotesString();
  if (modNotes.size() == 0)
    return NA_STRING;

  return modNotes;
}

// [[Rcpp::export]]
Rcpp::NumericVector getObjectiveFunction(SEXP model_ptr) {
  // Get the Model object
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  Rcpp::NumericVector objCoeff(model->getNumReactions());
  std::fill(objCoeff.begin(), objCoeff.end(), 0);

  FbcModelPlugin* mplugin = static_cast<FbcModelPlugin*>(model->getPlugin("fbc"));
  ListOfObjectives* objList = mplugin->getListOfObjectives();
  ListOfReactions* reaList = model->getListOfReactions();

  if(objList->getNumObjectives() > 0) {
    Objective* objFunc = objList->get(0); // only first objective is retrieved
    ListOfFluxObjectives* objFluxes = objFunc->getListOfFluxObjectives();
    for(unsigned int i = 0; i < objFluxes->getNumFluxObjectives(); i++) {
      FluxObjective* flxObj = objFluxes->get(i);
      std::string objrea = flxObj->getReaction();

      int targetReactionIndex = -1;

      for (unsigned int i = 0; i < reaList->size(); ++i) {
        if (reaList->get(i)->getId() == objrea) {
          targetReactionIndex = i;
          break; // Exit the loop once found
        }
      }

      if(targetReactionIndex >= 0) {
        double objCoeffVal = flxObj->getCoefficient();
        objCoeff[targetReactionIndex] = objCoeffVal;
      }
    }
  }

  return(objCoeff);
}

// [[Rcpp::export]]
Rcpp::List getSubsystems(SEXP model_ptr) {
  // Get the Model object
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  SBasePlugin* grpPlugin = model->getPlugin("groups");
  unsigned int num_reactions = model->getNumReactions();

  if(grpPlugin != NULL) {
    GroupsModelPlugin* mplugin = static_cast<GroupsModelPlugin*>(grpPlugin);
    ListOfGroups* grpList = mplugin->getListOfGroups();

    unsigned int num_groups = mplugin->getNumGroups();
    arma::sp_mat subSys(num_reactions, num_groups);

    Rcpp::CharacterVector subsysIds;
    Rcpp::CharacterVector subsysNames;

    // create char vector with reaction ids
    std::vector<std::string> rxn_id;
    for(unsigned int i = 0; i < num_reactions; i++) {
      Reaction* reaction = model->getReaction(i);
      rxn_id.push_back(reaction->getId());
    }
    // Create an unordered map for indexing
    std::unordered_map<std::string, size_t> indexMap;
    for(size_t i = 0; i < rxn_id.size(); ++i) {
      indexMap[rxn_id[i]] = i;
    }

    // iterate through subsystems
    for(unsigned int i = 0; i < num_groups; i++) {
      Group* iGrp = grpList->get(i);
      subsysIds.push_back(iGrp->getId());
      subsysNames.push_back(iGrp->getName());

      ListOfMembers* memberList = iGrp->getListOfMembers();
      for(unsigned int j = 0; j < memberList->getNumMembers(); j++) {
        Member* iRea = memberList->get(j);
        std::string iReaId = iRea->getIdRef();
        if (indexMap.find(iReaId) != indexMap.end()) {
          size_t targetIndex = indexMap[iReaId];
          subSys(targetIndex, i) = true;
        }
      }
    }
    Rcpp::List subSysList;
    subSysList["subSys"] = subSys;
    subSysList["subSys_ids"] = subsysIds;
    subSysList["subSys_names"] = subsysNames;
    return subSysList;
  }
  Rcpp::List subSysList;
  subSysList["subSys"] = arma::sp_mat(num_reactions,0);
  subSysList["subSys_ids"] = Rcpp::CharacterVector(0);
  subSysList["subSys_names"] = Rcpp::CharacterVector(0);
  return subSysList;
}

/*
 * Reaction fields
 */

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
Rcpp::CharacterVector getReactionAnnotation(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  unsigned int num_reactions = model->getNumReactions();
  Rcpp::CharacterVector rxn_anno;

  for(unsigned int i = 0; i < num_reactions; i++) {
    Reaction* reaction = model->getReaction(i);
    rxn_anno.push_back(reaction->getAnnotationString());
  }

  return rxn_anno;
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


/*
 * Metabolite fields
 */

// [[Rcpp::export]]
Rcpp::CharacterVector getMetaboliteIds(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  unsigned int num_metabolites = model->getNumSpecies();
  Rcpp::CharacterVector met_id;

  for(unsigned int i = 0; i < num_metabolites; i++) {
    Species* species = model->getSpecies(i);
    met_id.push_back(species->getId());
  }

  return met_id;
}


// [[Rcpp::export]]
Rcpp::CharacterVector getMetaboliteNames(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  unsigned int num_metabolites = model->getNumSpecies();
  Rcpp::CharacterVector met_name;

  for(unsigned int i = 0; i < num_metabolites; i++) {
    Species* species = model->getSpecies(i);
    met_name.push_back(species->getName());
  }

  return met_name;
}

// [[Rcpp::export]]
Rcpp::DataFrame getMetaboliteAnnotation(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  unsigned int num_metabolites = model->getNumSpecies();
  Rcpp::CharacterVector met_chemForm;
  Rcpp::DoubleVector met_charge;
  Rcpp::CharacterVector met_anno;

  for(unsigned int i = 0; i < num_metabolites; i++) {
    Species* species = model->getSpecies(i);

    FbcSpeciesPlugin* splugin = static_cast<FbcSpeciesPlugin*>(species->getPlugin("fbc"));

    met_chemForm.push_back(splugin->getChemicalFormula());
    met_charge.push_back(splugin->getCharge());
    met_anno.push_back(species->getAnnotationString());
  }

  Rcpp::List cols = Rcpp::List::create(Named("chemicalFormula") = met_chemForm,
                                       Named("charge") = met_charge,
                                       Named("annotation") = met_anno);
  Rcpp::DataFrame df(cols);

  return df;
}


// [[Rcpp::export]]
Rcpp::DataFrame getGeneProducts(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  FbcModelPlugin* mplugin = static_cast<FbcModelPlugin*>(model->getPlugin("fbc"));

  ListOfGeneProducts* lgenes = mplugin->getListOfGeneProducts();
  Rcpp::CharacterVector geneID;
  Rcpp::CharacterVector geneName;

  for(unsigned int i = 0; i < mplugin->getNumGeneProducts(); i++) {
    geneID.push_back(lgenes->get(i)->getId());
    geneName.push_back(lgenes->get(i)->getName());
  }

  Rcpp::List cols = Rcpp::List::create(Named("ID") = geneID,
                                       Named("name") = geneName);
  Rcpp::DataFrame df(cols);

  return df;
}

// small recursive helper function to retrieve GPR-association logical string
std::string getGPRString(FbcAssociation* fasso,
                         std::map<std::string, std::string>& variableMap,
                         int& variableCount) {
  if(fasso->isGeneProductRef()) {

    GeneProductRef* gp = static_cast<GeneProductRef*>(fasso);
    std::string gpID = gp->getGeneProduct();

    if (variableMap.find(gpID) != variableMap.end()) {
      return variableMap[gpID];
    }

    std::string newVariable = "x[" + std::to_string(variableCount++) + "]";
    variableMap[gpID] = newVariable;
    return newVariable;

  } else if(fasso->isFbcAnd()) {

    FbcAnd* fand = static_cast<FbcAnd*>(fasso);
    ListOfFbcAssociations* assoList = fand->getListOfAssociations();

    std::vector<std::string> chfasso;
    std::string gpString;
    for(unsigned int i = 0; i < assoList->getNumFbcAssociations(); i ++) {
      std::string tmpstr = getGPRString(assoList->get(i),
                                        variableMap, variableCount);
      if(i > 0)
        gpString += " & ";
      gpString += tmpstr;
    }
    return "( " + gpString + " )";

  } else {

    FbcOr* fand = static_cast<FbcOr*>(fasso);
    ListOfFbcAssociations* assoList = fand->getListOfAssociations();

    std::vector<std::string> chfasso;
    std::string gpString;
    for(unsigned int i = 0; i < assoList->getNumFbcAssociations(); i ++) {
      std::string tmpstr = getGPRString(assoList->get(i),
                                        variableMap, variableCount);
      if(i > 0)
        gpString += " | ";
      gpString += tmpstr;
    }
    return "( " + gpString + " )";

  }
}


// [[Rcpp::export]]
Rcpp::List getGPRs(SEXP model_ptr) {
  // Get the Model object from the SBMLDocument
  Model* model = Rcpp::XPtr<Model>(model_ptr);

  if (model == nullptr) {
    Rcpp::stop("Invalid Model pointer.");
  }

  Rcpp::CharacterVector igrp;
  Rcpp::List rgenes;

  for(unsigned int i = 0; i < model->getNumReactions(); i++) {
    Reaction* reaction = model->getReaction(i);

    std::string tmpstr = "";
    FbcReactionPlugin* rplugin = static_cast<FbcReactionPlugin*>(reaction->getPlugin("fbc"));
    GeneProductAssociation* fasso = rplugin->getGeneProductAssociation();

    if(fasso != nullptr) {
      std::map<std::string, std::string> variableMap;
      int variableCount = 1; // Initialize the variable count
      tmpstr = getGPRString(fasso->getAssociation(),variableMap,variableCount);

      std::vector<std::string> orderedKeys;
      for (const auto& pair : variableMap) {
        orderedKeys.push_back(pair.first);
      }

      rgenes.push_back(orderedKeys);
    }

    igrp.push_back(tmpstr);
  }

  Rcpp::List cols = Rcpp::List::create(Named("rules") = igrp,
                                       Named("genes") = rgenes);

  return cols;
}

RCPP_MODULE(sbml_module) {
  function("readSBMLfile", &readSBMLfile, "Read SBML document using libSBML");
  function("getModelObj", &getModelObj, "Get the Model object from the SBMLDocument");
}
