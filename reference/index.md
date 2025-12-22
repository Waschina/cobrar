# Package index

## Model import and export helpers

- [`readSBMLmod()`](https://waschina.github.io/cobrar/reference/readSBMLmod.md)
  : Reads an SBML file and constructs an object of class 'ModelOrg'
- [`readSybilmod()`](https://waschina.github.io/cobrar/reference/readSybilmod.md)
  : Reads a sybil model file and constructs an object of cobrar's class
  'ModelOrg'
- [`writeSBMLmod()`](https://waschina.github.io/cobrar/reference/writeSBMLmod.md)
  : Exports a Metabolic Network in SBML Format
- [`frog()`](https://waschina.github.io/cobrar/reference/frog-methods.md)
  : FROG: Reproducible fitness and robustness diagnostics for
  constraint-based models

## Model manipulation tools

- [`addCompartment()`](https://waschina.github.io/cobrar/reference/addCompartment.md)
  : Add compartments or update their data
- [`addConstraint()`](https://waschina.github.io/cobrar/reference/addConstraint-methods.md)
  : Add constraints to model
- [`addGene()`](https://waschina.github.io/cobrar/reference/addGene.md)
  : Add genes or update their data
- [`addMetabolite()`](https://waschina.github.io/cobrar/reference/addMetabolite.md)
  : Add metabolites or update their data
- [`addReact()`](https://waschina.github.io/cobrar/reference/addReact.md)
  : Add or modify a reaction
- [`addSubsystem()`](https://waschina.github.io/cobrar/reference/addSubsystem.md)
  : Add subsystems or update their data
- [`changeBounds()`](https://waschina.github.io/cobrar/reference/changeBounds.md)
  : Change flux bounds
- [`changeObjFunc()`](https://waschina.github.io/cobrar/reference/changeObjFunc.md)
  : Change the objective function
- [`setObjDir()`](https://waschina.github.io/cobrar/reference/setObjDir.md)
  : Set objective direction
- [`rmCompartment()`](https://waschina.github.io/cobrar/reference/rmCompartment.md)
  : Remove compartments from a model
- [`rmConstraint()`](https://waschina.github.io/cobrar/reference/rmConstraint.md)
  : Remove constraints
- [`rmGene()`](https://waschina.github.io/cobrar/reference/rmGene.md) :
  Remove genes from a model
- [`rmMetabolite()`](https://waschina.github.io/cobrar/reference/rmMetabolite.md)
  : Remove metabolites from a model
- [`rmReact()`](https://waschina.github.io/cobrar/reference/rmReact.md)
  : Remove reactions from a model
- [`rmSubsystem()`](https://waschina.github.io/cobrar/reference/rmSubsystem.md)
  : Remove subsystems from a model

## Model characteristics

- [`checkCompartmentId()`](https://waschina.github.io/cobrar/reference/checkCompartmentId.md)
  : Check compartment IDs and Indices
- [`checkGeneId()`](https://waschina.github.io/cobrar/reference/checkGeneId.md)
  : Check gene IDs and Indices
- [`checkMetId()`](https://waschina.github.io/cobrar/reference/checkMetId.md)
  : Check metabolite IDs and Indices
- [`checkReactId()`](https://waschina.github.io/cobrar/reference/checkReactId.md)
  : Check reaction IDs and Indices
- [`checkSubsystemId()`](https://waschina.github.io/cobrar/reference/checkSubsystemId.md)
  : Check subsystem IDs and Indices
- [`comp_num()`](https://waschina.github.io/cobrar/reference/comp_num-methods.md)
  : Number of compartments
- [`comp_pos()`](https://waschina.github.io/cobrar/reference/comp_pos-methods.md)
  : Index of compartment(s)
- [`constraint_num()`](https://waschina.github.io/cobrar/reference/constraint_num-methods.md)
  : Number of constraints
- [`met_num()`](https://waschina.github.io/cobrar/reference/met_num-methods.md)
  : Number of metabolites
- [`met_pos()`](https://waschina.github.io/cobrar/reference/met_pos-methods.md)
  : Index of metabolite(s)
- [`gene_num()`](https://waschina.github.io/cobrar/reference/gene_num-methods.md)
  : Number of genes
- [`gene_pos()`](https://waschina.github.io/cobrar/reference/gene_pos-methods.md)
  : Index of gene(s)
- [`geneDel()`](https://waschina.github.io/cobrar/reference/geneDel.md)
  : Identify reactions affected by gene knockouts
- [`deadEndMetabolites()`](https://waschina.github.io/cobrar/reference/deadEndMetabolites.md)
  : Identify dead end metabolites
- [`findExchReact()`](https://waschina.github.io/cobrar/reference/findExchReact-methods.md)
  : Find exchange reactions
- [`printConstraint()`](https://waschina.github.io/cobrar/reference/printConstraint.md)
  : Print Constraint(s)
- [`printReaction()`](https://waschina.github.io/cobrar/reference/printReaction.md)
  : Print reaction(s)
- [`react_num()`](https://waschina.github.io/cobrar/reference/react_num-methods.md)
  : Number of reactions
- [`react_pos()`](https://waschina.github.io/cobrar/reference/react_pos-methods.md)
  : Index of reaction(s)
- [`show(`*`<ModelOrg>`*`)`](https://waschina.github.io/cobrar/reference/show-ModelOrg-method.md)
  : Print a short summary of a metabolic network model
- [`subsys_num()`](https://waschina.github.io/cobrar/reference/subsys_num-methods.md)
  : Number of subsystems
- [`subsys_pos()`](https://waschina.github.io/cobrar/reference/subsys_pos-methods.md)
  : Index of subsystem(s)
- [`guessBMReaction()`](https://waschina.github.io/cobrar/reference/guessBMReaction-methods.md)
  : Guess biomass reaction(s)

## Flux prediction algorithms

- [`fba()`](https://waschina.github.io/cobrar/reference/fba.md) : Flux
  Balance Analysis
- [`fva()`](https://waschina.github.io/cobrar/reference/fva.md) : Flux
  Variability Analysis (FVA)
- [`pfba()`](https://waschina.github.io/cobrar/reference/pfba.md) :
  Parsimonious Flux Balance Analysis (pFBA)
- [`pfbaHeuristic()`](https://waschina.github.io/cobrar/reference/pfbaHeuristic.md)
  : Heuristic parsimonious Flux Balance Analysis (pFBA)

## Flux analysis tools

- [`getExchanges()`](https://waschina.github.io/cobrar/reference/getExchanges.md)
  : Get metabolite exchange rates
- [`show(`*`<FluxPrediction>`*`)`](https://waschina.github.io/cobrar/reference/show-FluxPrediction-method.md)
  : Print a short summary of a flux prediction result

## Community model functions

- [`joinModels()`](https://waschina.github.io/cobrar/reference/joinModels.md)
  : Join multiple metabolic models to form a community
- [`fixBMRatios()`](https://waschina.github.io/cobrar/reference/fixBMRatios.md)
  : Fixate biomass ratios
- [`fluxBMCoupling()`](https://waschina.github.io/cobrar/reference/fluxBMCoupling.md)
  : Couple reaction flux bounds to biomass production

## Chemoinformatics

- [`countElements()`](https://waschina.github.io/cobrar/reference/countElements.md)
  : Count elements in formulas
- [`elements()`](https://waschina.github.io/cobrar/reference/elements.md)
  : Data frame of elements and their average weights
- [`mass()`](https://waschina.github.io/cobrar/reference/mass.md) :
  Calculate molar mass of molecules

## Object classes

- [`Constraints-class`](https://waschina.github.io/cobrar/reference/Constraints-class.md)
  [`Constraints`](https://waschina.github.io/cobrar/reference/Constraints-class.md)
  : Structure of Constraints Class
- [`FluxPrediction-class`](https://waschina.github.io/cobrar/reference/FluxPrediction-class.md)
  [`FluxPrediction`](https://waschina.github.io/cobrar/reference/FluxPrediction-class.md)
  : Structure of FluxPrediction Class
- [`ModelOrg-class`](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)
  [`ModelOrg`](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)
  : Structure of ModelOrg Class
- [`ModelComm-class`](https://waschina.github.io/cobrar/reference/ModelComm-class.md)
  [`ModelComm`](https://waschina.github.io/cobrar/reference/ModelComm-class.md)
  : Structure of ModelComm Class

## Setup

- [`COBRAR_SETTINGS()`](https://waschina.github.io/cobrar/reference/COBRAR_SETTINGS.md)
  : Set and get central cobrar parameters

## Internal LP solver connectors

- [`addCols()`](https://waschina.github.io/cobrar/reference/addCols-methods.md)
  : Add columns to LP problem
- [`addRows()`](https://waschina.github.io/cobrar/reference/addRows-methods.md)
  : Add rows to LP problem
- [`addSingleConstraint()`](https://waschina.github.io/cobrar/reference/addSingleConstraint-methods.md)
  : Add single constraint
- [`deleteLP()`](https://waschina.github.io/cobrar/reference/deleteLP-methods.md)
  : Delete an LP problem
- [`fvaJob()`](https://waschina.github.io/cobrar/reference/fvaJob-methods.md)
  : Wrapper function for efficient FVA
- [`getColsPrimal()`](https://waschina.github.io/cobrar/reference/getColsPrimal-methods.md)
  : Retrieve column primal value
- [`getObjValue()`](https://waschina.github.io/cobrar/reference/getObjValue-methods.md)
  : Get the objective value of a solved LP problem
- [`getRedCosts()`](https://waschina.github.io/cobrar/reference/getRedCosts-methods.md)
  : Retrieve column reduced costs
- [`getSolStat()`](https://waschina.github.io/cobrar/reference/getSolStat-methods.md)
  : Get the solver status
- [`loadLPprob()`](https://waschina.github.io/cobrar/reference/loadLPprob-methods.md)
  : Initialize a LP problem
- [`loadMatrix()`](https://waschina.github.io/cobrar/reference/loadMatrix-methods.md)
  : Populate a constraint-X-variable matrix
- [`setColsBndsObjCoefs()`](https://waschina.github.io/cobrar/reference/setColsBndsObjCoefs-methods.md)
  : Set column bounds and objective coefficients
- [`setColsKind()`](https://waschina.github.io/cobrar/reference/setColsKind-methods.md)
  : Set column types
- [`setObjDirection()`](https://waschina.github.io/cobrar/reference/setObjDirection-methods.md)
  : Set objective direction
- [`setRowsBnds()`](https://waschina.github.io/cobrar/reference/setRowsBnds-methods.md)
  : Set row bounds
- [`solveLp()`](https://waschina.github.io/cobrar/reference/solveLp-methods.md)
  : Solve an LP problem
- [`LPproblem-class`](https://waschina.github.io/cobrar/reference/LPproblem-class.md)
  [`LPproblem`](https://waschina.github.io/cobrar/reference/LPproblem-class.md)
  : Structure of LPproblem Class
- [`LPproblem_glpk-class`](https://waschina.github.io/cobrar/reference/LPproblem_glpk-class.md)
  : Structure of LPproblem_glpk Class
