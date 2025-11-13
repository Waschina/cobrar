# Package index

## Model import and export helpers

- [`readSBMLmod()`](readSBMLmod.md) : Reads an SBML file and constructs
  an object of class 'ModelOrg'
- [`readSybilmod()`](readSybilmod.md) : Reads a sybil model file and
  constructs an object of cobrar's class 'ModelOrg'
- [`writeSBMLmod()`](writeSBMLmod.md) : Exports a Metabolic Network in
  SBML Format

## Model manipulation tools

- [`addCompartment()`](addCompartment.md) : Add compartments or update
  their data
- [`addConstraint()`](addConstraint-methods.md) : Add constraints to
  model
- [`addGene()`](addGene.md) : Add genes or update their data
- [`addMetabolite()`](addMetabolite.md) : Add metabolites or update
  their data
- [`addReact()`](addReact.md) : Add or modify a reaction
- [`addSubsystem()`](addSubsystem.md) : Add subsystems or update their
  data
- [`changeBounds()`](changeBounds.md) : Change flux bounds
- [`changeObjFunc()`](changeObjFunc.md) : Change the objective function
- [`setObjDir()`](setObjDir.md) : Set objective direction
- [`rmCompartment()`](rmCompartment.md) : Remove compartments from a
  model
- [`rmConstraint()`](rmConstraint.md) : Remove constraints
- [`rmGene()`](rmGene.md) : Remove genes from a model
- [`rmMetabolite()`](rmMetabolite.md) : Remove metabolites from a model
- [`rmReact()`](rmReact.md) : Remove reactions from a model
- [`rmSubsystem()`](rmSubsystem.md) : Remove subsystems from a model

## Model characteristics

- [`checkCompartmentId()`](checkCompartmentId.md) : Check compartment
  IDs and Indices
- [`checkGeneId()`](checkGeneId.md) : Check gene IDs and Indices
- [`checkMetId()`](checkMetId.md) : Check metabolite IDs and Indices
- [`checkReactId()`](checkReactId.md) : Check reaction IDs and Indices
- [`checkSubsystemId()`](checkSubsystemId.md) : Check subsystem IDs and
  Indices
- [`comp_num()`](comp_num-methods.md) : Number of compartments
- [`comp_pos()`](comp_pos-methods.md) : Index of compartment(s)
- [`constraint_num()`](constraint_num-methods.md) : Number of
  constraints
- [`met_num()`](met_num-methods.md) : Number of metabolites
- [`met_pos()`](met_pos-methods.md) : Index of metabolite(s)
- [`gene_num()`](gene_num-methods.md) : Number of genes
- [`gene_pos()`](gene_pos-methods.md) : Index of gene(s)
- [`geneDel()`](geneDel.md) : Identify reactions affected by gene
  knockouts
- [`deadEndMetabolites()`](deadEndMetabolites.md) : Identify dead end
  metabolites
- [`findExchReact()`](findExchReact-methods.md) : Find exchange
  reactions
- [`printConstraint()`](printConstraint.md) : Print Constraint(s)
- [`printReaction()`](printReaction.md) : Print reaction(s)
- [`react_num()`](react_num-methods.md) : Number of reactions
- [`react_pos()`](react_pos-methods.md) : Index of reaction(s)
- [`show(`*`<ModelOrg>`*`)`](show-ModelOrg-method.md) : Print a short
  summary of a metabolic network model
- [`subsys_num()`](subsys_num-methods.md) : Number of subsystems
- [`subsys_pos()`](subsys_pos-methods.md) : Index of subsystem(s)
- [`guessBMReaction()`](guessBMReaction-methods.md) : Guess biomass
  reaction(s)

## Flux prediction algorithms

- [`fba()`](fba.md) : Flux Balance Analysis
- [`fva()`](fva.md) : Flux Variability Analysis (FVA)
- [`pfba()`](pfba.md) : Parsimonious Flux Balance Analysis (pFBA)
- [`pfbaHeuristic()`](pfbaHeuristic.md) : Heuristic parsimonious Flux
  Balance Analysis (pFBA)

## Flux analysis tools

- [`getExchanges()`](getExchanges.md) : Get metabolite exchange rates
- [`show(`*`<FluxPrediction>`*`)`](show-FluxPrediction-method.md) :
  Print a short summary of a flux prediction result

## Community model functions

- [`joinModels()`](joinModels.md) : Join multiple metabolic models to
  form a community
- [`fixBMRatios()`](fixBMRatios.md) : Fixate biomass ratios
- [`fluxBMCoupling()`](fluxBMCoupling.md) : Couple reaction flux bounds
  to biomass production

## Chemoinformatics

- [`countElements()`](countElements.md) : Count elements in formulas
- [`elements()`](elements.md) : Data frame of elements and their average
  weights
- [`mass()`](mass.md) : Calculate molar mass of molecules

## Object classes

- [`Constraints-class`](Constraints-class.md)
  [`Constraints`](Constraints-class.md) : Structure of Constraints Class
- [`FluxPrediction-class`](FluxPrediction-class.md)
  [`FluxPrediction`](FluxPrediction-class.md) : Structure of
  FluxPrediction Class
- [`ModelOrg-class`](ModelOrg-class.md) [`ModelOrg`](ModelOrg-class.md)
  : Structure of ModelOrg Class
- [`ModelComm-class`](ModelComm-class.md)
  [`ModelComm`](ModelComm-class.md) : Structure of ModelComm Class

## Setup

- [`COBRAR_SETTINGS()`](COBRAR_SETTINGS.md) : Set and get central cobrar
  parameters

## Internal LP solver connectors

- [`addCols()`](addCols-methods.md) : Add columns to LP problem
- [`addRows()`](addRows-methods.md) : Add rows to LP problem
- [`addSingleConstraint()`](addSingleConstraint-methods.md) : Add single
  constraint
- [`deleteLP()`](deleteLP-methods.md) : Delete an LP problem
- [`fvaJob()`](fvaJob-methods.md) : Wrapper function for efficient FVA
- [`getColsPrimal()`](getColsPrimal-methods.md) : Retrieve column primal
  value
- [`getObjValue()`](getObjValue-methods.md) : Get the objective value of
  a solved LP problem
- [`getRedCosts()`](getRedCosts-methods.md) : Retrieve column reduced
  costs
- [`getSolStat()`](getSolStat-methods.md) : Get the solver status
- [`loadLPprob()`](loadLPprob-methods.md) : Initialize a LP problem
- [`loadMatrix()`](loadMatrix-methods.md) : Populate a
  constraint-X-variable matrix
- [`setColsBndsObjCoefs()`](setColsBndsObjCoefs-methods.md) : Set column
  bounds and objective coefficients
- [`setColsKind()`](setColsKind-methods.md) : Set column types
- [`setObjDirection()`](setObjDirection-methods.md) : Set objective
  direction
- [`setRowsBnds()`](setRowsBnds-methods.md) : Set row bounds
- [`solveLp()`](solveLp-methods.md) : Solve an LP problem
- [`LPproblem-class`](LPproblem-class.md)
  [`LPproblem`](LPproblem-class.md) : Structure of LPproblem Class
- [`LPproblem_glpk-class`](LPproblem_glpk-class.md) : Structure of
  LPproblem_glpk Class
