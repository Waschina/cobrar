# cobrar



#### Current TODOs

- [ ] Model modification functions
  - [x] addReact
  - [x] rmReact
  - [x] addMetabolite
  - [x] rmMetabolite
  - [ ] addComp
  - [ ] rmComp (no prio)
  - [x] addGene
  - [x] rmGene
  - [ ] addSubSystem
  - [ ] rmSubSystem
- [x] gene KOs
- [x] FVA (with adjustable proportions)
- [x] FVA documentation
- [ ] MOMA (iterative linearized objective?) -> quadratic opt not possible with glpk
- [ ] ROOM (mixed-integer optimization)
- [x] additional constraints as slot in `modelorg`
  - [x] Slot in modelorg
  - [x] rm if reaction of a constraint is removed
  - [x] addconstraint function
  - [x] rm constraint function
  - [x] constraints taking effect (glpk  LP formulation)
    - [x] FBA
    - [x] pFBA
    - [x] pFBA-heuristic
    - [x] FVA

  - [x] full documentation of function `addConstraint`, `rmConstraint`
- [x] `writeSBMLmod`
  - [x] CVTerms
    - [x] read reaction CVTerms
    - [x] write reaction CVTerms
    - [x] read metabolite CVTerms
    - [x] write metabolite CVTerms
    - [x] read geneProduct CVTerms
    - [x] write geneProduct CVTerms
    - [x] read model CVTerms
    - [x] write model CVTerms
  - [x] SBOTerms
    - [x] read reaction SBOTerms
    - [x] write reaction SBOTerms
    - [x] read metabolite SBOTerms
    - [x] write metabolite SBOTerms
    - [x] read geneProduct SBOTerms
    - [x] write geneProduct SBOTerms
    - [x] read model SBOTerm
    - [x] write model SBOTerm
  - [x] Model notes
  - [x] Subsystems/groups
  - [x] GPRs
  - [x] Objective function
- [x] export and documentation of position functions
- [x] export and documentation of count functions (genes, reactions, constraints, metabolites)
- [ ] read SBML documents of level 2 (necessary?)
- [x] define a class for FBA (and FBA-related) prediction results






#### Key differences to `sybil`

- cobrar is fully functional from reading SBML files until optimisation of linear programs, without the need of additional packages such as "sybilSBML" or "glpkAPI".
- The GNU glpk library is a system requirement. There are no links/references to iBM's cplex in the cobrar package, which were one reason why sybil was discontinued on CRAN.
- [roxygen2](https://roxygen2.r-lib.org/) is used for documenting functions and classes.
- Feature trim: A range of functionalities in sybil are not part of cobrar. 
- In cobrar, R's garbage collector handles memory management, including memory associated to C++-objects and the pointers to these.
- Simplifications in class and function architecture
  - No more Class "SysBiolAlg" nor sub-classes. Different algorithms have their own function and detailed documentation of their return variables.
- In sybil, columns names "annotation" were actually concatenated CVTerms (https://synonym.caltech.edu/software/libsbml/5.20.0/cpp-api/class_c_v_term.html). To avoid confusion with other levels of annotation, the columns (e.g. in `react_attr` or `met_attr`) are named 'CVTerms' in cobrar.
- Performance (i.e., computation time) improvements in certain procedures:
  - identification of dead-end metabolites
  - reading/exporting SBML files,
  - pFBA algorithm (a.k.a MTF 'Minimization of Total Flux').
  - FVA; also now allows relaxed constraints on optimal growth (e.g. flux variability with 90-100% optimal growth)


#### Notes

- cobrar exports SBML files level 3 version 2 with `fbc` version 3 and `groups` version 1. 
- Group assignments are only supported for reactions.
- currently only glpk is supported as solver. A plugin for IBM's ILOG cplex is planned.
- *Multiple objectives*. The SBML standard with its fbc extension allows to specify more than one objective (class `ListOfObjectives`). However, cobrar can only handle one current objective function per model, which is defined as an objective coefficient vector in slot `obj_coef` of an object of class `modelorg`. Note that when reading SBML models, cobrar will only use the first objective defined in the SBML document.
