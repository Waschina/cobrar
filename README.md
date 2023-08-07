# cobrar









#### Key differences to `sybil`

- cobrar is fully functional from reading SBML files until optimisation of linear programs, without the need of additional packages such as "sybilSBML" or "glpkAPI".
- The GNU glpk library is a system requirement. There are no links/references to iBM's cplex in the cobrar package, which were one reason why sybil was discontinued on CRAN.
- [roxygen2](https://roxygen2.r-lib.org/) is used for documenting functions and classes.
- Feature trim: A range of functionalities in sybil are not part of cobrar. 
- In cobrar, R's garbage collector handles memory management, including memory associated to C++-objects and the pointers to these.
- Simplifications in class and function architecture
  - No more Class "SysBiolAlg" nor derived classes. Different algorithms have their own function and detailed documentation of their return variables.
- Performance (i.e., computation time) improvements in certain procedures:
  - identification of dead-end metabolites
  - reading SBML files,
  - pFBA algorithm (a.k.a MTF 'Minimization of Total Flux').


#### Notes

*Multiple objectives*. The SBML standard with its fbc extension allows to specify more than one objective (class `ListOfObjectives`). However, cobrar can only handle one current objective function per model, which is defined as an objective coefficient vector in slot `obj_coef` of an object of class `modelorg`. Note that when reading SBML models, cobrar will only use the first objective defined in the SBML document.
