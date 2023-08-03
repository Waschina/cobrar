# cobrar

#### Notes
The SBML standard with its fbc extension allows to specify more than one objective (class `ListOfObjectives`). However, cobrar can only handle one current objective function per model, which is defined as an objective coefficient vector in slot `obj_coef` of an object of class `modelorg`. Note that when reading SBML models, cobrar will only use the first objective defined in the SBML document.
