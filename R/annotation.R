guessReactionSBOTerm <- function(id, met, Scoef, metComp, metChemicalFormula) {
  # List of relevant SBO terms for reactions
  #
  # Biological:
  # * SBO:0000167 - biochemical or transport reaction
  # * SBO:0000176 - biochemical reaction
  # * SBO:0000185 - translocation reaction (Movement of a physical entity without modification of the structure of the entity.)
  # * SBO:0000655 - transport reaction (The movement of an entity/entities across a biological membrane mediated by a transporter protein.)
  #
  # Technical:
  # * SBO:0000627 - Exchange reaction
  # * SBO:0000628 - Demand reaction
  # * SBO:0000632 - Sink reaction (SW: usually also just categorized as demand reaction)

  # Exchange/Demand/Sink reactions
  if(grepl("^EX_",id))
    return("SBO:0000627")
  if(grepl("^DM_",id))
    return("SBO:0000628")
  if(grepl("^SK_",id))
    return("SBO:0000632")

  # transport or transformation or both?
  if(!any(is.na(metComp))){
    isTransport <- FALSE
    isTranformation <- TRUE
    isMassBalanced <- TRUE

    # check for any biochemical transformation
    LHSm <- gsub("_[^_]*$|\\[.+\\]$","",met[Scoef < 0])
    RHSm <- gsub("_[^_]*$|\\[.+\\]$","",met[Scoef > 0])
    if(setequal(LHSm, RHSm))
      isTranformation <- FALSE


    # transport if: (i) transport of net mass (only works if reaction is mass
    # balanced) or (ii) transport of metabolite without transformation
    if(length(unique(metComp))>1) {

      # mass balanced?
      if(!any(is.na(metChemicalFormula))) {
        eleMU <- countElements(metChemicalFormula,
                               multiplier = Scoef)

        isMassBalanced <- all(abs(apply(eleMU,2,sum)) <= 1e-6) # some tolerance here.

        # Mass changes within compartments due to transport?
        compMassBalance <- lapply(unique(metComp),
                                  function(comp) {
                                    apply(eleMU[metComp == comp,, drop = FALSE],
                                          2,
                                          sum)
                                  })
        isTransport <- any(abs(unlist(compMassBalance)) > 1e-6) && isMassBalanced
      }

      # Transport of molecule not undergoing transformation
      if(!isTransport) {
        isTransport <- length(intersect(LHSm,RHSm)) > 0
      }
    }

    if(!isTransport && isTranformation)
      return("SBO:0000176")
    if(isTransport && isTranformation)
      return("SBO:0000655")
    if(isTransport && !isTranformation)
      return("SBO:0000185")
  }

  return("SBO:0000167")
}
