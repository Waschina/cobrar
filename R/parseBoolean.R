parseBoolean <- function(gprRule, tokens = "()&|~") {

  # quit, if there is no gene association
  if ( is.na(gprRule) || (gprRule == "") ) {
      return(list(gene = character(0L), rule = ""))
  }

  if( grepl("\\s*\\(\\s*\\)\\s*", gprRule) ){
  	warning("found empty expression rule: '( )'. check if this is intended.")
  	return(list(gene = "", rule = ""))
  }

  gpr <- gsub("and ", "& ", gprRule, ignore.case = TRUE)
  gpr <- gsub("or ",  "| ", gpr, ignore.case = TRUE)
  gpr <- gsub("not ", "~ ", gpr, ignore.case = TRUE)
  gpr <- gsub("[", "", gpr, fixed = TRUE)
  gpr <- gsub("]", "", gpr, fixed = TRUE)
  rule <- gpr

  # split the rule into the gene names
  genes_tmp <- strsplit(gpr, paste("[", tokens, "]", sep = ""))

  # remove trailing and leading whitespaces
  genes_tmp <- gsub("(^\\s+)|(\\s+$)", "", genes_tmp[[1]], perl = TRUE)

  # remove empty entries in genes_tmp
  not_empty <- which(nchar(genes_tmp) > 0 )
  genes     <- genes_tmp[not_empty]

  # number of entries
  num_genes <- length(genes)


  # vector with unique gene numbers like "x[1]", "x[2]", "x[1]", ...

  # a unique vector with all genes
  gene_uniq <- unique(genes)

  newTok    <- match(genes, gene_uniq)
  newTok    <- sapply(newTok, function(x) paste("x[", x, "]", sep = ""))

  # replace gene names in rule by their newTok string (which is x[gene_number]):
  for (i in 1:num_genes) {
      rule <- sub(genes[i], newTok[i], rule, fixed = TRUE)
  }

  # return vector with unique gene names and the rule where numbers correspond to unique gene names
  return(list(gene = gene_uniq, rule = rule))
}
