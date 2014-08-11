## Download trees

library (treebase)

data <- search_treebase(input = 100, by = 'ntax',
                max_trees = 1, branch_lengths = FALSE)

data <- search_treebase(input = 'mammal', by = 'taxon',
                        max_trees = 10, branch_lengths = FALSE,
                        only_metadata = TRUE)

nterms <- length(by)
search_term <- character(nterms)
section <- character(nterms)
for (i in 1:nterms) {
  search_term[i] <- switch(by[i], abstract = "dcterms.abtract", 
                           citation = "dcterms.bibliographicCitation", author = "dcterms.contributor", 
                           subject = "dcterms.subject", id.matrix = "tb.identifier.matrix", 
                           id.matrix.tb1 = "tb.identifer.matrix.tb1", ncbi = "tb.identifier.ncbi", 
                           id.study = "tb.identifier.study", id.study.tb1 = "tb.identifier.study.tb1", 
                           id.taxon = "tb.identifer.taxon", taxon.tb1 = "tb.identifier.taxon.tb1", 
                           taxonVariant.tb1 = "tb.identifier.taxonVarient.tb1", 
                           id.tree = "tb.identifier.tree", ubio = "tb.identifier.ubio", 
                           kind.tree = "tb.kind.tree", nchar = "tb.nchar.matrix", 
                           ntax = "tb.ntax.matrix", quality = "tb.quality.tree", 
                           matrix = "tb.title.matrix", study = "tb.title.study", 
                           taxon = "tb.title.taxon", taxonLabel = "tb.title.taxonLabel", 
                           taxonVariant = "tb.title.taxonVariant", tree = "tb.title.tree", 
                           type.matrix = "tb.type.matrix", type.tree = "tb.type.tree", 
                           doi = "prism.doi")
  section[i] <- switch(by[i], abstract = "study", citation = "study", 
                       author = "study", subject = "study", id.matrix = "matrix", 
                       id.matrix.tb1 = "matrix", ncbi = "taxon", id.study = "study", 
                       id.study.tb1 = "study", id.taxon = "taxon", taxon.tb1 = "taxon", 
                       taxonVariant.tb1 = "taxon", id.tree = "tree", ubio = "taxon", 
                       kind.tree = "tree", nchar = "matrix", ntax = "matrix", 
                       quality = "tree", matrix = "matrix", study = "study", 
                       taxon = "taxon", taxonLabel = "taxon", taxonVariant = "taxon", 
                       tree = "tree", type.matrix = "matrix", type.tree = "tree", 
                       doi = "study")
}
if (!all(section == section[1])) 
  stop("Multiple queries must belong to the same section (study/taxon/tree/matrix)")
search_type <- paste(section[1], "/find?query=", sep = "")
search_term[1] <- paste(search_term[1], "=", sep = "")
if (nterms > 1) {
  for (i in 2:nterms) {
    input <- sub("(and|or) ", paste("\\1%20", search_term[i], 
                                    "=", sep = ""), input)
  }
}
input <- gsub(" +", "%20\\1", input)
input <- gsub("\"", "%22", input)
input <- gsub("'", "%22", input)
if (by %in% c("doi")) 
  input <- paste("%22", input, "%22", sep = "")
if (exact_match) {
  search_term <- gsub("=", "==", search_term)
}
returns <- match.arg(returns)
schema <- switch(returns, tree = "tree", matrix = "matrix")
format <- "&format=rss1"
query <- paste("http://purl.org/phylo/treebase/phylows/", 
               search_type, search_term[1], input, format, "&recordSchema=", 
               schema, sep = "")
message(query)
if (max_trees == Inf) 
  max_trees <- "last()"
out <- get_nex(query, max_trees = max_trees, returns = returns, 
               curl = curl, pause1 = pause1, pause2 = pause2, attempts = attempts, 
               only_metadata = only_metadata)
if (schema == "tree" && only_metadata == FALSE) {
  out <- drop_nontrees(out)
  if (branch_lengths) {
    have <- have_branchlength(out)
    out <- out[have]
  }
}
out
}