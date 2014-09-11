# functions to produce summarised reports of SAPS data

#' @export
summarizeSaps <- function(geneSets) {

  extractSaps <- function(geneSet) {

      name <- geneSet["name"]
      size <- geneSet["size"]
      dir <- geneSet["direction"]
      saps <- geneSet$saps_unadjusted
      saps_adj <- geneSet$saps_adjusted

      mat <- matrix(data=c(size, dir, saps["p_pure"], saps["p_random"],
                    saps["p_enrich"], saps_adj["p_pure"],
                    saps_adj["p_random"], saps_adj["p_enrich"]),
             nrow=1, ncol=8)

      return(mat)

  }

  return(sapply(geneSets, extractSaps))

}
