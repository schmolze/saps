#' Implements Significance Analysis of Prognostic Signatures (SAPS), a
#' robust method for determining prognostically significant gene sets.
#'
#' \code{\link{saps}} will usually be the only function needed.
#'
#'
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
#' @name saps-package
#' @docType package
NULL


#' @export
#' @title Compute SAPS statistics
#' @description This is the main user interface to the \pkg{saps} package, and is
#' usually the only function needed.
#' @param candidateGeneSets A matrix with at least one row, where each row represents
#' a gene set, and the column values are gene identifiers. The row names should contain
#' unique names for the gene sets. The column values may contain \code{NA} values, since
#' in general gene sets will have differing lengths.
#' @param dataSet A matrix, where the column names are gene identifiers
#' (in the same format as the values in \code{candidateGeneSets})
#' and the values are gene expression levels. Each row should contain data for a single
#' patient.
#' @param survivalTimes A vector of survival times. The length must equal the number of
#' rows (i.e. patients) in \code{dataSet}. This parameter is used by
#' \code{\link[survcomp]{concordance.index}} to compute the concordance index, which
#' is in turn used by \code{\link{calculatePEnrichment}} to compute P_enrichment.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1). The length must equal the number of rows
#' (i.e. patients) in \code{dataSet}.
#' is in turn used by \code{\link{calculatePEnrichment}} to compute P_enrichment.
#' @param random.samples An integer that specifies how many random gene sets to sample
#' when computing P_random. Used by \code{\link{calculatePRandom}}.
#' @param cpus An integer that specifies the number of cpus/cores to be used when
#' calculating P_enrichment. If greater than 1 (the default), the \pkg{snowfall}
#' package must be installed or an error will occur.
#' @param verbose A boolean indicating whether to display status messages during
#' computation. Defaults to \code{TRUE}.
#' @return The function returns a matrix with the following columns:
#'
#' \code{Size P_pure  P_random  P_enrichment  direction saps_score}
#'
#' Each row represents a gene set, and the row names contain the gene set names.
#' \code{Size} is the size of the gene set. \code{P_pure}, \code{P_random}, and
#' \code{P_enrichment} are the respective SAPS statistics for the gene set, while
#' \code{direction} is the direction of association for the enrichment score.
#' saps_score is -log10 max(P_pure, P_random, P_enrichment) multiplied by
#' \code{direction}.
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
saps <- function(candidateGeneSets, dataSet, survivalTimes,
                 followup, random.samples=25, cpus=1, verbose=TRUE) {

  if ((cpus > 1) & (!is.installed("snowfall")))
    stop("'snowfall' package not found (required for multiple CPU support)")

  candidateSetCount <- nrow(candidateGeneSets)

  # prepare results
  results <- list("rankedGenes"=NA, "geneset.count"=candidateSetCount,
                  "genesets"=list())

  # get concordance index

  if (verbose)
    message("Calculating concordance index...", appendLF=FALSE)

  ci <- rankConcordance(dataSet, survivalTimes, followup)

  rankedGenes <- ci[, -1]

  results$rankedGenes <- rankedGenes

  if (verbose)
    message("done.")


  setNames = rownames(candidateGeneSets)
  geneNames <- colnames(dataSet)

  # compute saps statistics for each candidate gene set
  for (i in 1:candidateSetCount) {

    candidateGeneSet <- candidateGeneSets[i,,drop=FALSE]

    candidateSetName <- row.names(candidateGeneSet)

    setName <- setNames[i]

    # get candidate genes
    candidateGenes <- as.character(candidateGeneSet[!is.na(candidateGeneSet)])
    commonGenes <- intersect(geneNames, candidateGenes)

    candidateSetSize = length(commonGenes)

    # prepare results for this geneset
    set_results <- list("size" = candidateSetSize,
                    "saps_unadjusted" = vector(mode="numeric", length=3),
                    "saps_adjusted" = vector(mode="numeric", length=3),
                    "random_p_pures" = NA,
                    "direction" = NA,
                    "saps_score" = NA,
                    "saps_score_adj" = NA)

    names(set_results$saps_unadjusted) <- c("p_pure", "p_random", "p_enrich")
    names(set_results$saps_adjusted) <- c("p_pure", "p_random", "p_enrich")

    if(candidateSetSize == 0) {

      warning(c("No gene data found for gene set ", candidateSetName, ", cannot compute SAPS."))

    }
    else {

      if (verbose)
        message(c("Using gene set ", candidateSetName, ", size = ", candidateSetSize))

      set_results$size <- candidateSetSize

      scaledData <- scale(dataSet[,is.element(geneNames, commonGenes)])

      if (verbose)
        message("Calculating P_pure...", appendLF=FALSE)

      p_pure <- calculatePPure(scaledData, survivalTimes, followup)
      set_results$saps_unadjusted["p_pure"] <- p_pure

      if (verbose)
        message("done.")

      if (verbose)
        message("Calculating P_random...", appendLF=FALSE)

      p_random <- calculatePRandom(dataSet, candidateSetSize, p_pure, survivalTimes,
                                   followup, random.samples)

      set_results$saps_unadjusted["p_random"] <- p_random

      if (verbose)
        message("done.")

      if (verbose)
        message("Calculating P_enrichment...", appendLF=FALSE)

      gsa_results <- calculatePEnrichment(rankedGenes, candidateGeneSet, cpus)

      p_enrich <- gsa_results$P_enrichment
      direction <- gsa_results$direction

      set_results$saps_unadjusted["p_enrich"] <- p_enrich
      set_results["direction"] <- direction

      # calculate adjusted p-values
      set_results$saps_adjusted["p_pure"] <- p.adjust(p_pure, method="BH",
                                                      n=candidateSetCount)
      set_results$saps_adjusted["p_random"] <- p.adjust(p_random, method="BH",
                                                      n=candidateSetCount)
      set_results$saps_adjusted["p_enrich"] <- p.adjust(p_enrich, method="BH",
                                                      n=candidateSetCount)

      # calculate saps score, adjusted and unadjusted
      p_max <- max(p_pure, p_random, p_enrich)
      saps_score <- -log10(p_max) * direction
      saps_score_adj <- p.adjust(saps_score, method="BH", n=candidateSetCount)

      set_results["saps_score"] <- saps_score
      set_results["saps_score_adj"] <- saps_score_adj

      if (verbose)
        message("done.")

    }

    results$genesets[[setName]] <- set_results

  }

  # compute SAPS_score
  #results["saps_score"] <- apply(results, 1, calculateSAPSScore)

  return(results)

}


#' @export
#' @title Compute P_enrichment
#' @description This function performs a pre-ranked gene set enrichment
#' analysis (GSEA) to evaluate the degree to which a candidate gene set is
#' overrepresented at the top or bottom extremes of a ranked list of
#' concordance indices. This function is normally called by
#' \code{\link{saps}}.
#' @param rankedGenes An \emph{nx1} matrix of concordance indices for \emph{n}
#' genes. Generally this will be the z-score returned by
#' \code{\link{rankConcordance}}. The row names should contain gene identifiers.
#' @param candidateGeneSet A \emph{1xp} matrix of \emph{p} gene identifiers.
#' The row name should contain a name for the gene set.
#' @param cpus This value is passed to the \code{\link[piano]{runGSA}} function
#' in the \pkg{piano} package. For multi-core CPUs, this value should be set to
#' the number of cores (which will significantly improve the computational time).
#' @return The function returns a matrix with the following columns:
#'
#' \code{P_enrichment  direction}
#'
#' \code{P_enrichment} is the enrichment score, while \code{direction} is
#' either 1 or -1 depending on the direction of association.
#' @seealso \code{\link{saps}} \code{\link[piano]{runGSA}}
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
#'
#' Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, et al. (2005)
#' Gene set enrichment analysis: a knowledge-based approach for interpreting
#' genome-wide expression profiles. Proc Natl Acad Sci U S A 102: 15545–15550.
calculatePEnrichment <- function(rankedGenes, candidateGeneSet, cpus) {

  # reshape candidate gene set into long form for piano input
  candidateSetLong <- (melt(as.matrix(candidateGeneSet), na.rm=TRUE))[, c("value", "Var1")]

  gsc <- loadGSC(candidateSetLong, type="data.frame")

  # run GSEA using piano package
  gsa <- runGSA(rankedGenes, geneSetStat="gsea", signifMethod="geneSampling", adjMethod="fdr",
                gsc=gsc, gsSizeLim=c(1, 250), nPerm=1000, ncpus=cpus, verbose=FALSE)

  # get p-values for up and down-regulated gene sets
  gsa_results <- cbind(names(gsa$gsc), gsa$pDistinctDirUp, gsa$pDistinctDirDn)

  # build direction vector
  direction <- ifelse(is.na(gsa_results[, 2]), -1, 1)
  direction <- cbind(direction)
  rownames(direction) <- gsa_results[, 1]

  # merge the up and down p-values
  neg_pos <- cbind(gsa_results[, 2])

  neg_na <- is.na(neg_pos[, 1])

  neg_pos[neg_na] <- (gsa_results[, 3])[neg_na]

  gsa_results <- cbind(gsa_results[, 1], neg_pos)

  rownames(gsa_results) <- gsa_results[, 1]

  gsea_p <- cbind("P_enrichment"=as.numeric(gsa_results[, 2]))
  rownames(gsea_p) <- rownames(gsa_results)

  # merge p-values and direction vector
  results <- merge(gsea_p, direction, by="row.names")

  # restore row names
  rownames(results) <- results$Row.names
  results <- results[, -1]

  return(results)

}


#' @export
#' @title Compute P_pure
#' @description This function stratifies patients into two groups via k-means
#' clustering (k=2) on an \emph{nxp} matrix consisting of \emph{n} patients
#' and \emph{p} genes in the candidate prognostic set. It is normally called
#' by \code{\link{saps}}.
#' @param geneData An \emph{nxp} matrix consisting of \emph{n} patients
#' and \emph{p} genes in the candidate prognostic set.
#' @param survivalTimes A vector of survival times. The length must equal
#' the number of rows \emph{n} in \code{geneData}.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1). The length must equal the number of rows
#' (i.e. patients) in \code{geneData}.
#' @return A log-rank p-value indicating the probability that the two groups
#' show no survival difference.
#' @seealso \code{\link{saps}}
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
calculatePPure <- function(geneData, survivalTimes, followup) {

  # cluster the patients on the candidate genes
  cluster <- kmeans(geneData, 2)$cluster

  # compute probability of no survival difference
  survtest <- survdiff(Surv(survivalTimes, followup) ~ cluster)

  return (1 - pchisq(survtest$chisq, 1))

}


#' @export
#' @title Compute P_random
#' @description This function randomly samples gene sets, and calculates
#' P_pure (via \code{\link{calculatePPure}}) for each one. P_random is the
#' proportion of randomly sampled gene sets achieving a P_pure at least as
#' significant as the provided \code{p_pure}. This function is normally called
#' by \code{\link{saps}}.
#' @param dataSet A matrix, where the column names are gene identifiers
#' and the values are gene expression levels. Each row should contain data for a
#' single patient.
#' @param sampleSize The desired size for the randomly sampled gene sets.
#' @param p_pure The candidate P_pure against which to compare the P_pure
#' values for the randomly generated gene sets.
#' @param survivalTimes A vector of survival times. The length must equal
#' the number of rows in \code{dataSet}.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1). The length must equal the number of rows
#' (i.e. patients) in \code{dataSet}.
#' @param random.samples The number of random gene sets to sample.
#' @return The proportion of randomly sampled gene sets with a calculated
#' P_pure at least as significant as the provided \code{p_pure}.
#' @seealso \code{\link{saps}}
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
calculatePRandom <- function(dataSet, sampleSize, p_pure, survivalTimes, followup, random.samples=25) {

  geneNames = colnames(dataSet)

  if (sampleSize == 1)
    sampleSize <- 2

  # vector to hold P_pure values
  p_pures <- vector(length=random.samples)

  for (i in 1:random.samples) {

    randomGeneNames <- sample(geneNames, sampleSize)

    randomGeneSet <- scale(dataSet[, randomGeneNames])

    p_pures[i] <- calculatePPure(randomGeneSet, survivalTimes, followup)

  }

  return(sum(p_pures <= p_pure)/random.samples)

}


#' @export
#' @title Compute concordance indices
#' @description Computes concordance indices for a gene expression data set,
#' and returns the concordance index and the z-score.
#' @details This function is a wrapper for
#' \code{\link[survcomp]{concordance.index}} in the \pkg{survcomp}
#' package. It applies the latter over the columns of \code{dataset} to
#' calculate concordance indices and the corresponding z-score for each gene.
#' @param dataset A matrix, where the column names are gene identifiers
#' and the values are gene expression levels. Each row should contain data for
#' a single patient.
#' @param survivalTimes A vector of survival times. The length must equal the number of
#' rows (i.e. patients) in \code{dataset}.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1). The length must equal the number of rows
#' (i.e. patients) in \code{dataSet}.
#' @return The function returns a matrix with two columns:
#'
#' \code{cindex  z}
#'
#' and as many rows as \code{dataset}. \code{cindex} and \code{z}
#' are the concordance index and corresponding z-score. The
#' row names contain the gene identifiers.
#' @seealso \code{\link{saps}} \code{\link[survcomp]{concordance.index}}
rankConcordance <- function(dataset, survivalTimes, followup) {

  concordance_f <- function(x) {

    tt <- concordance.index(x, surv.time=survivalTimes, surv.event=followup, method="noether", na.rm=TRUE)
    return (c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))

  }

  return (t(apply(dataset, 2, concordance_f)))

}

is.installed <- function(pkg) {
  return (is.element(pkg, installed.packages()[,1]))
}
