#' @rdname gep_annotation
#'
#' @title Annotate GEPs 
#' 
#' @description Annotate the GEPs produced by GBCD using the patient
#'   identity information of cells
#' 
#' @param L cell x GEP matrix of GEP membership returned by
#'   \code{\link{fit_gbcd}}.
#' 
#' @param metadata data frame containing at least a column
#'   \dQuote{subject} indicating the patient identity of cells. Each
#'   row contains a cell and the rownames must match those of \code{L}.
#'  
#' @return A data frame containing GEPs in rows and the following two
#'   columns: the column \dQuote{patient_specific} reflects the degree
#'   of patient-specific expression of a GEP, and larger values indicate
#'   a GEP is more patient-specific rather than expressed in multiple
#'   patients. The column \dQuote{between_patient_variation} quantifies
#'   how strongly a GEP's memberships vary between patients as opposed
#'   to within patients, and larger values indicate a GEP has larger
#'   between- than within-patient variation.
#'  
#' @export
#' 
gep_annotation <- function (L, metadata) {
  ### check the input
  if(!all(rownames(L) %in% rownames(metadata)))
    stop("The rownames of L and metadata do not match!")
  if(!"subject" %in% colnames(metadata))
    stop("The column subject is missing from metadata!")
  GEPs <- colnames(L)
  L <- t(t(L)/apply(L, 2, max))
  colnames(L) <- GEPs
  metadata <- metadata[match(rownames(L), rownames(metadata)),]
  if(!is.factor(metadata$subject))
    subject <- as.factor(metadata$subject)
  
  ### quantify patient-specific vs shared expression
  res <- rep(0, ncol(L))
  for(k in 1:ncol(L)){
    table.k <- as.matrix(table(L[, k] > 0.05, subject))
    ratio.k <- table.k[2, ]/colSums(table.k)
    res[k] <- 1-max(ratio.k[ratio.k!=max(ratio.k)])/max(ratio.k)
  }
  
  ### quantify between vs within-patient variation
  res2 <- rep(NA, ncol(L))
  for(k in 1:ncol(L)){
    oneway.anova <- aov(L[, k] ~ subject)
    res2[k] <- summary(oneway.anova)[[1]][1,2]/(summary(oneway.anova)[[1]][1,2] + summary(oneway.anova)[[1]][2,2])
  }
  
  ### return the results
  result <- data.frame(patient_specific = res, between_patient_variation = res2)
  rownames(result) <- colnames(L)
  return(result)
}
