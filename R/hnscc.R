#' @name hnscc
#'
#' @title Single-cell RNA-seq from HNSCC Tumor Cells
#'
#' @docType data
#'  
#' @description This dataset contains transcriptome profiles of
#' malignant cells generated via single-cell RNA sequencing. The cells
#' were collected from primary tumors in 10 HNSCC patients and
#' matching lymph node (LN) metastases from 5 of these patients.
#' Puram \emph{et al} (2017) found that each of these 10 patients
#' clearly mapped to a molecular subtype of HNSCC, whose signatures
#' were previously defined by analysis of bulk expression data of 279
#' TCGA HNSCC tumors.
#'
#' The data are normalized, log-transformed counts for 17,113 genes in
#' 2,176 cells: \eqn{y_{ij} = \log_2(1 + TPM_{ij}/10)}, where
#' \eqn{TPM_{ij}} is the transcript-per-million (TPM) value for gene
#' \eqn{j} in cell \eqn{i}. The counts are stored as a 2,176 x 17,113
#' sparse matrix.
#'
#' These data are used in the vignette to illustrate how gbcd
#' can be used to analyze single-cell RNA-seq data derived from
#' multiple tumor samples.
#'
#' @format \code{hnscc} is a list with the following elements:
#' 
#' \describe{
#'
#' \item{Y}{2,176 x 17,113 sparse matrix of normalized,
#'   log-transformed counts, with rows corresponding to cells and
#'   columns corresponding to genes.}
#'
#' \item{info}{Data frame containing information about the cells,
#'   including patient identity, primary/metastatic status and
#'   molecular subtype of the tumor sample.}
#'
#' \item{sample_col}{Color used in plots to indicate the different
#'   tumor samples.}
#' 
#' \item{subtype_col}{Color used in plots to indicate the different
#'   molecular subtypes.}}
#' 
#' @references
#' 
#' S.V. Puram \emph{et al} (2017). Single-cell transcriptomic analysis of 
#'   primary and metastatic tumor ecosystems in head and neck cancer
#'   \emph{Cell} \bold{171}, 1611-1624. \doi{10.1016/j.cell.2017.10.044}
#' 
#' @keywords data
#'
#' @examples
#' library(Matrix)
#' data(hnscc)
#' 
NULL
