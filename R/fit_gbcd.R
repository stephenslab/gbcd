#' @rdname fit_gbcd
#'
#' @title Fit Generalized Binary Covariance Decomposition
#' 
#' @description Fit a generalized binary covariance decomposition (GBCD)
#' to single cell RNA-seq data containing multiple tumors.
#' 
#' @param Y Cell x gene matrix of normalized and log-transformed gene
#'   expression data.
#' 
#' @param Kmax Integer 2 or greater used to determine the number of
#'   final number of factors in the GBCD. The final number of factors is
#'   often close to \code{Kmax}, and is never greater than \code{2*Kmax
#'   - 1}.
#'  
#' @param prior Nonnegative prior for GEP memberships, usually the
#'   generalized binary prior. The ratio \eqn{\rho = \sigma_k/\mu_k} in
#'   the generalized binary prior can be controlled via the
#'   \code{\link[flashier]{flash_ebnm}} function from flashier,
#'   \code{flash_ebnm(prior_family = "generalized_binary", scale =
#'   rho)}. See also the \code{\link[ebnm]{ebnm_generalized_binary}}
#'   function from the ebnm package for more information.
#'   
#' @param maxiter1 Positive integer specifying the maximum number of
#'   backfit iterations during the GEP membership matrix L
#'   initialization.
#'  
#' @param maxiter2 Positive integer specifying the maximum number of
#' backfit iterations during the GEP membership matrix L estimation.
#'   
#' @param maxiter3 Positive integer specifying the maximum number of
#' backfit iterations during the GEP signature matrix F estimation.
#'   
#' @param control List of control parameters with the following
#' elements: \dQuote{warmstart}, a logical indicator specifying
#' whether to use warmstart to initialize the prior g when solving
#' EBNM subproblems (see ebnm for details); \dQuote{extrapolate}, a
#' logical indicator specifying whether to use extrapolation to
#' accelerate backfitting GEP memberships (see flashier for details);
#' \dQuote{corr_thres}, a numeric scalar between 0 and 1 such that we
#' only keep columns \eqn{l_k} whose Pearson correlation with
#' \eqn{\tilde{l}_k} exceeds \dQuote{corr_thres}.
#' 
#' @param verbose Integer specifying whether and how to display
#'   progress updates, as described in flashier.
#'   
#' @return A list including the following elements:
#' 
#' \item{L}{cell x GEP matrix containing the posterior estimates of
#' the GEP membership matrix L.}
#' 
#' \item{F}{List containing the posterior summaries of the GEP
#'   signature matrix F.}
#'
#' @examples
#' # Please see the vignettes for examples.
#'
#' @import Matrix
#' @importFrom Matrix tcrossprod
#' @importFrom utils modifyList
#' @importFrom ebnm ebnm_point_laplace
#' @importFrom ebnm ebnm_generalized_binary
#' @importFrom flashier flash_init
#' @importFrom flashier flash_greedy
#' @importFrom flashier flash_backfit
#' @importFrom flashier flash_factors_init
#' @importFrom flashier flash_factors_remove
#' 
#' @export
#' 
fit_gbcd <- function (Y, Kmax, prior = ebnm::ebnm_generalized_binary, 
                      maxiter1 = 500, maxiter2 = 200, maxiter3 = 500,
                      control = list(), verbose = 1) {

  control <- modifyList(fit_gbcd_control_default(), control, keep.null = TRUE)
  extrapolate <- control$extrapolate
  warmstart <- control$warmstart
  corr_thres <- control$corr_thres
    
  ### form the covariance matrix from the cell by gene matrix of gene expression data
  print("Form cell by cell covariance matrix...")
  start_time = proc.time()
  if (2 * ncol(Y) * mean(Y > 0) < nrow(Y)) {
    # Use lowrank plus sparse representation:
    dat <- list(U = Y, D = rep(1 / ncol(Y), ncol(Y)), V = Y)
  } else {
    # Form covariance matrix:
    dat <- Matrix::tcrossprod(Y) / ncol(Y)
  }
  fit.init <- flash_init(dat, var_type = 0)
  runtime = proc.time() - start_time
  print(runtime)

  ### fit EBMF with point laplace prior to covariance matrix without
  ### considering the diagonal component
  print("Initialize GEP membership matrix L...")
  start_time = proc.time()
  fit.cov <- fit.init |>
    flash_greedy(Kmax = 1, ebnm_fn = ebnm_point_laplace) |>
    flash_greedy(Kmax = Kmax - 1, ebnm_fn = ebnm_point_laplace) |>
    flash_backfit(maxiter = 25, verbose = verbose)

  ### fit EBMF with point laplace prior to covariance matrix with the
  ### diagonal component
  fit.cov <- fit_ebmf_to_YY(dat = dat, fl = fit.cov, maxiter = maxiter1, verbose = verbose)$fl
  runtime = proc.time() - start_time
  print(runtime)

  ### initialize EB-NMF fit from the EBMF fit with point laplace prior
  print("Estimate GEP membership matrix L...")
  start_time = proc.time()
  cov.init <- init_cov_ebnmf(fit.cov)
  kmax <- which.max(colSums(cov.init[[1]]))
  fit.cov <- fit.init |>
    flash_factors_init(
      init = lapply(cov.init, function(x) x[, kmax, drop = FALSE]),
      ebnm_fn = ebnm_point_laplace
    ) |>
    flash_factors_init(
      init = lapply(cov.init, function(x) x[, -c(kmax), drop = FALSE]),
      ebnm_fn = prior
    ) |>
    flash_backfit(extrapolate = FALSE, warmstart = warmstart, maxiter = 25,
                  verbose = verbose)

  ### Remove factors with PVE of zero then refit EB-NMF to covariance
  ### matrix.
  kset <- (fit.cov$pve > 0)
  kall <- 1:fit.cov$n_factors
  if(!all(kset))
    fit.cov <- flash_factors_remove(fit.cov, kset=kall[!kset])
  fit.cov <- fit_ebmf_to_YY(dat = dat, fl = fit.cov,
                            extrapolate = extrapolate, warmstart = warmstart,
                            maxiter = maxiter2, verbose = verbose)$fl
  runtime = proc.time() - start_time
  print(runtime)
  
  ### estimate GEP signatures by fitting EB-SNMF to gene expression
  ### data matrix with fixed L estimated from covariance decomposition
  ### above
  print("Estimate GEP signature matrix F...")
  start_time = proc.time()
  res <- fit_ebmf_to_Y(Y, fit.cov, corr_thres, maxiter3)
  runtime = proc.time() - start_time
  print(runtime)

  return(res)
}

#' @rdname fit_gbcd
#'
#' @export
#' 
fit_gbcd_control_default <- function()
  list(warmstart = TRUE,
       extrapolate = FALSE,
       corr_thres = 0.8)
