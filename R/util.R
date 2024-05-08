### initialize the EB-NMF fit to covariance matrix YY' s.t. E[YY'] = LL'+ D from an estimate of L without nonnegative constraints
init_cov_ebnmf <- function(fl, kset = 1:ncol(fl$flash_fit$EF[[1]])) {
  LL <- fl$flash_fit$EF[[1]][, kset, drop = FALSE]
  FF <- fl$flash_fit$EF[[2]][, kset, drop = FALSE]
  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  LL <- cbind(fl$flash_fit$EF[[1]][, -kset, drop = FALSE], LL)
  FF <- cbind(pmax(FF, 0), pmax(-FF, 0))
  FF <- cbind(fl$flash_fit$EF[[2]][, -kset, drop = FALSE], FF)
  to.keep <- (colSums(LL) > .Machine$double.eps) & (colSums(FF) > .Machine$double.eps)
  LL <- LL[, to.keep, drop = FALSE]
  FF <- FF[, to.keep, drop = FALSE]
  return(list(LL, FF))
}


# Fit EBMF to covariance matrix YY' s.t. E[YY'] = LL'+ D, where D =
# s2*I and I is an identity matrix.
#
#' @importFrom Matrix rowSums
fit_ebmf_to_YY <- function (dat, fl, extrapolate = TRUE, warmstart = TRUE,
                            maxiter = 500, tol = NULL, epsilon = 2e-2,
                            verbose = 1) {
  if (is.matrix(dat) || inherits(dat, "Matrix")) {
    data_diag <- diag(dat)
  } else {
    data_diag <- Matrix::rowSums(dat$U * dat$D * dat$V)
  }
  
  s2 <- max(0, mean(data_diag - rowSums(fl$L_pm * fl$F_pm)))
  s2_diff <- Inf
  old_s2 <- 0
  
  ### alternate between estimating s2 and backfitting until convergence
  while(s2 > 0 && abs(s2_diff - 1) > epsilon) {
    if (is.matrix(dat) || inherits(dat, "Matrix")) {
      dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    }
    else{
      dat_minuss2 <- list(U = dat$U, D = dat$D, V = dat$V,
                          S = Matrix::Diagonal(nrow(dat$U), -s2))
    }
    Y2_diff <- sum((data_diag - s2)^2 - (data_diag - old_s2)^2)
    fl <- fl %>%
      flash_update_data(dat_minuss2, Y2_diff = Y2_diff) %>%
      flash_backfit(extrapolate = extrapolate, warmstart = warmstart, maxiter = maxiter, tol = tol, verbose = verbose)
    
    old_s2 <- s2
    s2 <- max(0, mean(data_diag - rowSums(fl$L_pm * fl$F_pm)))
    s2_diff <- s2 / old_s2
  }
  
  return(list(dat = dat, fl = fl, s2 = s2))
}

# Fit EB-SNMF to Y to estimate F, with L previously calculated from
# GBCD.
#
#' @importFrom stats cor
#' @importFrom stats lm
#' @importFrom stats aov
fit_ebmf_to_Y <- function(Y, fit.cov, corr_thres, maxiter){
  ### scale GEP membership estimates to 0-1 scale, and calculate Pearson correlations between L and L-tilde
  k.order <- order(fit.cov$pve, decreasing = TRUE)
  fit.L <- fit.cov$L_pm[, k.order]
  fit.L <- t(t(fit.L)/apply(fit.L, 2, max))
  corr <- diag(cor(fit.cov$L_pm[, k.order], fit.cov$F_pm[, k.order]))
  
  ### estimate GEP signatures by fitting EB-SNMF to gene expression data matrix with fixed L estimated from covariance decomposition above
  init.F <- t(solve(crossprod(fit.L), as.matrix(Matrix::crossprod(fit.L, Y))))
  fit.snmf <- flash_init(Y, S = 1/sqrt(nrow(Y)), var_type = 2) %>%
    flash_factors_init(
      init = list(as.matrix(fit.L), as.matrix(init.F)),
      ebnm_fn = c(ebnm::ebnm_generalized_binary, ebnm::ebnm_point_laplace)
    ) %>%
    flash_factors_fix(kset = 1:ncol(fit.L), which_dim = "loadings") %>%
    flash_backfit(extrapolate = FALSE, maxiter = maxiter, verbose = 1)
  
  ### calculate the z-score and lfsr of GEP signatures by running linear regression followed by ash
  J <- ncol(Y)
  genes <- colnames(Y)
  F.est <- matrix(0, J, ncol(fit.snmf$L_pm))
  rownames(F.est) <- genes
  colnames(F.est) <- colnames(fit.snmf$L_pm)
  F.se <- F.est
  F.z.ash <- F.est
  F.lfsr.ash <- F.est
  
  for (j in 1:J) {
    y     <- Y[,j]
    dat   <- as.data.frame(cbind(y, fit.snmf$L_pm))
    fit   <- lm(y ~ 0 + ., dat)
    coefs <- summary(fit)$coefficients
    F.est[j,] <- coefs[, "Estimate"]
    F.se[j,]  <- coefs[, "Std. Error"]
  }
  
  for(k in 1:ncol(F.est)){
    fit <- ashr::ash(F.est[,k], F.se[,k], mixcompdist = "normal", method = "shrink")
    F.z.ash[,k] <- fit$result$PosteriorMean/fit$result$PosteriorSD
    F.lfsr.ash[,k] <- fit$result$lfsr
  }
  
  ### return the estimated memberships and signatures only for GEPs whose L and Ltilde from covariance decomposition are highly concordant
  k.idx <- corr > corr_thres
  L.pm <- fit.snmf$L_pm[, k.idx]
  F.lfc <- fit.snmf$F_pm[, k.idx]/log(2)
  F.z <- F.z.ash[, k.idx]
  F.lfsr <- F.lfsr.ash[, k.idx]
  colnames(L.pm) <- c("Baseline", paste0("GEP", 1:(ncol(L.pm)-1)))
  colnames(F.lfc) <- colnames(L.pm)
  colnames(F.z) <- colnames(L.pm)
  colnames(F.lfsr) <- colnames(L.pm)
  return(list(L = L.pm, F = list(lfc = F.lfc, z_score = F.z, lfsr = F.lfsr)))
}
