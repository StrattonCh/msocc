#'@title Posterior samples of derived site-level occurence probability
#'
#'@description This function converts posterior samples from site-level
#'  regression coefficients into posterior samples from the derived site-level
#'  occurence probability using the logit link function. If psi was sampled
#'  directly using a beta-binomial sampler, if restructures these samples into
#'  an appropriate format.
#'
#'@details This function returns one column for each replicate and arranges the
#'  columns by site, sample, and replicate.
#'
#'@param msocc_mod output from \code{\link{msocc_mod}}
#'
#'@return an object of class \code{matrix} of dimension \code{num.mcmc} by
#'  \code{sum(K)} of posterior samples of psi.
#'
#'@example examples/psi_mcmc_ex.R
#'@export

psi_mcmc <- function(msocc_mod){
  num.mcmc <- msocc_mod$model.info$num.mcmc
  beta <- msocc_mod$beta
  X <- msocc_mod$model.info$X
  M <- msocc_mod$model.info$M
  K <- msocc_mod$model.info$K

  if(msocc_mod$model.info$site_mod == as.formula("~1") & msocc_mod$model.info$beta_bin){
    psi <- matrix(rep(msocc_mod$psi, sum(K)), ncol = sum(K))
    colnames(psi) <- NULL
  } else {
    rep.X <- with(msocc_mod$model.info$df, tapply(rep, site, length))
    X.full <- matrix(0, sum(K), ncol(X))
    for(j in 1:ncol(X)){
      X.full[,j] <- rep(X[,j], rep.X)
    }

    eta <- beta %*% t(X.full)
    psi <- exp(eta) / (1 + exp(eta))
  }

  return(psi)
}
