#'@title Posterior samples of derived site-level occurence probability
#'
#'@description This function converts posterior samples from site-level
#'  regression coefficients into posterior samples from the derived site-level
#'  occurence probability using the logit link function. If psi was sampled
#'  directly using a beta-binomial sampler, if restructures these samples into
#'  an appropriate format.
#'
#'@param MSOcc_mod output from \code{\link{MSOcc_mod}}
#'
#'@return an object of class \code{matrix} of dimension \code{num.mcmc} by \code{M} of
#'  posterior samples of psi
#'
#'@example examples/psi_mcmc_ex.R

psi_mcmc <- function(MSOcc_mod){
  site_mod_char <- as.character(paste0(MSOcc_mod$model.info$site_mod, collapse = ""))
  site_hier <- grepl("[|]", site_mod_char)
  beta <- MSOcc_mod$beta
  X <- MSOcc_mod$model.info$X

  if(site_hier){
    psi <- matrix(0, dim(beta)[1], dim(beta)[3])
    for(i in 1:dim(beta)[3]){
      psi[,i] <- exp(beta[,,i] %*%  matrix(X[,,1], ncol = 1)) / (1 + exp(beta[,,i] %*%  matrix(X[,,1], ncol = 1)))
    }
    colnames(psi) <- dimnames(X)[[3]]
  } else if(!site_hier) {
    if(MSOcc_mod$model.info$site_mod == as.formula("~1") & MSOcc_mod$model.info$beta_bin){
      M <- MSOcc_mod$model.info$M
      psi <- matrix(rep(MSOcc_mod$psi, M), ncol = M)
      colnames(psi) <- rownames(X)
    } else {
      eta <- beta %*% t(X)
      psi <- as.matrix(exp(eta) / (1 + exp(eta)), ncol = M)
    }
  }

  return(psi)
}
