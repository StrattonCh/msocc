#'@title Posterior samples of derived sample-level occurence probability
#'
#'@description This function converts posterior samples from sample-level
#'  regression coefficients into posterior samples from the derived sample-level
#'  occurence probability using the logit link function. If theta was sampled
#'  directly using a beta-binomial sampler, if restructures these samples into
#'  an appropriate format.
#'
#'@details This function returns one column for each replicate and arranges the
#'  columns by site, sample, and replicate.
#'
#'@param msocc_mod output from \code{\link{msocc_mod}}
#'
#'@return an object of class \code{matrix} of dimension \code{num.mcmc} by \code{sum(K)}
#'  of posterior samples of theta
#'
#'@example examples/theta_mcmc_ex.R
#'@export

theta_mcmc <- function(msocc_mod){
  num_mcmc <- msocc_mod$model.info$num.mcmc
  J <- msocc_mod$model.info$J
  K <- msocc_mod$model.info$K
  alpha <- msocc_mod$alpha
  W <- msocc_mod$model.info$W

  if(msocc_mod$model.info$beta_bin & msocc_mod$model.info$samp_mod == as.formula("~1")){
    out <- matrix(rep(msocc_mod$theta, sum(K)), ncol = sum(K))
    colnames(out) <- NULL
  } else if(msocc_mod$model.info$beta_bin & msocc_mod$model.info$samp_mod == as.formula("~site")){
    rep.theta <- with(msocc_mod$model.info$df, tapply(rep, site, length))
    out <- msocc_mod$theta[,rep(1:ncol(msocc_mod$theta), rep.theta)]
    colnames(out) <- NULL
  } else{
    W.full <- matrix(0, sum(K), ncol(W))
    for(j in 1:ncol(W)){
      W.full[,j] <- rep(W[,j], K)
    }

    nu <- alpha %*% t(W.full)
    out <- exp(nu) / (1 + exp(nu))
  }
  return(out)
}

