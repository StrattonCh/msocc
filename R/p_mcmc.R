#'@title Posterior samples of derived rep-level occurence probability
#'
#'@description This function converts posterior samples from rep-level
#'  regression coefficients into posterior samples from the derived rep-level
#'  occurence probability using the logit link function. If p was sampled
#'  directly using a beta-binomial sampler, if restructures these samples into
#'  an appropriate format.
#'
#'@details This function returns one column for each replicate and arranges the
#'columns by site, sample, and replicate.
#'
#'@param msocc_mod output from \code{\link{msocc_mod}}
#'
#'@return an object of class \code{matrix} of dimension \code{num.mcmc} by \code{sum(K)}
#'  of posterior samples of p
#'
#'@example examples/p_mcmc_ex.R
#'@export

p_mcmc <- function(msocc_mod){
  num.mcmc <- msocc_mod$model.info$num.mcmc
  M <- msocc_mod$model.info$M
  J <- msocc_mod$model.info$J
  K <- msocc_mod$model.info$K

  if('p' %in% names(msocc_mod)){
    out <- matrix(rep(msocc_mod$p, each = sum(K)), nrow = num.mcmc, ncol = sum(K), byrow = T)
  } else if(dim(msocc_mod$model.info$V)[1] == sum(J)){
    V <- msocc_mod$model.info$V
    V.exp <- matrix(0, sum(K), ncol(V))
    for(j in 1:ncol(V)){
      V.exp[,j] <- rep(V[,j], K)
    }
    nu_rep <- msocc_mod$delta %*% t(V.exp)
    out <- exp(nu_rep) / (1 + exp(nu_rep))
    colnames(out) <- NULL
  } else {
    nu_rep <- msocc_mod$delta %*% t(msocc_mod$model.info$V)
    out <- exp(nu_rep) / (1 + exp(nu_rep))
    colnames(out) <- NULL
  }

  return(out)
}

