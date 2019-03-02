#'@title Posterior samples of derived sample-level occurence probability
#'
#'@description This function converts posterior samples from sample-level
#'  regression coefficients into posterior samples from the derived sample-level
#'  occurence probability using the logit link function. If theta was sampled
#'  directly using a beta-binomial sampler, if restructures these samples into
#'  an appropriate format.
#'
#'@param MSOcc_mod output from \code{\link{MSOcc_mod}}
#'
#'@return an object of class \code{matrix} of dimension \code{num.mcmc} by \code{sum(J)}
#'  of posterior samples of theta
#'
#'@example examples/theta_mcmc_ex.R

theta_mcmc <- function(MSOcc_mod){
  sample_mod_char <- as.character(paste0(MSOcc_mod$model.info$samp_mod, collapse = ""))
  sample_hier <- grepl("[|]", sample_mod_char)
  J <- MSOcc_mod$model.info$J

  alpha <- MSOcc_mod$alpha
  W <- MSOcc_mod$model.info$W
  num_mcmc <- MSOcc_mod$model.info$num.mcmc

  if(sample_hier){
    names <- list()
    theta <- list()
    for(i in 1:dim(alpha)[3]){
      theta[[i]] <- exp(alpha[,,i] %*%  t(W[[i]])) / (1 + exp(alpha[,,i] %*% t(W[[i]])))
      names[[i]] <- paste(names(W)[i], 1:nrow(W[[i]]), sep = "_")
    }
    out <- matrix(unlist(theta), nrow = num_mcmc)
    colnames(out) <- unlist(names)
  } else if(!sample_hier){
    if(MSOcc_mod$model.info$beta_bin & MSOcc_mod$model.info$samp_mod == as.formula("~1")){
      out <- matrix(rep(MSOcc_mod$theta, sum(J)), ncol = sum(J))
      colnames(out) <- rownames(MSOcc_mod$model.info$W)
    } else if(MSOcc_mod$model.info$beta_bin & MSOcc_mod$model.info$samp_mod == as.formula("~site")){
      out <- t(apply(MSOcc_mod$theta, 1, FUN = function(x) rep(x, J)))
      colnames(out) <- rownames(MSOcc_mod$model.info$W)
    } else{
      nu <- alpha %*% t(W)
      out <- exp(nu) / (1 + exp(nu))
    }
  }
  return(out)
}
