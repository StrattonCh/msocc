#'@title Compute type 1 waic penalty
#'
#'@description This function computes the type 1 penalty for waic (pwaic1)
#'  described in Gelman et al. (2013) for multi-scale occupancy models.
#'
#'@param msocc_mod output from \code{\link{msocc_mod}}
#'
#'@return numeric value that is the pwaic1
#'@export

compute_pwaic1 <- function(msocc_mod){
  #first piece of the difference
  #convert MCMC samples to probability
  psi <- psi_mcmc(msocc_mod)
  theta <- theta_mcmc(msocc_mod)
  p <- p_mcmc(msocc_mod)

  #pull model info
  num.mcmc <- msocc_mod$model.info$num.mcmc
  J <- msocc_mod$model.info$J
  K <- msocc_mod$model.info$K
  z <- msocc_mod$model.info$z.vec[, rep(1:ncol(msocc_mod$model.info$z.vec), K)]
  A <- msocc_mod$model.info$A[, rep(1:ncol(msocc_mod$model.info$A), K)]
  y <- matrix(rep(msocc_mod$model.info$y, num.mcmc), nrow = num.mcmc, byrow = T)


  site <- psi^z * (1 - psi)^(1 - z)
  sample <- (z*theta)^A * (1 - z*theta)^(1 - A)
  rep <- (A*p)^y * (1 - A*p)^(1 - y)

  tmp <- (site * sample * rep)[2:num.mcmc,]
  tmp1 <- log(colMeans(tmp))

  #second piece of the difference
  tmp2 <- colMeans(log(tmp))

  #calculate measure
  pwaic1 <- 2 * sum(tmp1 - tmp2)
  return(pwaic1)
}
