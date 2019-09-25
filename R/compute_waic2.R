#'@title Compute type 2 waic penalty
#'
#'@description This function computes the type 1 penalty for waic (pwaic2)
#'  described in
#'  \href{https://link.springer.com/article/10.1007/s11222-013-9416-2}{Gelman et
#'  al. (2013) information criterion paper} for multi-scale occupancy models.
#'
#'@param msocc_mod output from \code{\link{msocc_mod}}
#'
#'@return numeric value that is the pwaic2
#'@export

compute_pwaic2 <- function(msocc_mod){
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

  tmp <- log(site * sample * rep)[2:num.mcmc,]
  pwaic2 <- sum(apply(tmp, 2, var))
  return(pwaic2)
}
