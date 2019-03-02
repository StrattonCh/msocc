#'@title Compute log pointwise predictive density
#'
#'@description This function computes the log pointwise predictive density (lppd)
#'  described in Gelman et al. (2013) for multi-scale occupancy models.
#'
#'@param MSOcc_mod output from \code{\link{MSOcc_mod}}
#'
#'@return numeric value that is the lppd

compute_lppd <- function(MSOcc_mod){
  #convert MCMC samples to probability
  psi.mcmc <- as.matrix(psi_mcmc(MSOcc_mod), ncol = M)
  theta.mcmc <- theta_mcmc(MSOcc_mod)
  p.mcmc <- MSOcc_mod$p

  #pull model info
  num.mcmc <- dim(psi.mcmc)[1]
  J <- MSOcc_mod$model.info$J
  K <- MSOcc_mod$model.info$K
  z.mcmc <- MSOcc_mod$model.info$z
  A.mcmc <- MSOcc_mod$model.info$A
  Y <- MSOcc_mod$model.info$Y

  p.mat <- matrix(rep(p.mcmc, sum(J)), nrow = num.mcmc, byrow = F)
  y.mat <- matrix(rep(Y, num.mcmc), nrow = num.mcmc, byrow = T)
  y.choose.mat <- matrix(rep(choose(K, Y), num.mcmc), nrow = num.mcmc, byrow = T)
  z.vec.mcmc <- matrix(0, num.mcmc, sum(J))
  psi.vec.mcmc <- matrix(0, num.mcmc, sum(J))
  for(i in 1:num.mcmc){
    z.vec.mcmc[i,] <- rep(z.mcmc[i,], J)
    psi.vec.mcmc[i,] <- rep(psi.mcmc[i,], J)
  }

  site <- psi.vec.mcmc ^ z.vec.mcmc * (1 - psi.vec.mcmc)^(1 - z.vec.mcmc)
  sample <- (z.vec.mcmc*theta.mcmc)^A.mcmc * (1 - z.vec.mcmc*theta.mcmc)^(1 - A.mcmc)
  rep <- y.choose.mat * (A.mcmc*p.mat)^y.mat * (1 - A.mcmc*p.mat)^(K - y.mat)

  tmp <- (site * sample * rep)[2:num.mcmc,]
  return(sum(log(colMeans(tmp))))
}
