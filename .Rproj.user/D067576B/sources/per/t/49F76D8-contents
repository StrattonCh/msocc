#'@title Summarize posterior samples
#'
#'@description This function allows for easy summaries of the posterior samples
#'  drawn by \code{\link{MSOcc_mod}}.
#'
#'@param MSOcc_mod output from \code{\link{MSOcc_mod}}
#'@param burnin number of samples to discard as burnin when summarizing the
#'  posterior
#'@param print should summaries be printed in the console?
#'@param level one of \code{c('overall', 'site', 'sample', 'rep')}; which level of the
#'  design should be summarized?
#'@param quantiles vector of quantiles for credibility intervals
#'@param n number of rows to print in the console for the posterior summary
#'
#'@return object of class \code{tbl} providing summaries of the posterior
#'
#'@example examples/posterior_summary_ex.R

posterior_summary <- function(MSOcc_mod, burnin = 0, print = T, level = 'overall', quantiles = c(0.025, 0.975),
                              n = 'all'){
  old <- options(tibble.print_max = 20, tibble.print_min = 10)

  options(tibble.print_max = n, tibble.print_min = 10)
  #pull model info
  num.mcmc <- MSOcc_mod$model.info$num.mcmc
  M <- MSOcc_mod$model.info$M
  J <- MSOcc_mod$model.info$J
  K <- MSOcc_mod$model.info$K

  #convert MCMC samples to probability, accounting for burnin
  psi.mcmc <- as.matrix(psi_mcmc(MSOcc_mod)[(burnin+1):num.mcmc,], ncol = M)
  theta.mcmc <- theta_mcmc(MSOcc_mod)[(burnin+1):num.mcmc,]
  p.mcmc <- MSOcc_mod$p[(burnin+1):num.mcmc]

  if(level == "overall"){
    #build summary table
    site_sample <- colnames(theta.mcmc)
    psi <- rep(apply(psi.mcmc, 2, median), J)
    theta <- apply(theta.mcmc, 2, median)
    p <- rep(median(p.mcmc), sum(J))
    sum_tbl <- tibble::tibble(site_sample, psi, theta, p)

    #print table
    if(print){
      cat(paste('Overall summary of occupancy given by posterior medians: \n'))
    }
    return(sum_tbl)
  }
  if(level == "site"){
    #build table
    site_sample <- colnames(theta.mcmc)
    mean <- rep(apply(psi.mcmc, 2, mean), J)
    median <- rep(apply(psi.mcmc, 2, median), J)
    lwr <- rep(apply(psi.mcmc, 2, quantile, probs = quantiles[1]), J)
    upr <- rep(apply(psi.mcmc, 2, quantile, probs = quantiles[2]), J)
    sum_tbl <- tibble::tibble(site_sample, median, mean, lwr, upr)
    names(sum_tbl)[4:5] <- c(noquote(as.character(quantiles)))

    #print table
    if(print){
      cat(paste('Posterior summary of occupancy at the site: \n'))
    }
    return(sum_tbl)
  }
  if(level == "sample"){
    #build table
    site_sample <- colnames(theta.mcmc)
    mean <- apply(theta.mcmc, 2, mean)
    median <- apply(theta.mcmc, 2, median)
    lwr <- apply(theta.mcmc, 2, quantile, probs = quantiles[1])
    upr <- apply(theta.mcmc, 2, quantile, probs = quantiles[2])
    sum_tbl <- tibble::tibble(site_sample, median, mean, lwr, upr)
    names(sum_tbl)[4:5] <- c(noquote(as.character(quantiles)))

    #print table
    if(print){
      cat(paste('Posterior summary of occupancy at the sample: \n'))
    }
    return(sum_tbl)
  }
  if(level == "rep"){
    p.mcmc <- matrix(p.mcmc, ncol = 1) #currently assumes fixed p
    #build table
    site_sample <- colnames(theta.mcmc)
    mean <- rep(apply(p.mcmc, 2, mean), sum(J))
    median <- rep(apply(p.mcmc, 2, median), sum(J))
    lwr <- rep(apply(p.mcmc, 2, quantile, probs = quantiles[1]), sum(J))
    upr <- rep(apply(p.mcmc, 2, quantile, probs = quantiles[2]), sum(J))
    sum_tbl <- tibble::tibble(site_sample, median, mean, lwr, upr)
    names(sum_tbl)[4:5] <- c(noquote(as.character(quantiles)))

    #print table
    if(print){
      cat(paste('Posterior summary of detection in a replicate: \n'))
    }
    return(sum_tbl)
  }
  if(!(level %in% c("overall", "site", "sample", "rep"))){
    stop('level should be set to one of: overall, site, sample, rep.')
  }
  on.exit(options(old), add = TRUE)
}
