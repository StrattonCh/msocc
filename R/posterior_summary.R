#'@title Summarize posterior samples
#'
#'@description This function allows for easy summaries of the posterior samples
#'  drawn by \code{\link{msocc_mod}}.
#'
#'@param msocc_mod output from \code{\link{msocc_mod}} (object of class \code{msocc})
#'@param burnin number of samples to discard as burnin when summarizing the
#'  posterior
#'@param print should summaries be printed in the console?
#'@param level one of \code{c('overall', 'site', 'sample', 'rep')}; which level of the
#'  design should be summarized?
#'@param quantiles vector of quantiles for credibility intervals
#'@param n number of rows to print in the console for the posterior summary
#'@param unique should only unique rows be printed?
#'
#'@return object of class \code{data.frame} providing summaries of the posterior
#'  organized by site, sample, and replicate.
#'
#'@example examples/posterior_summary_ex.R
#'@export

posterior_summary <- function(msocc_mod, burnin = 0, print = F, level = 'overall', quantiles = c(0.025, 0.975),
                              unique = T){
  # error check
  if(!(level %in% c("overall", "site", "sample", "rep"))){
    stop('level should be set to one of: overall, site, sample, rep.')
  }
  if(!('msocc' %in% class(msocc_mod))) stop('msocc_mod should be an object of class msocc.')

  #pull model info
  num.mcmc <- msocc_mod$model.info$num.mcmc
  M <- msocc_mod$model.info$M
  J <- msocc_mod$model.info$J
  K <- msocc_mod$model.info$K

  #convert MCMC samples to probability, accounting for burnin
  psi.mcmc <- psi_mcmc(msocc_mod)[(burnin+1):num.mcmc,]
  theta.mcmc <- theta_mcmc(msocc_mod)[(burnin+1):num.mcmc,]
  p.mcmc <- p_mcmc(msocc_mod)[(burnin+1):num.mcmc,]

  if(level == "overall"){
    #build summary table
    psi <- apply(psi.mcmc, 2, median)
    theta <- apply(theta.mcmc, 2, median)
    p <- apply(p.mcmc, 2, median)
    sum_tbl <- msocc_mod$model.info$df %>%
      dplyr::mutate(psi = psi, theta = theta, p = p)
    sum_tbl.unique <- sum_tbl %>% dplyr::distinct(site, psi, theta, p, .keep_all = TRUE)

    #print table
    if(print){
      cat(paste('Overall summary of occupancy given by posterior medians: \n'))
    }
  }
  if(level == "site"){
    #build table
    mean <- apply(psi.mcmc, 2, mean)
    median <- apply(psi.mcmc, 2, median)
    lwr <- apply(psi.mcmc, 2, quantile, probs = quantiles[1])
    upr <- apply(psi.mcmc, 2, quantile, probs = quantiles[2])

    sum_tbl <- msocc_mod$model.info$df %>%
      dplyr::mutate(median = median,
                    mean = mean,
                    lwr = lwr,
                    upr = upr) %>%
      dplyr::select(-sample, -rep) %>%
      dplyr::distinct()
    sum_tbl.unique <- sum_tbl %>% dplyr::distinct(site, median, mean, lwr, upr, .keep_all = TRUE)
    names(sum_tbl)[4:5] <- names(sum_tbl.unique)[4:5] <- c(noquote(as.character(quantiles)))


    #print table
    if(print){
      cat(paste('Posterior summary of occupancy at the site: \n'))
    }
  }
  if(level == "sample"){
    #build table
    mean <- apply(theta.mcmc, 2, mean)
    median <- apply(theta.mcmc, 2, median)
    lwr <- apply(theta.mcmc, 2, quantile, probs = quantiles[1])
    upr <- apply(theta.mcmc, 2, quantile, probs = quantiles[2])

    sum_tbl <- msocc_mod$model.info$df %>%
      dplyr::mutate(median = median,
                    mean = mean,
                    lwr = lwr,
                    upr = upr)
    sum_tbl.unique <- sum_tbl %>% dplyr::distinct(site, mean, median, lwr, upr, .keep_all = TRUE)
    names(sum_tbl)[6:7] <- names(sum_tbl.unique)[6:7] <- c(noquote(as.character(quantiles)))

    #print table
    if(print){
      cat(paste('Posterior summary of occupancy at the sample: \n'))
    }
  }
  if(level == "rep"){
    mean <- apply(p.mcmc, 2, mean)
    median <- apply(p.mcmc, 2, median)
    lwr <- apply(p.mcmc, 2, quantile, probs = quantiles[1])
    upr <- apply(p.mcmc, 2, quantile, probs = quantiles[2])

    sum_tbl <- msocc_mod$model.info$df %>%
      dplyr::mutate(median = median,
                    mean = mean,
                    lwr = lwr,
                    upr = upr)
    sum_tbl.unique <- sum_tbl %>% dplyr::distinct(site, mean, median, lwr, upr, .keep_all = TRUE)
    names(sum_tbl)[6:7] <- names(sum_tbl.unique)[6:7] <- c(noquote(as.character(quantiles)))

    #print table
    if(print){
      cat(paste('Posterior summary of detection in a replicate: \n'))
    }
  }

  if(unique){
    out <- sum_tbl.unique
  } else{
    out <- sum_tbl
  }

  return(out)
}
