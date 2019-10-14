#'@title Fit multi-scale occupancy model.
#'
#'@description This function fits the Bayesian multi-scale occupancy model
#'  described by \href{https://pubs.er.usgs.gov/publication/70194000}{Dorazio
#'  and Erickson (2017)} using the polya-gamma data augmentation strategy
#'  described by \href{https://arxiv.org/pdf/1205.0310.pdf}{Polson et al.
#'  (2012)}. Note that this documentation assumes there are \eqn{M} sites,
#'  \eqn{J_i} samples within each site, and \eqn{K_{ij}} replicates from each
#'  sample.
#'
#'@details This function fits the multi-scale occupancy model described by
#'  \href{https://pubs.er.usgs.gov/publication/70194000}{Dorazio and Erickson
#'  (2017)}. However, this function implements a fully Bayesian sampler based on
#'  the data augmentation strategy described by
#'  \href{https://arxiv.org/pdf/1205.0310.pdf}{Polson et al. (2012)}
#'
#'@param wide_data object of class \code{data.frame} containing site, sample,
#'  and PCR replicates in wide format. Column names should be \code{site},
#'  \code{sample}, \code{PCR1}, \code{PCR2}, ... and contain no other columns.
#'
#'@param site object of class \code{list} containing the following elements: \cr
#'  * \code{model} formula describing the within site model. \cr
#'  * \code{cov_tbl} object of class \code{data.frame} containing site specific
#'  covariates; \code{cov_tbl} should have exactly one row for each site and have a column named 'site'.
#'
#'@param sample object of class \code{list} containing the following elements: \cr
#'  * \code{model} formula describing the within sample model. \cr
#'  * \code{cov_tbl} object of class \code{data.frame} containing sample specific
#'  covariates; \code{cov_tbl} should have exactly one row for each sample and
#'  have columns named 'site' and 'sample'.
#'
#'@param rep object of class \code{list} containing the following elements: \cr
#'  * \code{model} formula describing the within replicate model. \cr
#'  * \code{cov_tbl} object of class \code{data.frame} containing replicate specific
#'  covariates; \code{cov_tbl} should have exactly one row for each sample if
#'  replicate specific covariates are aggregated at the sample level and contain
#'  columns named 'site' and 'sample'. Otherwise, \code{cov_tbl} should have
#'  exactly one row for each replicate and contain columns name 'site',
#'  'sample', and 'rep'.
#'
#'@param priors object of class \code{list} containing the following elements: \cr
#'  * \code{site} \verb{ } object of class \code{list} containing the
#'  following elements: \cr
#'    * \code{mu0} prior mean for site-level regression coefficients \cr
#'    * \code{Sigma0} prior covariance matrix for site-level regression
#'    coefficients
#'  * \code{sample} \verb{ } object of class \code{list} containing the
#'  following elements: \cr
#'    * \code{mu0} prior mean for sample-level regression coefficients \cr
#'    * \code{Sigma0} prior covariance matrix for sample-level regression
#'    coefficients
#'  * \code{rep} \verb{ } object of class \code{list} containing the
#'  following elements: \cr
#'    * \code{mu0} prior mean for replicate-level regression coefficients \cr
#'    * \code{Sigma0} prior covariance matrix for replicate-level regression
#'    coefficients
#'  * \code{a0,b0} numeric shape parameters for beta prior on probability
#'  parameters if sampled directly (see \code{beta_bin} argument)
#'@param num.mcmc number of MCMC samples
#'@param progress should sampling progress be printed?
#'@param print interval for printing; defaults to 5 percent of \code{num.mcmc}
#'  of left \code{NULL}
#'@param seed optional seed for reproducible samples
#'@param beta_bin optional; should a beta-binomial sampler be used when
#'  possible? This option is considerably faster.
#'
#'@example examples/msocc_mod_ex.R
#'
#'@return object of class \code{list} containing the following elements: \cr
#'  * \code{beta} an object of class \code{matrix} of samples from the joint
#'  posterior distribution of the regression coefficients at the site level
#'  * \code{psi} an object of class \code{numeric} of samples from the joint
#'  posterior distribution of the site-level presence probability, psi. Note
#'  that this is only returned if \code{beta_bin} is \code{TRUE} and
#'  \code{site$model = ~ 1}
#'  * \code{alpha} an object of class \code{matrix} of
#'  samples from the joint posterior distribution of the regression coefficients
#'  at the sample level
#'  * \code{theta} an object of class \code{matrix} of
#'  samples from the joint posterior distribution of the sample-level presence
#'  probability, theta. Note that this is only returned if \code{beta_bin} is
#'  \code{TRUE} and \code{sample$model = ~ 1} or \code{sample$model = ~ site}
#'  * \code{delta} an object of class \code{matrix} of samples from the joint
#'  posterior distribution of the regression coefficients at the replicate level
#'  * \code{p} an object of class \code{numeric} of sample from the posterior
#'  distribution of the replicate level detection probability. NOte that this is
#'  only returned if \code{beta_bin} is \code{TRUE} and \code{rep$model = ~ 1}.
#'  * \code{model.info} an object of class list containing the following
#'  elements: \cr
#'    * \code{X} design matrix for site level predictors
#'    * \code{W} design matrix for sample level predictors
#'    * \code{V} design matrix for replicate level predictors
#'    * \code{M} number of sites
#'    * \code{J} vector of number of samples per site
#'    * \code{K} vector of number of replicates per site-sample combination
#'    * \code{z} matrix of posterior samples of latent site-presence z
#'    * \code{z.vec} matrix of posterior samples of latent site-presence
#'    stretched across samples
#'    * \code{A} matrix of posterior samples of latent sample-presence A
#'    * \code{Y} vector of binomial responses (aggregated at sample level)
#'    * \code{y} vector of Bernoulli responses (stretched across site-sample
#'    combination)
#'    * \code{site_mod} model statement for site-level predictors
#'    * \code{samp_mod} model statement for sample-level predictors
#'    * \code{rep_mod} model statement for replicate-level predictors
#'    * \code{num.mcmc} number of MCMC samples run
#'    * \code{beta_bin} was beta-binomial sampler implemented if possible?
#'    * \code{df} empty \code{data.frame} of design in long format
#'
#'@importFrom magrittr %>%
#'@export
#'
#'@md

msocc_mod <- function(wide_data,
                      site = list(model = ~ 1, cov_tbl),
                      sample = list(model = ~ 1, cov_tbl),
                      rep = list(model = ~ 1, cov_tbl),
                      priors = list(site = list(mu0 = 0, Sigma0 = 9),
                                    sample = list(mu0 = 0, Sigma0 = 9),
                                    rep = list(mu0 = 0, Sigma0 = 9), a0 = 1, b0 = 1),
                      num.mcmc = 1000, progress = T, print = NULL, seed = NULL, beta_bin = T){

  # set optional seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # fix column names and arrange tables
  if(names(wide_data)[1] != 'site' | names(wide_data)[2] != 'sample') stop("First two columns of wide_data are not 'site' and 'sample' respectively. Please rename columns.")
  wide_data <- dplyr::arrange(wide_data, site, sample)

  if(!('site' %in% names(site$cov_tbl))) stop("There should be a column named 'site' in site$cov_tbl.")
  site$cov_tbl <- dplyr::arrange(site$cov_tbl, site)

  if(!('site' %in% names(sample$cov_tbl))) stop("There should be a column named 'site' in sample$cov_tbl.")
  if(!('sample' %in% names(sample$cov_tbl))) stop("There should be a column named 'sample' in sample$cov_tbl.")
  sample$cov_tbl <- dplyr::arrange(sample$cov_tbl, site, sample)

  if(nrow(rep$cov_tbl) == nrow(sample$cov_tbl)){
    if(!('site' %in% names(rep$cov_tbl))) stop("There should be a column named 'site' in rep$cov_tbl.")
    if(!('sample' %in% names(rep$cov_tbl))) stop("There should be a column named 'sample' in rep$cov_tbl.")
    rep$cov_tbl <- dplyr::arrange(rep$cov_tbl, site, sample)
  } else {
    if(!('site' %in% names(rep$cov_tbl))) stop("There should be a column named 'site' in rep$cov_tbl.")
    if(!('sample' %in% names(rep$cov_tbl))) stop("There should be a column named 'sample' in rep$cov_tbl.")
    if(!('rep' %in% names(rep$cov_tbl))) stop("There should be a column named 'rep' in rep$cov_tbl.")
    rep$cov_tbl <- dplyr::arrange(rep$cov_tbl, site, sample, rep)
  }

  # remove samples with all NAs
  Y.mat <- as.matrix(wide_data[,-c(1:2)])
  if(any(apply(Y.mat, 1, FUN = function(x) all(is.na(x))))){
    remove.ndx <- which(apply(Y.mat, 1, FUN = function(x) all(is.na(x))))
    if(nrow(rep$cov_tbl) == nrow(sample$cov_tbl)){
      rep$cov_tbl <- rep$cov_tbl[-remove.ndx,]
    } else{
      remove.df <- wide_data %>% dplyr::select(site, sample)
      keep.df <- remove.df[-remove.ndx,]

      rep$cov_tbl <- dplyr::left_join(keep.df, rep$cov_tbl, by = c('site', 'sample'))
    }
    sample$cov_tbl <- sample$cov_tbl[-remove.ndx,]
    wide_data <- wide_data[-remove.ndx,]
  }

  # build predictors
  site_mod <- as.formula(site$model)
  sample_mod <- as.formula(sample$model)
  rep_mod <- as.formula(rep$model)

  # define number of sites, samples, reps
  M <- length(unique(wide_data$site))
  J <- unname(with(wide_data, tapply(sample, site, length)))
  K <- wide_data %>%
    dplyr::select(-c(1:2)) %>%
    is.na(.) %>%
    `!` %>%
    rowSums(.)

  # factor site and sample variables
  wide_data$site <- as.factor(wide_data$site)
  wide_data$sample <- as.factor(wide_data$sample)

  # initialize z and A vectors based on data
  Y.mat <- as.matrix(wide_data[,-c(1:2)])
  y <- na.omit(c(t(Y.mat)))
  wide_data$pcr.total <- rowSums(matrix(wide_data[,3:ncol(wide_data)], nrow = nrow(wide_data)), na.rm = T)
  z.count <- unname(with(wide_data, tapply(pcr.total, site, sum)))
  z <- ifelse(z.count == 0, 0, 1)
  A <- ifelse(wide_data$pcr.total != 0, 1, 0)
  Y <- wide_data$pcr.total

  # define models
  ## site
  if(!('site' %in% names(site$cov_tbl))) stop("There should be a column named 'site' in site$cov_tbl.")
  site$cov_tbl <- dplyr::arrange(site$cov_tbl, site)
  X <- model.matrix(object = site_mod, data = site$cov_tbl)

  ## sample
  if(!('site' %in% names(sample$cov_tbl))) stop("There should be a column named 'site' in sample$cov_tbl.")
  if(!('sample' %in% names(sample$cov_tbl))) stop("There should be a column named 'sample' in sample$cov_tbl.")
  sample$cov_tbl <- dplyr::arrange(sample$cov_tbl, site, sample)
  W <- model.matrix(object = sample_mod, data = sample$cov_tbl)

  ## rep
  if(nrow(rep$cov_tbl) == sum(J)){
    if(!('site' %in% names(rep$cov_tbl))) stop("There should be a column named 'site' in rep$cov_tbl.")
    if(!('sample' %in% names(rep$cov_tbl))) stop("There should be a column named 'sample' in rep$cov_tbl.")
    rep$cov_tbl <- dplyr::arrange(rep$cov_tbl, site, sample)

    V <- model.matrix(object = rep_mod, data = rep$cov_tbl)
  } else if(nrow(rep$cov_tbl == sum(K))){
    if(!('site' %in% names(rep$cov_tbl))) stop("There should be a column named 'site' in rep$cov_tbl.")
    if(!('sample' %in% names(rep$cov_tbl))) stop("There should be a column named 'sample' in rep$cov_tbl.")
    if(!('rep' %in% names(rep$cov_tbl))) stop("There should be a column named 'rep' in rep$cov_tbl.")
    rep$cov_tbl <- dplyr::arrange(rep$cov_tbl, site, sample, rep)

    V <- model.matrix(object = rep_mod, data = rep$cov_tbl)
  } else{
    stop('The number of rows in rep$cov_tbl should be either equal to the number total number of samples or the total number of replicates.')
  }

  # setup storage for Gibbs sampler
  ## site
  ### pg
  p.site <- dim(X)[2]
  beta.mcmc <- matrix(0, num.mcmc, p.site);colnames(beta.mcmc) <- colnames(X)
  z.mcmc <- matrix(0, num.mcmc, M); z.mcmc[1,] <- z
  z.vec.mcmc <- matrix(0, num.mcmc, sum(J)); z.vec.mcmc[1,] <- rep(z, J)

  ### beta-bin
  psi.mcmc <- rep(.5, num.mcmc)

  ## sample
  ### pg
  p.samp <- dim(W)[2]
  alpha.mcmc <- matrix(0, num.mcmc, p.samp);colnames(alpha.mcmc) <- colnames(W)
  A.mcmc <- matrix(0, num.mcmc, sum(J)); A.mcmc[1,] <- A

  ### beta-bin
  if(sample$model == as.formula("~1")){
    theta.mcmc <- rep(.5, num.mcmc)
  } else if(sample$model == as.formula("~site")){
    theta.mcmc <- matrix(.5, num.mcmc, length(J))
  } else{
    theta.mcmc <- NULL
  }

  ## rep
  ### pg
  p.rep <- dim(V)[2]
  delta.mcmc <- matrix(0, num.mcmc, p.rep);colnames(delta.mcmc) <- colnames(V)

  ### beta-bin
  if(rep$model == as.formula("~1")){
    p.mcmc <- rep(.5, num.mcmc)
  } else{
    p.mcmc <- NULL
  }

  # define priors
  ## site
  mu0_site <- matrix(priors$site$mu0, p.site, 1)
  Sigma0_site <- priors$site$Sigma0* diag(p.site)
  Sigma0_site_inv <- solve(Sigma0_site)
  prior_prod_site <- Sigma0_site_inv %*% mu0_site

  ## sample
  mu0_samp <- matrix(priors$sample$mu0, p.samp, 1)
  Sigma0_samp <- priors$sample$Sigma0 * diag(p.samp)
  Sigma0_samp_inv <- solve(Sigma0_samp)
  prior_prod_samp <- Sigma0_samp_inv %*% mu0_samp

  ## rep
  mu0_rep <- matrix(priors$rep$mu0, p.rep, 1)
  Sigma0_rep <- priors$rep$Sigma0 * diag(p.rep)
  Sigma0_rep_inv <- solve(Sigma0_rep)
  prior_prod_rep <- Sigma0_rep_inv %*% mu0_rep

  ## for uniform priors
  a0 <- priors$a0
  b0 <- priors$b0

  # initialize
  ## site
  if(beta_bin & site$model == as.formula('~1')){
    psi <- rep(psi.mcmc[1], M)
  } else{
    eta <- X %*% beta.mcmc[1,]
    psi <- exp(eta) / (1 + exp(eta))
  }

  ## sample
  if(beta_bin & sample$model == as.formula("~1")){
    theta <- rep(theta.mcmc[1], sum(J))
  } else if(beta_bin & (sample$model == as.formula("~site") | sample$model == as.formula("~Site"))){
    theta <- rep(theta.mcmc[1,], J)
  } else{
    nu <- W %*% alpha.mcmc[1,]
    theta <- exp(nu) / (1 + exp(nu))
  }

  ## rep
  if(beta_bin & rep$model == as.formula("~1")){
    p <- rep(p.mcmc[1], length(y))
  } else if(dim(V)[1] == sum(J)){
    nu_rep <- V %*% delta.mcmc[1,]
    p.tmp <- exp(nu_rep) / (1 + exp(nu_rep))
    p <- rep(p.tmp, K)
  } else{
    nu_rep <- V %*% delta.mcmc[1,]
    p <- exp(nu_rep) / (1 + exp(nu_rep))
  }

  # conduct Gibbs sampler
  for(i in 2:num.mcmc){
    if(i == 2){begin <- proc.time()}

    # sample z
    A.list <- split(A, rep(1:length(J), J))
    theta.list <- split(c(theta), rep(1:length(J), J))
    z.prob <- rep(0, M)
    for(ndx in 1:M){
      if(sum(A.list[[ndx]]) > 0) {
        z.prob[ndx] <- 1
      } else{
        num <- c(psi)[ndx] * prod(1 - theta.list[[ndx]])
        denom <- num + 1 - c(psi)[ndx]
        z.prob[ndx] <- num / denom
      }
    }
    z <- rbinom(M, 1, z.prob)
    z.mcmc[i,] <- c(z)
    z.vec <- rep(z, J)
    z.vec.mcmc[i,] <- c(z.vec)

    # sample A
    y.list <- split(y, rep(1:length(K), K))
    p.list <- split(c(p), rep(1:length(K), K))
    A.prob <- rep(0, length(y.list))
    for(ndx in 1:length(y.list)){
      if(sum(y.list[[ndx]]) > 0){
        A.prob[ndx] <- 1
      } else{
        num <- z.vec[ndx] * theta[ndx] * prod(1 - p.list[[ndx]])
        denom <- num + 1 - z.vec[ndx] * theta[ndx]
        A.prob[ndx] <- num/denom
      }
    }
    A <- rbinom(sum(J), 1, A.prob)
    A.mcmc[i, ] <- A
    A.sum <- unname(sapply(split(A, rep(1:length(J), J)), sum))

    # sample beta/psi
    if(beta_bin & site$model == as.formula("~1")){
      psi.mcmc[i] <- rbeta(1, a0 + sum(z), b0 + M - sum(z))
    } else{
      # sample PG latents variables
      eta <- X %*% beta.mcmc[i-1, ]
      omega.site <- pgdraw::pgdraw(1, c(eta))

      # sample betas
      kappa.z <- z - .5
      Omega.site <- diag(omega.site, nrow = M, ncol = M)
      V.inv <- solve(t(X) %*% Omega.site %*% X + Sigma0_site_inv)
      m <- V.inv %*% (t(X) %*% kappa.z + prior_prod_site)

      beta.mcmc[i,] <- mvtnorm::rmvnorm(1, mean = m, sigma = V.inv)
    }

    # sample alpha/theta
    if(beta_bin & sample$model == as.formula("~1")){
      theta.mcmc[i] <- rbeta(1, a0 + sum(z*A.sum), b0 + sum(z*(J - A.sum)))
    } else if(beta_bin & sample$model == as.formula("~site")){
      theta.mcmc[i,] <- rbeta(M, a0 + z*A.sum, b0 + z*(J - A.sum))
    } else{
      # restrict to only where z = 1
      z.ndx <- which(z.vec == 1)
      if(dim(W)[2] == 1){
        W.red <- matrix(W[z.ndx,], ncol = 1)
        colnames(W.red) <- '(Intercept)'
      } else{
        W.red <- W[z.ndx,]
      }

      # sample PG latents variables
      nu <- W.red %*% alpha.mcmc[i-1, ]
      omega.samp <- pgdraw::pgdraw(1, c(nu))

      # sample betas
      a <- A[z.ndx]
      kappa.a <- a - .5
      Omega.samp <- diag(omega.samp)
      V.inv <- solve(t(W.red) %*% Omega.samp %*% W.red + Sigma0_samp_inv)
      m <- V.inv %*% (t(W.red) %*% kappa.a + prior_prod_samp)

      alpha.mcmc[i,] <- mvtnorm::rmvnorm(1, mean = m, sigma = V.inv)
    }

    # sample delta/p
    if(beta_bin & rep$model == as.formula("~1")){
      p.mcmc[i] <- rbeta(1, a0 + sum(A*Y), b0 + sum(A * (K - Y)))
    } else{

      # restrict to only where A = 1
      if(dim(V)[1] == sum(J)){
        A.ndx <- which(A == 1)
        if(dim(V)[2] == 1){
          V.red <- matrix(V[A.ndx,], ncol = 1)
          colnames(V.red) <- '(Intercept)'
        } else{
          V.red <- V[A.ndx,]
        }

        # sample PG latents variables
        nu_rep <- V.red %*% delta.mcmc[i-1, ]
        omega.rep <- pgdraw::pgdraw(K[A.ndx], nu_rep)

        # sample betas
        Y.vec <- Y[A.ndx]
        kappa.y <- Y.vec - K[A.ndx] / 2
        Omega.rep <- diag(omega.rep)
        V.inv <- solve(t(V.red) %*% Omega.rep %*% V.red + Sigma0_rep_inv)
        m <- V.inv %*% (t(V.red) %*% kappa.y + prior_prod_rep)

        delta.mcmc[i,] <- mvtnorm::rmvnorm(1, mean = m, sigma = V.inv)
      } else{
        A.vec <- rep(A, K)
        A.ndx <- which(A.vec == 1)
        if(dim(V)[2] == 1){
          V.red <- matrix(V[A.ndx,], ncol = 1)
          colnames(V.red) <- '(Intercept)'
        } else{
          V.red <- V[A.ndx,]
        }

        # sample PG latents variables
        nu_rep <- V.red %*% delta.mcmc[i-1, ]
        omega.rep <- pgdraw::pgdraw(1, nu_rep)

        # sample betas
        y.vec <- y[A.ndx]
        kappa.y <- y.vec - .5
        Omega.rep <- diag(omega.rep)
        V.inv <- solve(t(V.red) %*% Omega.rep %*% V.red + Sigma0_rep_inv)
        m <- V.inv %*% (t(V.red) %*% kappa.y + prior_prod_rep)

        delta.mcmc[i,] <- mvtnorm::rmvnorm(1, mean = m, sigma = V.inv)
      }

    }

    # calculate linear predictors/probs
    ## site
    if(beta_bin & site$model == as.formula('~1')){
      psi <- rep(psi.mcmc[i], M)
    } else{
      eta <- X %*% beta.mcmc[i,]
      psi <- exp(eta) / (1 + exp(eta))
    }

    ## sample
    if(beta_bin & sample$model == as.formula("~1")){
      theta <- rep(theta.mcmc[i], sum(J))
    } else if(beta_bin & (sample$model == as.formula("~site") | sample$model == as.formula("~Site"))){
      theta <- rep(theta.mcmc[i,], J)
    } else{
      nu <- W %*% alpha.mcmc[i,]
      theta <- exp(nu) / (1 + exp(nu))
    }

    ## rep
    if(beta_bin & rep$model == as.formula("~1")){
      p <- rep(p.mcmc[i], length(y))
    } else if(dim(V)[1] == sum(J)){
      nu_rep <- V %*% delta.mcmc[i,]
      p.tmp <- exp(nu_rep) / (1 + exp(nu_rep))
      p <- rep(p.tmp, K)
    } else{
      nu_rep <- V %*% delta.mcmc[i,]
      p <- exp(nu_rep) / (1 + exp(nu_rep))
    }

    if(is.null(print)){
      print <- round(.05*num.mcmc)
    }

    ## print progress
    if(progress){
      if(i %% print == 0){
        cat("\014")
        current.runtime <- unname((proc.time() - begin)[1])
        runtime <- round(current.runtime/60, 2)
        percent <- round(i/num.mcmc * 100, 2)
        message <- paste('\nIteration ', i, ' of ', num.mcmc, '; ', percent, '% done. Current runtime of ', runtime, ' minutes.\n', sep = "")
        cat(message)
        txtProgressBar(min = 2, max = num.mcmc, initial = i, style = 3)
      }
    }

  }

  empty.df <- data.frame(
    site = unique(wide_data$site) %>% rep(., J) %>% rep(., K) %>% factor(),
    sample = J %>% sapply(., FUN = function(x) 1:x, simplify = FALSE) %>% unlist() %>% rep(., K) %>% factor(),
    rep = K %>% sapply(., FUN = function(x) 1:x, simplify = FALSE) %>% unlist() %>% factor()
  )
  out <- list(beta = beta.mcmc, alpha = alpha.mcmc, delta = delta.mcmc, psi = psi.mcmc, theta = theta.mcmc, p = p.mcmc,
              model.info = list(X = X, W = W, V = V, M = M, J = J, K = K, z = z.mcmc,
                                z.vec = z.vec.mcmc, A = A.mcmc, Y = Y, y = y,
                                site_mod = site$model, samp_mod = sample$model, rep_mod = rep$model,
                                num.mcmc = num.mcmc, beta_bin = beta_bin, df = empty.df))

  # remove unused things
  if(all(out$beta == 0)) out[['beta']] <- NULL
  if(all(out$alpha == 0)) out[['alpha']] <- NULL
  if(all(out$delta == 0)) out[['delta']] <- NULL
  if(all(out$psi == .5)) out[['psi']] <- NULL
  if(all(out$theta == .5)) out[['theta']] <- NULL
  if(all(out$p == .5)) out[['p']] <- NULL

  class(out) <- c('msocc', class(out))
  return(out)
}
