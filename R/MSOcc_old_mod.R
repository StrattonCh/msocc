#'@title Fit multi-scale occupancy model.
#'
#'@description This function fits the Bayesian multi-scale occupancy model
#'  described by Dorazio and Erickson (2017) using the polya-gamma data
#'  augmentation strategy described by Polson et al. (2012). A hierarchical
#'  implementation is also available. Note that this documentation assumes there
#'  are \code{M} sites, \code{J} samples within each site, and \code{K}
#'  replicates from each sample.
#'
#'@details This function fits variations of the multi-scale occupancy model
#'  described by Dorazio and Erickson (2017). However, this function implements
#'  a fully Bayesian sampler based on the data augmentation strategy described
#'  by Polson et al. (2012). \cr \cr It also supports hierarchical regression at
#'  the site and sample levels. Currently, only hierarchies defined by
#'  \code{site} are supported. To implement a hierarchical regression model,
#'  define the model using the typical conditional syntax. For example, to fit
#'  an intercept only regression model at the site level with random effects for
#'  each site, define \code{site$model = ~ 1 | site}. Think carefully about your
#'  prior specification, as the regression coefficients are related to
#'  probabilities through the logit link function; small changes in the
#'  magnitude of these coefficients can result in large changes in estimated
#'  probability. \cr \cr Currently, only the site and sample level occurence
#'  probabilities can be modelled by covariates. The replicate level detection
#'  probability is assumed constant.
#'
#'@param wide_data object of class \code{data.frame} containing site, sample,
#'  and PCR replicates in wide format. Column names should be \code{site},
#'  \code{sample}, \code{PCR1}, \code{PCR2}, ...
#'
#'@param site object of class \code{list} containing the following elements: \cr
#'  * \code{model} object of class formula describing the within site
#'  model. A hierarchical regression model is also available using \code{lmer}
#'  syntax; currently, only a hierarchy by site is implemented. See details. \cr
#'  * \code{cov_tbl} object of class \code{data.frame}
#'  containing site specific covariates; \code{cov_tbl} should have \code{M}
#'  rows.
#'
#'@param sample object of class \code{list} containing the following elements: \cr
#'  * \code{model} object of class formula describing the within
#'  sample model. A hierarchical regression model is also available using
#'  \code{lmer} syntax; currently, only a hierarchy by site is implemented. See details. \cr
#'  * \code{cov_tbl} object of class \code{data.frame} containing sample
#'  specific covariates; \code{cov_tbl} should have \code{sum(J)} rows.
#'
#'@param priors object of class \code{list} containing the following elements: \cr
#'  * \code{site} \verb{ } object of class \code{list} containing the
#'  following elements: \cr
#'    * \code{mu0} prior mean for site-level regression coefficients
#'  (if non-hierarchical) or site-level mu (if hierarchical) \cr
#'    * \code{Lambda0} prior covariance matrix for site-level mu \cr
#'    * \code{eta0} degrees of freedom for IW prior on Sigma \cr
#'    * \code{S0} sums of squares for IW prior on Sigma \cr
#'    * \code{Sigma0} prior covariance matrix for site-level regression coefficients (non-hierarchical)
#'  * \code{sample} \verb{ } object of class \code{list} containing the
#'  following elements: \cr
#'    * \code{mu0} prior mean for sample-level regression coefficients
#'  (if non-hierarchical) or sample-level mu (if hierarchical) \cr
#'    * \code{Lambda0} prior covariance matrix for sample-level mu \cr
#'    * \code{eta0} degrees of freedom for IW prior on Sigma \cr
#'    * \code{S0} sums of squares for IW prior on Sigma \cr
#'    * \code{Sigma0} prior covariance matrix for sample-level regression coefficients (non-hierarchical)
#'  * \code{rep} \verb{ } object of class \code{list} containing the
#'  following elements: \cr
#'    * \code{a0,b0} shape parameters for beta prior on p
#'@param num.mcmc number of MCMC samples
#'@param progress should sampling progress be printed?
#'@param print interval for printing; defaults to 5 percent of \code{num.mcmc}
#'@param seed optional seed for reproducible samples
#'@param beta_bin optional; should a beta-binomial sampler be used when possible? This option is considerably faster.
#'
#'@example examples/MSOcc_mod_ex.R
#'
#'@return object of class \code{list} containing the following elements: \cr
#'  * \code{beta} an object of class \code{matrix} of samples from the joint
#'  posterior distribution of the regression coefficients at the site level
#'  * \code{psi} an object of class \code{matrix} of samples from the joint
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
#'  * \code{p} an object of class \code{numeric} of sample from the posterior
#'  distribution of the replicate level detection probability
#'  * \code{model.info} an object of class list containing the following elements: \cr
#'    * \code{X} design matrix for site level predictors
#'    * \code{W} design matrix for sample level predictors
#'    * \code{M} number of sites
#'    * \code{J} vector of number of samples per site
#'    * \code{K} vector of number of replicates per site-sample combination, or numeric value if K is constant
#'    * \code{z} matrix of posterior samples of latent site-presence z
#'    * \code{z.vec} matrix of posterior samples of latent site-presence stretched across samples
#'    * \code{A} matrix of posterior samples of latent sample-presence A
#'    * \code{Y} observed response matrix
#'    * \code{site_mod} model statement for site-level predictors
#'    * \code{samp_mod} model statement for sample-level predictors
#'    * \code{num.mcmc} number of MCMC samples run
#'    * \code{beta_bin} was beta-binomial sampler implemented if possible?
#'
#'
#'@md

MSOcc_mod <- function(wide_data,
                      site = list(model = ~ 1, cov_tbl),
                      sample = list(model = ~ 1, cov_tbl),
                      priors = list(site = list(mu0 = 0, Lambda0 = 10, eta0 = .01, S0 = 1,
                                                Sigma0 = 25),
                                    sample = list(mu0 = 0, Lambda0 = 10, eta0 = .01, S0 = 1,
                                                  Sigma0 = 25),
                                    rep = list(a0 = 1, b0 = 1)),
                      num.mcmc = 10000, progress = T, print = NULL, seed = NULL, beta_bin = F){

  # set optional seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # fix column names
  names(wide_data)[c(1,2)] <- c('site', 'sample')

  # build predictors
  site_mod <- as.formula(site$model)
  sample_mod <- as.formula(sample$model)

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
  wide_data$pcr.total <- rowSums(wide_data[,3:ncol(wide_data)], na.rm = T)
  z.count <- unname(with(wide_data, tapply(pcr.total, site, sum)))
  z <- ifelse(z.count == 0, 0, 1)
  A <- ifelse(wide_data$pcr.total != 0, 1, 0)
  Y <- wide_data$pcr.total

  # define models
  ## site
  site_mod_char <- as.character(paste0(site$model, collapse = ""))
  site_hier <- grepl("[|]", site_mod_char)
  if(site_hier){
    ### hierarchical model
    within_model <- strsplit(site_mod_char, '[|]')[[1]][1] %>%
      stringr::str_trim(.) %>%
      as.formula(.)
    X.array <- array(0, dim = c(1, ncol(model.matrix(within_model, data = site$cov_tbl)), M))
    for(i in 1:nrow(site$cov_tbl)){
      X.array[,,i] <- model.matrix(within_model, data = site$cov_tbl[i,])
    }
    dimnames(X.array)[[3]] <- unique(wide_data$site)
  } else {
    ### non-hierarchical model
    X <- model.matrix(object = site_mod, data = site$cov_tbl)
    # rownames(X) <- paste('site', 1:M, sep = '')
    rownames(X) <- unique(wide_data$site)
  }

  ## sample
  sample_mod_char <- as.character(paste0(sample$model, collapse = ""))
  sample_hier <- grepl("[|]", sample_mod_char)
  if(sample_hier){
    samp_within_model <- strsplit(sample_mod_char, '[|]')[[1]][1] %>%
      stringr::str_trim(.) %>%
      as.formula(.)

    samp_ndx <- unique(sample$cov_tbl$site)
    W.list <- list()
    for(samp in 1:length(samp_ndx)){
      W.list[[samp]] <- model.matrix(samp_within_model, data = filter(sample$cov_tbl, site == samp_ndx[samp]))
    }
    names(W.list) <- unique(wide_data$site)
  } else {
    ### nonhierarchical model
    W <- model.matrix(object = sample_mod, data = sample$cov_tbl)
    Wnames <- paste(sample$cov_tbl$site, sample$cov_tbl$sample, sep = "_")
    rownames(W) <- Wnames
  }

  # setup storage for Gibbs sampler
  if(site_hier){
    p.site <- dim(X.array)[2]
    beta.mcmc <- array(0, dim = c(num.mcmc, p.site, M))

    z.mcmc <- matrix(0, num.mcmc, M)
    z.mcmc[1,] <- z
    z.vec.mcmc <- matrix(0, num.mcmc, sum(J))
    z.vec.mcmc[1,] <- rep(z, J)

    # beta-bin
    psi.mcmc <- rep(.5, num.mcmc)
  } else{
    # pg
    p.site <- dim(X)[2]
    beta.mcmc <- matrix(0, num.mcmc, p.site);colnames(beta.mcmc) <- colnames(X)
    z.mcmc <- matrix(0, num.mcmc, M); z.mcmc[1,] <- z
    z.vec.mcmc <- matrix(0, num.mcmc, sum(J)); z.vec.mcmc[1,] <- rep(z, J)

    # beta-bin
    psi.mcmc <- rep(.5, num.mcmc)
  }

  if(sample_hier){
    p.samp <- dim(W.list[[1]])[2]
    alpha.mcmc <- array(0, dim = c(num.mcmc, p.samp, M))

    A.mcmc <- matrix(0, num.mcmc, sum(J))
    A.mcmc[1,] <- A

    # beta-bin
    if(sample$model == as.formula("~1")){
      theta.mcmc <- rep(.5, num.mcmc)
    } else if(sample$model == as.formula("~site")){
      theta.mcmc <- matrix(.5, num.mcmc, length(J))
    } else{
      theta.mcmc <- NULL
    }
  } else {
    # pg
    p.samp <- dim(W)[2]
    alpha.mcmc <- matrix(0, num.mcmc, p.samp);colnames(alpha.mcmc) <- colnames(W)
    A.mcmc <- matrix(0, num.mcmc, sum(J)); A.mcmc[1,] <- A

    # beta-bin
    if(sample$model == as.formula("~1")){
      theta.mcmc <- rep(.5, num.mcmc)
    } else if(sample$model == as.formula("~site")){
      theta.mcmc <- matrix(.5, num.mcmc, length(J))
    } else{
      theta.mcmc <- NULL
    }
  }

  p.mcmc <- rep(.5, num.mcmc)

  # define priors
  if(site_hier){
    mu0_site <- matrix(priors$site$mu0, ncol = 1, nrow = p.site)
    Lambda0_site <- priors$site$Lambda0*diag(p.site)
    Lambda0_site_inv <- solve(Lambda0_site)
    prior_prod_site <- Lambda0_site_inv %*% mu0_site

    eta0_site <- priors$site$eta0
    S0_site <- priors$site$S0*diag(p.site)
  } else{
    mu0_site <- matrix(priors$site$mu0, p.site, 1)
    Sigma0_site <- priors$site$Sigma0* diag(p.site)
    Sigma0_site_inv <- solve(Sigma0_site)
    prior_prod_site <- Sigma0_site_inv %*% mu0_site
  }

  if(sample_hier){
    mu0_sample <- matrix(priors$sample$mu0, ncol = 1, nrow = p.samp)
    Lambda0_sample <- priors$sample$Lambda0*diag(p.samp)
    Lambda0_sample_inv <- solve(Lambda0_sample)
    prior_prod_sample <- Lambda0_sample_inv %*% mu0_sample

    eta0_sample <- priors$sample$eta0
    S0_sample <- priors$sample$S0*diag(p.samp)
  } else {
    mu0_samp <- matrix(priors$sample$mu0, p.samp, 1)
    Sigma0_samp <- priors$sample$Sigma0 * diag(p.samp)
    Sigma0_samp_inv <- solve(Sigma0_samp)
    prior_prod_samp <- Sigma0_samp_inv %*% mu0_samp
  }

  ### rep
  a0 <- priors$rep$a0
  b0 <- priors$rep$b0

  # define cutoffs to handle imbalanced designs
  samp.cutoffs <- c(0, cumsum(J))

  # initialize sampler
  theta_site <- matrix(0, ncol = 1, nrow = p.site)
  Sigma_site <- diag(p.site)
  Sigma_site_inv <- solve(Sigma_site)

  theta_sample <- matrix(0, ncol = 1, nrow = p.site)
  Sigma_sample <- diag(p.samp)
  Sigma_sample_inv <- solve(Sigma_sample)

  # conduct Gibbs sampler
  for(i in 2:num.mcmc){
    if(i == 2){begin <- proc.time()}

    #######################
    ### UPDATE LATENT Z ###
    #######################

    # linear predictor for site
    if(site_hier){
      tmp <- X.array %>%
        c(.) %>%
        matrix(., ncol = M) * beta.mcmc[i-1,,]
      eta <- colSums(tmp)
      psi <- exp(eta) / (1 + exp(eta))
    } else{
      if(beta_bin & site$model == as.formula('~1')){
        psi <- rep(psi.mcmc[i-1], M)
      } else{
        eta <- X %*% beta.mcmc[i-1,]
        psi <- exp(eta) / (1 + exp(eta))
      }
    }

    # linear predictor for sample
    if(sample_hier){
      theta <- c()
      nu <- c()
      for(samp in 1:length(samp_ndx)){
        lin_pred <- W.list[[samp]] %*% matrix(alpha.mcmc[i-1,,samp], ncol = 1)
        tmp <- exp(lin_pred) / (1 + exp(lin_pred))
        theta <- c(theta, tmp)
        nu <- c(nu, lin_pred)
      }
    } else {
      if(beta_bin & sample$model == as.formula("~1")){
        theta <- rep(theta.mcmc[i-1], sum(J))
      } else if(beta_bin & sample$model == as.formula("~site")){
        theta <- rep(theta.mcmc[i-1,], J)
      } else{
        nu <- W %*% alpha.mcmc[i-1,]
        theta <- exp(nu) / (1 + exp(nu))
      }
    }

    # determine z sampling probability
    num <- NULL;A.sum <- NULL
    for(q in 2:length(samp.cutoffs)){
      lwr <- samp.cutoffs[q-1] + 1
      upr <- samp.cutoffs[q]
      num[q-1] <- prod((1 - theta)[lwr:upr])

      A.sum[q-1] <- sum(A[lwr:upr])
    }

    z.prob.num <- psi * num ## numerator of full conditional
    z.prob.denom <- z.prob.num - psi + 1 ## denominator of full conditional
    z.prob <- z.prob.num/z.prob.denom ## calculate probability for full conditional bernoulli distr.

    z.tmp <- NULL
    for(q in 2:length(samp.cutoffs)){
      lwr <- samp.cutoffs[q-1] + 1
      upr <- samp.cutoffs[q]
      z.tmp[q-1] <- as.numeric(all(A[lwr:upr] == 0))
    }
    z.sample.prob <- ifelse(z.prob * z.tmp == 0, 1, z.prob * z.tmp) ## calculate sampling probabilities for latent z's
    # sample z
    z <- rbinom(M, size = 1, prob = z.sample.prob)
    z.mcmc[i,] <- z

    #######################
    ### UPDATE LATENT A ###
    #######################

    # create vector of presence for each sample
    z.vec <- rep(z, times = J)

    # create vector of probability of detection for each sample
    p.vec <- rep(p.mcmc[i-1], sum(J))

    A.prob.num <- z.vec * theta * (1 - p.vec)^K
    A.prob.denom <- A.prob.num + 1 - z.vec*theta
    A.prob <- A.prob.num/A.prob.denom

    y.tmp <- ifelse(Y == 0, 0, 1)
    A.prob[which(y.tmp == 1)] <- 1

    #sample A
    A <- rbinom(sum(J), size = 1, prob = A.prob)
    A.mcmc[i,] <- A

    #############################################
    ### UPDATE REGRESSION COEFFICIENTS - BETA ###
    #############################################

    if(site_hier){
      # sample hyper-parameters
      ## theta
      ### theta varcov
      Lambda_site <- Lambda0_site_inv + M*Sigma_site_inv
      Lambda_site_inv <- solve(Lambda_site)

      ### theta mean
      beta_bar <- beta.mcmc[i-1,,] %>%
        matrix(., ncol = M) %>%
        rowMeans(.) %>%
        matrix(., ncol = 1)
      mu_site <- Lambda_site_inv %*% (prior_prod_site + M*Sigma_site_inv %*% beta_bar)

      ### sample theta
      theta_site <- mvtnorm::rmvnorm(1, mean = mu_site, sigma = Lambda_site_inv)

      ## sigma
      eta_site <- eta0_site + M
      tmp <- beta.mcmc[i-1,,] %>%
        matrix(., ncol = M) %>%
        apply(., 2, FUN = function(x) x - theta_site) %>%
        matrix(., ncol = M) %>%
        rowSums(.) %>%
        matrix(., ncol = 1)
      S_theta_site <- tmp %*% t(tmp)

      Sigma_site <- MCMCpack::riwish(v = eta_site, S = solve(S0_site + S_theta_site))
      Sigma_site_inv <- solve(Sigma_site)

      # sample betas
      ## sample PG latents variables
      omega.site <- BayesLogit::rpg(M, rep(1, M), eta)

      ## define mean and variance
      z.star <- (1/omega.site)*(z - .5)
      z.star.array <- array(z.star, dim = c(1,1,M))

      for(loc in 1:M){
        X.site <- matrix(X.array[,,loc], 1, p.site)
        Omega.site <- matrix(omega.site[loc], 1, 1)
        V.inv <- solve(t(X.site) %*% Omega.site %*% X.site + Sigma_site_inv)
        m <- V.inv %*% (t(X.site) %*% Omega.site %*% z.star.array[,,loc] + Sigma_site_inv %*% mu_site)
        beta.mcmc[i,,loc] <- mvtnorm::rmvnorm(1, mean = m, sigma = V.inv)
      }
    } else{
      if(beta_bin & site$model == as.formula("~1")){
        psi.mcmc[i] <- rbeta(1, a0 + sum(z), b0 + M - sum(z))
      } else{
        # sample PG latents variables
        eta <- X %*% beta.mcmc[i-1, ]
        omega.site <- BayesLogit::rpg(M, rep(1, M), eta)

        # sample betas
        z.star <- (1/omega.site)*(z - .5)
        Omega.site <- diag(omega.site, nrow = M, ncol = M)
        V.inv <- solve(t(X) %*% Omega.site %*% X + Sigma0_site_inv)
        m <- V.inv %*% (t(X) %*% Omega.site %*% z.star + prior_prod_site)

        beta.mcmc[i,] <- mvtnorm::rmvnorm(1, mean = m, sigma = V.inv)
      }
    }

    ##############################################
    ### UPDATE REGRESSION COEFFICIENTS - ALPHA ###
    ##############################################

    if(sample_hier){
      # sample hyper-parameters
      ## theta
      ### theta varcov
      Lambda_sample <- Lambda0_sample_inv + M*Sigma_sample_inv
      Lambda_sample_inv <- solve(Lambda_sample)

      ### theta mean
      alpha_bar <- alpha.mcmc[i-1,,] %>%
        matrix(., ncol = M) %>%
        rowMeans(.) %>%
        matrix(., ncol = 1)
      mu_sample <- Lambda_sample_inv %*% (prior_prod_sample + M*Sigma_sample_inv %*% alpha_bar)

      ### sample theta
      theta_sample <- mvtnorm::rmvnorm(1, mean = mu_sample, sigma = Lambda_sample_inv)

      ## sigma
      eta_sample <- eta0_sample + M
      tmp <- alpha.mcmc[i-1,,] %>%
        matrix(., ncol = M) %>%
        apply(., 2, FUN = function(x) x - theta_sample) %>%
        matrix(., ncol = M) %>%
        rowSums(.) %>%
        matrix(., ncol = 1)
      S_theta_sample <- tmp %*% t(tmp)

      Sigma_sample <- MCMCpack::riwish(v = eta_sample, S = solve(S0_sample + S_theta_sample))
      Sigma_sample_inv <- solve(Sigma_sample)

      # sample alphas
      ## sample PG latents variables
      omega.samp <- BayesLogit::rpg(sum(J), rep(1, sum(J)), nu)
      a.star <- (1/omega.samp)*(A - .5)

      ## sample alphas
      a.star.tbl <- tibble(
        site = sample$cov_tbl$site,
        a.star = a.star,
        omega.samp = omega.samp
      )

      tmp_prod <- Sigma_sample_inv %*% mu_sample
      for(loc in 1:length(unique(a.star.tbl$site))){
        W.samp <- matrix(W.list[[loc]], ncol = p.samp)
        a.star.site <- filter(a.star.tbl, site == samp_ndx[loc]) %>%
          dplyr::select(a.star) %>%
          unlist(.) %>%
          unname(.) %>%
          matrix(., ncol = 1)
        Omega.samp <- filter(a.star.tbl, site == samp_ndx[loc]) %>%
          dplyr::select(omega.samp) %>%
          unlist(.) %>%
          unname(.) %>%
          diag(.)

        V.inv <- solve(t(W.samp) %*% Omega.samp %*% W.samp + Sigma_sample_inv)
        m <- V.inv %*% (t(W.samp) %*% Omega.samp %*% a.star.site + tmp_prod)

        alpha.mcmc[i,,loc] <- mvtnorm::rmvnorm(1, mean = m, sigma = V.inv)
      }

    } else {
      if(beta_bin & sample$model == as.formula("~1")){
        theta.mcmc[i] <- rbeta(1, a0 + sum(z.vec*A), b0 + sum(z*(J - A.sum)))
      } else if(beta_bin & sample$model == as.formula("~site") | sample$model == as.formula("~Site")){
        # z*A.sum problematic if z = 0???
        theta.mcmc[i,] <- rbeta(M, a0 + z*A.sum, b0 + z*(J - A.sum))
      } else{
        # restrict to only where z = 1
        z.ndx <- which(z.vec == 1)
        if(dim(W)[2] == 1){
          W.red <- matrix(W[z.ndx,], ncol = 1)
        } else{
          W.red <- W[z.ndx,]
        }

        # sample PG latents variables
        nu <- W.red %*% alpha.mcmc[i-1, ]
        omega.samp <- BayesLogit::rpg(dim(W.red)[1], rep(1, dim(W.red)[1]), nu)

        # sample betas
        a <- A[z.ndx]
        a.star <- (1/omega.samp)*(a - .5)
        Omega.samp <- diag(omega.samp)
        V.inv <- solve(t(W.red) %*% Omega.samp %*% W.red + Sigma0_samp_inv)
        m <- V.inv %*% (t(W.red) %*% Omega.samp %*% a.star + prior_prod_samp)

        alpha.mcmc[i,] <- mvtnorm::rmvnorm(1, mean = m, sigma = V.inv)
      }
    }

    ######################################
    ### UPDATE REGRESSION COEFFICIENTS ### for rep to rep covariates - assumed to be constant
    ######################################
    p.mcmc[i] <- rbeta(1, a0 + sum(A*Y), b0 + sum(A * (K - Y)))

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

  if(site_hier & sample_hier){
    out <- list(beta = beta.mcmc, alpha = alpha.mcmc, p = p.mcmc, psi = psi.mcmc, theta = theta.mcmc,
                model.info = list(X = X.array, W = W.list, M = M, J = J, K = K, z = z.mcmc,
                                  z.vec = z.vec.mcmc, A = A.mcmc, Y = Y,
                                  site_mod = site_mod_char, samp_mod = sample_mod_char,
                                  num.mcmc = num.mcmc, beta_bin = beta_bin))
  } else if(site_hier & !sample_hier){
    out <- list(beta = beta.mcmc, alpha = alpha.mcmc, p = p.mcmc, psi = psi.mcmc, theta = theta.mcmc,
                model.info = list(X = X.array, W = W, M = M, J = J, K = K, z = z.mcmc,
                                  z.vec = z.vec.mcmc, A = A.mcmc, Y = Y,
                                  site_mod = site_mod_char, samp_mod = sample$model,
                                  num.mcmc = num.mcmc, beta_bin = beta_bin))

  } else if(!site_hier & sample_hier){
    out <- list(beta = beta.mcmc, alpha = alpha.mcmc, p = p.mcmc, psi = psi.mcmc, theta = theta.mcmc,
                model.info = list(X = X, W = W.list, M = M, J = J, K = K, z = z.mcmc,
                                  z.vec = z.vec.mcmc, A = A.mcmc, Y = Y,
                                  site_mod = site$model, samp_mod = sample_mod_char,
                                  num.mcmc = num.mcmc, beta_bin = beta_bin))

  } else if(!site_hier & !sample_hier){
    out <- list(beta = beta.mcmc, alpha = alpha.mcmc, p = p.mcmc, psi = psi.mcmc, theta = theta.mcmc,
                model.info = list(X = X, W = W, M = M, J = J, K = K, z = z.mcmc,
                                  z.vec = z.vec.mcmc, A = A.mcmc, Y = Y,
                                  site_mod = site$model, samp_mod = sample$model,
                                  num.mcmc = num.mcmc, beta_bin = beta_bin))
  }

  # remove unused things
  if(all(out$beta == 0)) out[['beta']] <- NULL
  if(all(out$alpha == 0)) out[['alpha']] <- NULL
  if(all(out$psi == .5)) out[['psi']] <- NULL
  if(all(out$theta == .5)) out[['theta']] <- NULL

  return(out)
}
