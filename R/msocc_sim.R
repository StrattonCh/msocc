#'@title Simulate data from a multi-scale occupancy model.
#'
#'@description This function simulates data from the multi-scale occupancy model
#'  described by \href{https://pubs.er.usgs.gov/publication/70194000}{Dorazio
#'  and Erickson (2017)}. Note that this documentation assumes there are \eqn{M}
#'  sites, \eqn{J_i} samples within each site, and \eqn{K_{ij}} replicates from
#'  each sample.
#'
#'@details This function supports both balanced an unbalanced designs as well as
#'  constant probabilities or those modeled by covariates. For unbalanced
#'  designs, \code{J} should be specified as a vector of length \code{M}
#'  representing the samples taken from each site. Similarly, \code{K} should be
#'  specified as a vector of length \code{sum(J)} representing the number of
#'  replicates taken from each sample. For balanced designs, each of \code{M},
#'  \code{J}, and \code{K} should be specified as integers. Any combination of
#'  balance or imbalance is accepted at all three levels. \cr For constant
#'  probability at any level, specify \code{psi}, \code{theta}, or \code{p}
#'  appropriately. Alternatively, the data frame, model and appropriate
#'  parameter vector must be supplied to model a probability by covariates. See
#'  the examples for more detail.
#'
#'@param M integer; number of sites to simulate
#'@param J may be either an integer or vector of length \code{M} representing
#'  the number of samples taken from each site. In the former case, \code{J} is
#'  set to \code{rep(J, M)}.
#'@param K may be either an integer or vector of length \code{sum(J)}
#'  representing the number of replicates taken from each sample. In the former
#'  case, \code{K} is set to \code{rep(K, sum(J))}.
#'@param psi may be either a single numeric value in (0,1) representing the constant
#'  probability of presence at every site or a vector of length \code{M}
#'  representing the probability of presence at each site.
#'@param theta may be either a single numeric value in (0,1) representing the
#'  constant probability of occurence in every sample or a vector of length
#'  \code{sum(J)} representing the probability of occurence in each sample.
#'@param p may be either a single numeric value in (0,1) representing the
#'  constant probability of detection in every replicate or a vector of length
#'  \code{sum(K)} representing the probability of detection in each replicate.
#'@param seed optional seed for reproducibility
#'@param site.df optional \code{data.frame} used to create the model matrix
#'  implied by \code{site.mod} if \code{psi} is a function of covariates. This
#'  \code{data.frame} must have \code{M} rows and include a column named 'site'.
#'  If specified, \code{site.mod} and \code{beta} must also be specified.
#'@param site.mod optional model statement to produce \code{psi} when it is a
#'  function of covariates. If specified, \code{site.df} and \code{beta} must
#'  also be specified.
#'@param beta optional vector of parameters use to calculate \code{psi} when it
#'  is modeled with covariates. It must have length equal to the number of
#'  columns in the model matrix implied by \code{site.mod}.
#'@param sample.df optional \code{data.frame} used to create the model matrix
#'  implied by \code{sample.mod} if \code{theta} is a function of covariates.
#'  This \code{data.frame} must have \code{sum(J)} rows and include columns
#'  named 'site' and 'sample'. If specified, \code{sample.mod} and \code{alpha}
#'  must also be specified.
#'@param sample.mod optional model statement to produce \code{theta} when it is
#'  a function of covariates. If specified, \code{sample.df} and \code{alpha}
#'  must also be specified.
#'@param alpha optional vector of parameters use to calculate \code{theta} when it
#'  is modeled with covariates. It must have length equal to the number of
#'  columns in the model matrix implied by \code{sample.mod}.
#'@param rep.df optional \code{data.frame} used to create the model matrix
#'  implied by \code{rep.mod} if \code{p} is a function of covariates. This
#'  \code{data.frame} must have either \code{sum(J)} rows (if the covariates
#'  used to model \code{p} are aggregated at the sample level) or \code{sum(K)}
#'  rows (if the covariates used to model \code{p} are at the individual
#'  replicate level). In the former case, it must include columns named 'site'
#'  and 'sample'; in the latter case, it must include columns named 'site',
#'  'sample', and 'rep'. If specified, \code{rep.mod} and \code{delta} must
#'  also be specified.
#'@param rep.mod optional model statement to produce \code{p} when it is
#'  a function of covariates. If specified, \code{rep.df} and \code{delta}
#'  must also be specified.
#'@param delta optional vector of parameters use to calculate \code{p} when it
#'  is modeled with covariates. It must have length equal to the number of
#'  columns in the model matrix implied by \code{rep.mod}.
#'
#'@example examples/msocc_sim_ex.R
#'
#'@return an object of class \code{list} containing the following elements: \cr
#' * \code{resp} a \code{data.frame} with \code{sum(J)} rows containing columns
#' 'site', 'sample', prc1, ... \cr
#' * \code{site} a \code{data.frame} with \code{M} rows containing a column
#' named 'site' and any site level covariates (if applicable) \cr
#' * \code{sample} a \code{data.frame} with \code{sum(J)} rows containing
#' columns 'site' and 'sample' and any sample level covariates (if applicable) \cr
#' * \code{rep} a \code{data.frame} with \code{sum(K)} rows containing columns
#' 'site', 'sample', and 'rep' and any rep level covariates (if applicable) \cr
#' * \code{params} a \code{list} containing the following elements: \cr
#'   * \code{psi} site level presence probabilities
#'   * \code{theta} sample level occurence probabilities
#'   * \code{p} replicate level detection probabilities
#'   * \code{z} vector of latent variables of length \code{M} denoting presence
#'   at the site
#'   * \code{a} vector of latent variables of length \code{sum(J)} denoting
#'   occurence in the sample
#'   * \code{z.vec} vector of latent variables of length \code{sum(J)} denoting
#'   site level presence stretch across samples
#' * \code{beta} matrix of regression coefficients used to generate \code{psi} if
# modeled by covariates (it is only returned in this case)
#' * \code{alpha} matrix of regression coefficients used to generate \code{theta} if
# modeled by covariates (it is only returned in this case)
#' * \code{delta} matrix of regression coefficients used to generate \code{p} if
# modeled by covariates (it is only returned in this case)
#'
#'@importFrom magrittr %>%
#'@export
#'
#'@md

msocc_sim <- function(M = 10, J = 5, K = 5, psi = .8, theta = .75, p = .9, seed = NULL,
                      site.df = NULL, site.mod = NULL, beta = NULL,
                      sample.df = NULL, sample.mod = NULL, alpha = NULL,
                      rep.df = NULL, rep.mod = NULL, delta = NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }

  if(length(J) == 1) {
    J <- rep(J, M)
  } else if(length(J) == M){
    J <- J
  } else{
    stop(paste0('J should either be vector of length 1 (denoting the number of samples at each site) or ', M))
  }

  if(length(K) == 1) {
    K <- rep(K, sum(J))
  } else if(length(K) == sum(J)){
    K <- K
  } else{
    stop(paste0('J should either be vector of length 1 (denoting the number of rep at each site and sample) or ', sum(J)))
  }

  if(!is.null(site.df)){
    if(nrow(site.df) != M) stop(paste0('site.df should have ', M, ' rows.'))
    if(!('site' %in% names(site.df))) stop("site.df should have a column named 'site'.")
    site.df <- dplyr::arrange(site.df, site) %>%
      dplyr::mutate(site = factor(site))
  }
  if(!is.null(sample.df)){
    if(nrow(sample.df) != sum(J)) stop(paste0('sample.df should have ', sum(J), ' rows.'))
    if(!('site' %in% names(sample.df))) stop("sample.df should have a column named 'site'.")
    if(!('sample' %in% names(sample.df))) stop("sample.df should have a column named 'sample'.")
    sample.df <- dplyr::arrange(sample.df, site, sample) %>%
      dplyr::mutate(site = factor(site),
                    sample = factor(sample))
  }

  if(!is.null(rep.df)){
    if(all(nrow(rep.df) != sum(J), nrow(rep.df) != sum(K)))  stop(paste0('rep.df should have ', sum(J), ' or ', sum(K), ' rows.'))
    if(nrow(rep.df) == sum(J)){
      if(!('site' %in% names(rep.df))) stop("rep.df should have a column named 'site'.")
      if(!('sample' %in% names(rep.df))) stop("rep.df should have a column named 'sample'.")
      rep.df <- dplyr::arrange(rep.df, site, sample) %>%
        dplyr::mutate(site = factor(site),
                      sample = factor(sample))
    }
    if(nrow(rep.df) == sum(K)){
      if(!('site' %in% names(rep.df))) stop("rep.df should have a column named 'site'.")
      if(!('sample' %in% names(rep.df))) stop("rep.df should have a column named 'sample'.")
      if(!('rep' %in% names(rep.df))) stop("rep.df should have a column named 'rep'.")
      rep.df <- dplyr::arrange(rep.df, site, sample, rep) %>%
        dplyr::mutate(site = factor(site),
                      sample = factor(sample),
                      rep = factor(rep))
    }
  }

  if(!is.null(site.df) & !is.null(beta) & !is.null(site.mod)){
    site.mod <- as.formula(site.mod)
    X <- model.matrix(object = site.mod, data = site.df)
    beta <- matrix(beta, ncol = 1)
    psi <- c(exp(X %*% beta) / (1 + exp(X %*% beta)))
  } else{
    psi <- psi
  }

  if(!is.null(sample.df) & !is.null(alpha) & !is.null(sample.mod)){
    sample.mod <- as.formula(sample.mod)
    W <- model.matrix(object = sample.mod, data = sample.df)
    alpha <- matrix(alpha, ncol = 1)
    theta <- c(exp(W %*% alpha) / (1 + exp(W %*% alpha)))
  } else{
    theta <- theta
  }

  if(!is.null(rep.df) & !is.null(delta) & !is.null(rep.mod)){
    if(nrow(rep.df) == sum(K)) {
      rep.mod <- as.formula(rep.mod)
      V <- model.matrix(object = rep.mod, data = rep.df)
      delta <- matrix(delta, ncol = 1)
      p <- c(exp(V %*% delta) / (1 + exp(V %*% delta)))
    } else if(nrow(rep.df) == sum(J)){
      rep.mod <- as.formula(rep.mod)
      V <- model.matrix(object = rep.mod, data = rep.df)
      delta <- matrix(delta, ncol = 1)
      p.tmp <- c(exp(V %*% delta) / (1 + exp(V %*% delta)))
      p <- rep(p.tmp, K)
    }
  } else{
    p <- p
  }

  # presence at site
  z <- rbinom(M, 1, psi)
  z.vec <- rep(z, J)

  # presence in sample
  a <- rbinom(sum(J), 1, z.vec * theta)
  a.vec <- rep(a, K)

  # presence in rep
  y <- rbinom(sum(K), 1, a.vec * p)

  # dfs
  if(is.null(site.df)) site.df <- data.frame(site = factor(1:M))
  if(is.null(sample.df)){
    sample.df <- data.frame(site = factor(rep(1:M, J)),
                            sample = J %>% sapply(., FUN = function(x) 1:x, simplify = FALSE) %>% unlist() %>% factor())
  }
  if(is.null(rep.df)){
    rep.df <- data.frame(site = rep(sample.df$site, K) %>% factor(),
                         sample = rep(sample.df$sample, K) %>% factor(),
                         rep = K %>% sapply(., FUN = function(x) 1:x, simplify = FALSE) %>% unlist() %>% factor())
  } else if(nrow(rep.df) == sum(J)){
    rep.df <- rep.df[rep(1:nrow(rep.df), K),] %>%
      dplyr::mutate(rep = K %>% sapply(., FUN = function(x) 1:x, simplify = FALSE) %>% unlist() %>% factor()) %>%
      dplyr::select(site, sample, rep, dplyr::everything())
  }

  # resp df
  y.list <- split(y, paste(rep.df$site, rep.df$sample, sep = "_") %>% factor(., levels = unique(.)))
  resp.df <- lapply(y.list, FUN = function(x) as.data.frame(t(x))) %>% dplyr::bind_rows()
  colnames(resp.df) <- paste0('pcr', 1:(ncol(resp.df)))
  resp.df <- resp.df %>%
    dplyr::mutate(site = sample.df$site,
                  sample = sample.df$sample) %>%
    dplyr::select(site, sample, dplyr::everything())

  # return
  out <- list(resp = resp.df, site = site.df, sample = sample.df, rep = rep.df,
              params = list(psi = psi, theta = theta, p = p, z = z, a = a, z.vec = z.vec))
  if(!is.null(beta)) out$params$beta <- beta
  if(!is.null(alpha)) out$params$alpha <- alpha
  if(!is.null(delta)) out$params$delta <- delta
  return(out)
}
