#'@title Compute waic
#'
#'@description This function computes the waic described in Gelman et al. (2013)
#'  for multi-scale occupancy models.
#'
#'@param msocc_mod output from \code{\link{msocc_mod}}
#'@param type one of \code{c(1, 2)} denoting the type of penalty to use when
#'  calculating the waic
#'
#'@return numeric value that is the waic
#'
#'@details The authors of Gelman et al. (2013) note that the type 2 penalty is a
#'  better representation of leave one out cross-validation, and therefore
#'  recommend its use. \cr In the case of hierarchical models, they also note
#'  that there are two ways to think the likelihood; one that incorporates the
#'  hyper-parameters and one that does not. Both are arguably justifiable
#'  depending upon the situation. We do not incorporate the hyper-parameters in
#'  our calculations here.
#'
#'@example examples/waic_ex.R
#'@export

waic <- function(msocc_mod, type = 2){
  if(type == 1){
    waic <- -2 * (compute_lppd(msocc_mod) - compute_pwaic1(msocc_mod))
    return(waic)
  } else{
    if(type == 2){
      waic <- -2 * (compute_lppd(msocc_mod) - compute_pwaic2(msocc_mod))
      return(waic)
    } else{
      stop('WAIC type should either be 1 or 2. See help page for more information.')
    }
  }
}
