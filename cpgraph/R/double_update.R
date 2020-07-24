#'@title Double iteration of the split, fit and minimize procedure.
#'@description Finds the optimal estimator of the change point by applying
#'twice the split, fit and minimize procedure on a sequence of observations.
#'@details This function is a double application of the split, fit and
#'minimize procedure.
#'If the user inputs \code{start_grid} as a vector of values then a split value
#'is chosen that has minimal LS cost with parameters
#'fitted at this split value. The parameters are obtained using the
#'\code{glmnet} regularization path with parameters specified in
#'\code{lambda_vec}.
#'If \code{lambda_vec} is omitted the path is set to
#'\deqn{seq(lam_sep * ( sqrt( log(max(p,n) ) / n ), l_end, length.out = length_lam)}
#'If \code{par=TRUE} the coefficient estimation is done in parallel using
#'the library \code{parallel} with the number of CPU cores
#'set to \code{NumCores}.
#'@author Konstantinos Tsampourakis <kostas.tsampourakis@@gmail.com>
#'@inheritParams single_update
#'@return
#'\item{upd1}{First estimate (sub-optimal)}
#'\item{upd2}{Second estimate (optimal)}
#'@export
#'@examples
#' p = 50; n = 225; t0 = 0.2
#' params = gen_param(p = p, reg_param = TRUE)
#' sig1 = params$pre_cov; sig2 = params$post_cov
#' image(as.matrix(sig1));image(as.matrix(sig2))
#' #generate data
#' x=gen.data(n = n, p = p, t0 = t0, sig1, sig2)$data
#' params$jump_size
#' params$normalized_jump
#'
#' #parallel
#' stp = proc.time()
#' mod_full = double_update(x, par = TRUE)
#' proc.time() - stp
#' mod_full
#'
#' #sequential
#' stp = proc.time()
#' mod_full = double_update(x)
#' proc.time() - stp
#' mod_full


double_update = function(x, type = "mean", edge = 0.05, lambda_vec = NULL,
                         length_lam = 50, lam_sep = 1.75, par = F, cv = 5,
                         bic_const = 0.25, NumCores = 4, start_grid = 0.5,
                         l_end = 2) {

  tau_hat = single_update(x, type, edge, lambda_vec, length_lam,
                          lam_sep, par, cv, bic_const, NumCores, start_grid,
                          l_end)$est
  # To check if step 1 estimate at the boundary
  if (tau_hat <= edge | tau_hat >= (1-edge)) {
    tau_tilde = tau_hat
  } else {
    n = dim(x)[1] ; p = dim(x)[2]; pvn = max(p, n)
    if ( type == "graph") {
      if (is.null(lambda_vec) == T) {
        gsp = lam_sep * ( sqrt(log(pvn) / n) )
        lambda_vec = seq(gsp, l_end, length.out = length_lam)
      }
    }
    snf = split_and_fit(x, tau_hat, type, lambda_vec, cv, NumCores, par)
    param_pre = snf$param_pre;
    param_post = snf$param_post
    loss_grid = ls_eval(x, tau = "all", type, param_pre, param_post)
    prel_est = ( which.min(loss_grid) ) / n
    tau_tilde = prel_est
    # To check if estimate is beyond edge
    if (prel_est < edge) { tau_tilde = edge }
    if (prel_est > (1-edge)) { tau_tilde = 1 - edge }
    }
  list( upd1 = tau_hat,
        upd2 = tau_tilde)
}







