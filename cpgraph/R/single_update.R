#'@title Single iteration of the split, fit and minimize procedure.
#'@description Finds the minimizer of the least squares (LS) cost given a data
#'matrix and initial split location or grid thereof.
#'It applies an information criterion with user-chosen penalty.
#'It returns the estimated value with the corresponding cost as well as the
#'initial split location.
#'@details This function is a single application of the split, fit and minimize
#'procedure. The procedure initially splits the data at an initial estimate
#'\code{start_grid}. Then in each segment it fits a sequence of linear models
#'to estimate precision matrices for each segment. Finally it minimizes
#'the least squares cost, that is defined throught these estimates, to obtain
#'an updated estimate of the change point. If the user inputs \code{start_grid}
#'as a vector of values then a split value is chosen that has minimal LS cost
#'with parameters fitted at this split value. The parameters are obtained using
#' the \code{glmnet} regularization path with parameters specified in
#'\code{lambda_vec}. If \code{lambda_vec} is omitted the path is set to
#'\deqn{seq(lam_sep * ( sqrt( log(max(p, n) ) / n ), l_end, length.out = length_lam)}
#'If \code{par=TRUE} the coefficient estimation is done in parallel using the
#'library \code{parallel} with the number of CPU cores set to \code{NumCores}.
#'@author Konstantinos Tsampourakis <kostas.tsampourakis@@gmail.com>
#'@param x an n by p matrix of observations.
#'@param type a string of characters denoting the detection type : "mean" or
#' "graph"
#'@param edge percentage of throaway data from the edges.
#'@param lambda_vec a grid of values of the regularization parameter to be used
#' for cross-validation.
#'@param length_lam length of grid.
#'@param lam_sep numeric value. Used in the regularizer grid if not provided by
#' the user.
#'@param par If TRUE parallelization of the lasso estimation is performed.
#'@param cv number of cross-validation folds.
#'@param bic_const penalty used in the BIC evaluation.
#'@param NumCores number of processor cores to be used
#'@param start_grid grid of values for the initial split location
#'@param l_end end value of the lambda sequence
#'@return
#'\item{est}{numerical value in (0,1]. If est = 1 then no
#' change-point was detected.
#' If est<1 the value corresponds to the LS minimizer.}
#'\item{tau_check}{the value used for the initial data split.}
#'\item{loss.grid}{LS cost values along the time grid.}
#'@export
#'@examples
#' p = 50; n = 300; t0 = 0.2
#' params = gen_param(p = p, reg_param = TRUE)
#' sig1 = params$pre_cov; sig2 = params$post_cov
#' image(sig1); image(sig2)
#' # Generate data
#' x = gen_data(n = n, p = p, t0 = t0, sig1, sig2)$data
#' params$jump_size
#' params$normalized_jump
#'
#' # Sequential
#' mod = single_update(x)
#' plot(mod$loss_grid,type="l")
#' mod$est
#'
#' # Parallel
#' mod = single_update(x, par = T)
#' plot(mod$loss_grid, type="l")
#' mod$est
#'
#' # Parallel and search grid
#' mod = single_update(x, par = T, start_grid = c(0.3,0.6))
#' plot(mod$loss_grid, type="l")
#' mod$est

single_update = function(x, type = "mean", edge = 0.05, lambda_vec = NULL,
                         length_lam = 50, lam_sep = 1.75, par = F, cv = 5,
                         bic_const = 0.25, NumCores = 4, start_grid = 0.5,
                         l_end = 2) {
  n = dim(x)[1]; p = dim(x)[2]; pvn = max(p,n)
  # To set regularizer grid if it is not defined in function
  if ( type == "graph" ) {
    if (is.null(lambda_vec) == T) {
      gsp = lam_sep * ( sqrt(log(pvn) / n) )
      lambda_vec = seq(gsp, l_end, length.out = length_lam)
    }
  }

  # To search in start_grid for best starting value
  ls_val = c(); param_list_pre = list(); param_list_post = list()
  lsg = length(start_grid)
  for(j in 1:lsg) {
    t_split = start_grid[j]
    if (t_split >= 1 || t_split <= 0) {
      warning("split value out of bounds")
    }
    snf = split_and_fit(x, t_split, type, lambda_vec, cv, NumCores,
                              par, length_lam, lam_sep, l_end)
    param_pre = snf$param_pre
    param_post = snf$param_post
    param_list_pre[[j]] = param_pre
    param_list_post[[j]] = param_post
    ls_val[j] = ls_eval(x, t_split, type, param_pre, param_post)
  }
  min_val = which.min(ls_val)
  sp = start_grid[min_val]
  param_pre = param_list_pre[[min_val]]
  param_post = param_list_post[[min_val]]
  loss_grid = ls_eval(x, tau="all", type, param_pre, param_post)
  prel_est = ( which.min(loss_grid) ) / n
  # To check if estimate is too close to edge and sdeclare no change point
  if ((prel_est < edge) | (prel_est > (1 - edge))) { ret_est = 1 }
  # To check if not and use BIC
  if ( (prel_est >= edge) & (prel_est <= 1 - edge) ) {
    ret_est = bic(x, prel_est, type, lambda_vec, cv, bic_const, NumCores)
  }
  list( est = ret_est,
        tau_check = sp,
        loss_grid = loss_grid )
}

