

bic = function(x, prel_est, type, lambda_vec = NULL, cv = 5,
               bic_const = 0.25, NumCores = NULL) {
    n = dim(x)[1]; p = dim(x)[2]
    if ( type == "graph") {
    if (is.null(lambda_vec) == T) {
      gsp = 1.5 * ( sqrt(log(p)/n) )
      lambda_vec = seq(gsp, 1, length.out = 25)
    }
  }
  snf = split_and_fit(x, prel_est, type, lambda_vec, cv)
  param_pre = snf$param_pre
  param_post = snf$param_post
  nnz = sum(param_pre!=0) + sum(param_post!=0) + 1
  bic = n * ls_eval(x, prel_est, type, param_pre,
                    param_post) + bic_const * log(n) * nnz
  # For end (no change point)
  snf_end = split_and_fit(x, t_split = 1, type, lambda_vec, cv)
  param_end = snf_end$param_pre
  nnz_end = sum(param_end != 0)
  bic_end = n * ls_eval(x, tau = 1, type,
                        param_end) + bic_const * log(n) * nnz_end
  if (bic_end < bic) { est = 1 } else { est = prel_est }
  est
}
