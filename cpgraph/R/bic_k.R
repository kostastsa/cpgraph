

bic_k = function(x, t_split, k_seq = seq(0, 0.5, length.out = 50)) {
  n = dim(x)[1]; p = dim(x)[2]
  bic_list = c()
  cut = floor(n * t_split)
  x_pre = x[(1:cut), ]
  if ( t_split != 1 ){
    x_post = x[((cut+1):n), ]
    for (k in k_seq) {
      param_pre = soft_thresh(x_pre, k)
      param_post = soft_thresh(x_post, k)
      nnz = sum(param_pre!=0) + sum(param_post!=0)
      bic = n * ls_eval(x, t_split, type = "mean",
                        param_pre, param_post) + log(n) * nnz
      bic_list = c(bic_list, bic)
    }
    ind = which.min(bic_list)
  } else {
      for (k in k_seq) {
        param = soft_thresh(x_pre, k)
        nnz = sum(param!=0)
        bic = n * ls_eval(x, t_split, type = "mean", param) + log(n) * nnz
        bic_list = c(bic_list, bic)
      }
    ind = which.min(bic_list)
    }
  k_seq[ind]
  }
