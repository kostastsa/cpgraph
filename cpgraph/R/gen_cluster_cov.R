#' #create single cluster type covariance matrix with fixed sparsity.
#' #s is number of non zero components in each column of covariance
#' @importFrom Matrix bdiag


cluster_toeplitz = function(p, rho = 0.75, root = 2, s = 6){
  ##length of each cluster preserving sparsity (approximately)
  nclus = ceiling(p / (s + 1)); cut_clus = c(0, (s+1) * (1:(nclus - 1)), p)
  #create cluster structure
  diff = cut_clus[2] - cut_clus[1]
  A = matrix(1, diff, diff)
  for(j in 2:nclus) {
    diff = cut_clus[j + 1] - cut_clus[j]
    Asub = matrix(1, diff, diff)
    A = bdiag(A, Asub)
  }
  # To obtain toeplitz
  mod = gen_toeplitz(p = p, rho = rho, root = root)
  # To convert to cluster structure
  cov = mod$cov;  cluster_cov = A * cov
  # To check for positive definite
  egv = eigen(cluster_cov)$values;
  if (min(egv) <= 0) { print("Toep matrix is not positive definite") }
  list( cov = cluster_cov)
  }


##example
# mod=cluster.toeplitz(p=200,s=10)
# cov=mod$cov
# image(cov)

