####generate a toeplitz type single covariance structure


gen_toeplitz = function(p, rho = 0.75, root = 2) {
  toep = ( toeplitz(0:(p-1)) ) ** ( 1 / root )
  sig = rho ** toep
  # To check if positive definite
  egv = eigen(sig)$values
  if (min(egv) <= 0) { print("Toep matrix is not positive definite") }
  list(cov = sig)
}



###example

# mod=gen.toeplitz(p=10)
# mod=gen.toeplitz(p=100,root=3)
# image(mod$cov)
