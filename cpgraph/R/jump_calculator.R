##jump size calculator (calculate \xi_{2,2})

jump_calc = function(sig1, sig2, NumCores = 2, par = F) {
  p = dim(sig1)[1]
  if (p != dim(sig2)[2] ) { stop("matrices incompatible") }
  # To run in parallel
  if (par == T) {
    cl <- makeCluster(NumCores);
    registerDoSNOW(cl)
    pb = txtProgressBar(max = p, style = 3)
    progress <- function(r) setTxtProgressBar(pb, r)
    opts <- list(progress = progress)
    f=foreach (j = 1:p,.options.snow = opts) %dopar% {
      mu = (solve(sig1[-j, -j])) %*% (sig1[-j, j])
      gam = (solve(sig2[-j, -j])) %*% (sig2[-j, j])
      eta = (mu - gam)
      out = list(mu = mu, gam = gam, eta = eta)
      }
    close(pb); stopCluster(cl)
    eta = matrix(0, (p-1), p)
    mu_mat = matrix(0, (p-1), p)
    gamma_mat = matrix(0, (p-1), p)
    for(j in 1:p){
      eta[ , j] = f[[j]]$eta
      mu_mat[ , j] = f[[j]]$mu
      gamma_mat[ , j] = f[[j]]$gam
      }
    eta_sq = eta ** 2
    col_sum_eta_sq = colSums(eta_sq)
    xi22 = sqrt(sum(col_sum_eta_sq))
    normalized_jump = xi22 / sqrt(p)
    }
  # To run in sequence
  if (par == F) {
    eta = matrix(0, (p-1), p);
    mu_mat = matrix(0, (p-1), p)
    gamma_mat = matrix(0, (p-1), p)
    for (j in 1:p) {
      mu_mat[ , j] = solve(sig1[-j, -j]) %*% sig1[-j, j]
      gamma_mat[ , j] = solve(sig2[-j, -j]) %*% (sig2[-j, j])
      eta[ , j] = mu_mat[ , j] - gamma_mat[ , j]
      }
    eta_sq = eta ** 2
    col_sum_eta_sq = colSums(eta_sq)
    xi22 = sqrt(sum(col_sum_eta_sq))
    normalized_jump = xi22 / sqrt(p)
  }
  list(jump_size = xi22,
       normalized_jump = normalized_jump,
       mu_mat = mu_mat,
       gamma_mat = gamma_mat )
}



####example

# source("gen_toeplitz.R")
# sig1=gen.toeplitz(rho=0.25, p=300)$cov
# sig2=gen.toeplitz(rho=0.75,p=300)$cov

# stp=proc.time()
# a=jump.calc(sig1,sig2,par = T)
# proc.time()-stp
#
# stp=proc.time()
# b=jump.calc(sig1,sig2,par = F)
# proc.time()-stp
# b$normalized.jump

