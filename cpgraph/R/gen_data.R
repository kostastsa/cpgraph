
gen_data = function(n, p, t0, sig1, sig2) {
  if (t0 == 1) { x = mvrnorm(n, rep(0, p), sig1) }
  if (t0 != 1) {
    n1 = round(t0 * n);
    n2 = n-n1
    u = mvrnorm(n1, rep(0, p), sig1)
    v = mvrnorm(n2, rep(0, p), sig2)
    x = rbind(u, v)
    }
  list( data = x )
  }





# example
# source("gen_paramV2.R")
#  p=250;
#  params=gen.paramV2(p=p)
# sig1=params$pre.cov; sig2=params$post.cov
#  #generate data
#  x=gen.data(n=500,p=250,t0=0.3,sig1,sig2)$data
