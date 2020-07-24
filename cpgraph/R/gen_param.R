###create two distinct cluster type toeplitz covariances

#s: Number of non-zero components in each column of each covariance matrix

gen_param = function(p, rho1 = 0.5, rho2 = 0.5, root1 = 4, root2 = 4, s1 = 5,
                     s2 = 5, par = F, reg_param = F, NumCores =2 , cvar = 2) {
        mod1 = cluster_toeplitz(p = p, rho = rho1, root = root1, s = s1)
        cov1 = cvar * mod1$cov
        s3 = 2 * s2
        mod2 = cluster_toeplitz(p = p, rho = rho2, root = root2, s = s3)
        cov2 = mod2$cov
        # Create structural changes withing covariances by alternating signs
        A1 = (-1) ** toeplitz(0:(p-1))
        sig1 = A1 * cov1
        A2 = ifelse((A1 == -1), 0, 1)
        sig2 = A2 * cov2
        # To check if both matrices are positive definite
        egv = eigen(sig1)$values
        if (min(egv)<=0) { print("sig1 is not positive definite") }
        egv = eigen(sig2)$values
        if (min(egv)<=0) { print("sig2 is not positive definite") }
        jump_size = NULL
        normalized_jump = NULL
        mu_mat = NULL
        gamma_mat = NULL
        if (reg_param == T) {
                jump = jump_calc(as.matrix(sig1), as.matrix(sig2), par = par,
                                 NumCores = NumCores)
                jump_size = jump$jump_size
                normalized_jump = jump$normalized_jump
                mu_mat = jump$mu_mat
                gamma_mat = jump$gamma_mat
                }
        list( pre_cov = sig1,
              post_cov = sig2,
              jump_size = jump_size,
              normalized_jump = normalized_jump,
              mu_mat = mu_mat,
              gamma_mat = gamma_mat
              )
}


###example

 # mod=gen.paramV2(p=100)
 # mod$normalized.jump
 # image(mod$pre.cov)
 # image(mod$post.cov)
 #
 #mod=gen.paramV2(200,reg.param = T)
 #mod$normalized.jump
 # mu.mat=mod$mu.mat
 # mu.mat[,1]
 #
 # mod=gen.paramV2(50,rho1=0.75,rho2=0.35,root1=4,root2=2,reg.param = T,par=T)
 # mod$normalized.jump
 # mu.mat=mod$mu.mat
 # mu.mat[,1]
 # mod$gamma.mat[,1]









