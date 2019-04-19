source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")
source("check-iv-z-arbitrary.R") 
library(MCMCpack)

load(file = "bayesian-simulation/true_ace_and_ace-bound.Rdata")

do.one = function(){
        count.vec = c(rmultinom(1, 25, p_observed[1:4]), 
                      rmultinom(1, 25, p_observed[5:8]), 
                      rmultinom(1, 25, p_observed[9:12]),
                      rmultinom(1, 25, p_observed[13:16]))
        
        names(count.vec) = names(p_observed)
        
        # Specification of prior for p(x,y|z)
        prior.z = rep(1, 4) # for each arm
        
        ## post
        post.z = rep(prior.z, 4) + count.vec
        
        num.sims <- 10000  # Number of posterior simulations for each arm
        
        theta.sims = t(replicate(num.sims, dirichlet(post.z[1:4])))
        
        colnames(theta.sims) = c("X0,Y0|Z0","X0,Y1|Z0","X1,Y0|Z0","X1,Y1|Z0")
        for(i in 1: 3){
                theta.sim.z = t(replicate(num.sims, dirichlet(post.z[(i*4+1):((i+1)*4)])))
                name = paste(c("X0,Y0", "X0,Y1", "X1,Y0", "X1,Y1"), "|Z", i, sep = "")
                colnames(theta.sim.z) = name
                theta.sims = cbind(theta.sims,theta.sim.z)
        }
        
        summary(theta.sims)
        
        ## Remove sampled distributions that violate the IV inequalites:
        all.ivs.ok <- rep(NA,num.sims)
        
        for (i in 1:num.sims){
                all.ivs.ok[i] <- check.iv.ineq.arbz(theta.sims[i,])
        }
        
        summary(all.ivs.ok)
        
        mean(all.ivs.ok) 
        
        posterior.theta.sims.iv <- theta.sims[all.ivs.ok,]
        
        # likelihood for the posterier
        log_likelihood = apply(posterior.theta.sims.iv, 1, 
                               function(x) log(ddirichlet(x[1:4], post.z[1:4]) * 
                                                       ddirichlet(x[5:8], post.z[5:8]) * 
                                                       ddirichlet(x[9:12], post.z[9:12]) *
                                                       ddirichlet(x[13:16], post.z[13:16])))
        # remove bottom 5%
        posterior.theta.sims.iv.top95 = posterior.theta.sims.iv[which(rank(log_likelihood) >= nrow(posterior.theta.sims.iv) * 0.05),]
        
        ace.bound.mat = matrix(ncol = 2, nrow = nrow(posterior.theta.sims.iv.top95))
        for(k in 1: nrow(posterior.theta.sims.iv.top95)){
                ace.bound.mat[k,] = ace.bounds(posterior.theta.sims.iv.top95[k,])
        }
        
        # coverage of true
        n_cover = length(which(ace_true < ace.bound.mat[,2] & ace_true > ace.bound.mat[,1]))
        n_hpd95 = nrow(ace.bound.mat)
        # union of the interval
        return(c(min(ace.bound.mat[,1]), max(ace.bound.mat[,2]), n_cover, n_hpd95, n_cover/n_hpd95))
}
set.seed(17)
mat = replicate(1000,do.one())
rownames(mat) = c("min_LB", "max_UB", "NO.ace_covered", "NO.top95", "Prop")

save(mat, file = "sim1000.R")
