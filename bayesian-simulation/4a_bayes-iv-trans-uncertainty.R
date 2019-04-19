rm(list = ls())
#source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")
# the following case only for arm size as 4
# assume the prior for each arm is the same

source("check-iv-z-arbitrary.R")
library(MCMCpack)

load(file = "general_true_ace_and_ace-bound.Rdata")

do.one.tranp.v2 = function(sample_size, prior_z = rep(1, 4)){
        if(length(prior_z) != 4){
                stop("the length of prior vector should be 4!")
        }
        arm_size = floor(sample_size / 4)
        count.vec = c(rmultinom(1, arm_size, p_observed[1:4]),
                      rmultinom(1, arm_size, p_observed[5:8]),
                      rmultinom(1, arm_size, p_observed[9:12]),
                      rmultinom(1, sample_size - 3 * arm_size, p_observed[13:16]))
        
        names(count.vec) = names(p_observed)
        
        # Specification of prior for p(x,y|z)
        # prior.z = rep(1, 4) # for each arm
        
        ## post
        post.z = rep(prior_z, 4) + count.vec
        
        num.sims <- 10000  # Number of posterior simulations for each arm
        
        theta.sims = t(replicate(num.sims, dirichlet(post.z[1:4])))
        
        colnames(theta.sims) = c("X0,Y0|Z0","X0,Y1|Z0","X1,Y0|Z0","X1,Y1|Z0")
        for(i in 1: 3){
                theta.sim.z = t(replicate(num.sims, dirichlet(post.z[(i*4+1):((i+1)*4)])))
                name = paste(c("X0,Y0", "X0,Y1", "X1,Y0", "X1,Y1"), "|Z", i, sep = "")
                colnames(theta.sim.z) = name
                theta.sims = cbind(theta.sims,theta.sim.z)
        }
        
        
        ## Remove sampled distributions that violate the IV inequalites:
        all.ivs.ok <- rep(NA,num.sims)
        
        for (i in 1:num.sims){
                all.ivs.ok[i] <- check.iv.ineq.arbz(theta.sims[i,])
        }
        
        
        posterior.theta.sims.iv <- theta.sims[all.ivs.ok,]
        
        ace.bound.mat = matrix(ncol = 2, nrow = nrow(posterior.theta.sims.iv))
        
        for(k in 1: nrow(posterior.theta.sims.iv)){
                ace.bound.mat[k,] = ace.bounds(posterior.theta.sims.iv[k,])
        }
        # save all bounds
        save(ace.bound.mat, file = paste("bounds_sample_size_", sample_size, "priorz_", prior_z[1], ".Rdata", sep = ""))
        
        N_pair = nrow(ace.bound.mat)
        
        L.vec = ace.bound.mat[,1]
        L_sort.vec = sort(L.vec)
        U.vec = ace.bound.mat[,2]
        U_sort.vec = sort(U.vec)
        
        t.vec = sort(c(L.vec, U.vec)) #increasing order
        
        start = L_sort.vec[floor(0.95 * N_pair)]
        end = U_sort.vec[ceiling(0.05 * N_pair)]
        
        start_index = which(t.vec == start)[1]
        end_index = which(t.vec == end)[1]
        
        if(start_index > end_index){
                return(c("NA", "NA"))
        }else{
                T.vec = c()
                for(i in start_index:end_index){
                        t_i = t.vec[i]
                        
                        if(sum(L.vec <= t_i) - sum(U.vec <= t_i) == ceiling(0.95 * N_pair)){
                                T.vec = c(T.vec, t_i)  
                                print(paste("i=",i))
                        }
                }
                if(length(T.vec) == 0){
                        return(c("NA", "NA")) 
                }else{
                        return(c(min(T.vec), max(T.vec)))
                }
        }
}
seed = 17
priorz = rep(1/4, 4)
set.seed(seed)

sample_size.vec = c(0, 100, 200, 400, 800, 1600)
bounds.mat = matrix(NA, ncol = 2, nrow = length(sample_size.vec))
colnames(bounds.mat) = c("LowerBounds", "UpperBounds")
rownames(bounds.mat) = sample_size.vec
for(i in 1: length(sample_size.vec)){
        size_i = sample_size.vec[i]
        print(paste("sample size = ", size_i, sep = ""))
        
        bounds.mat[i, ] = do.one.tranp.v2(size_i, prior_z = priorz)
}
save(bounds.mat, file = paste("bounds_sample_size_vary_prior_z", priorz[1], "_seed_", seed,".Rdata", sep = ""))

# thinking, as sample size increase, overlap wider, varibities for lower and upper decreases
