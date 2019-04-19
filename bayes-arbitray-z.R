rm(list = ls())
library(dplyr)
library(ggplot2)
library(xtable)
library(hexbin)
source("check-iv-z-arbitrary.R")
#source("http://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew.R")
#---------------------------------------------------#
# Bayesion Method of IV model, arbitray level of z
#
# Notation:
# qij = P(Y=i, D=j|Z=z), p_ij = P(Y0=i, Y1=j), i,j = 0,1
# IV Inequality:
#     q11.bar + q01.bar <= 1;
#     q10.bar + q00.bar <= 1;
#     q10.bar <= q1.lowbar + q10-01.lowbar;
#     q00.bar <= q0.lowbar + q00-11.lowbar;
#     q11.bar <= q1.lowbar + q00-11.lowbar;
#     q01.bar <= q0.lowbar + q10-01.lowbar;
#     1 <= q1.lowbar + q10-01.lowbar + q0.lowbar + q00-11.lowbar
#---------------------------------------------------#

#---------------------------------------------------#
# Simulation, RCT
# X has no effect on Y
#---------------------------------------------------#
n_sample = 5000
z_level = 3

set.seed(1)
u = rnorm(n = n_sample, 0, 1)
z_true = sample(x = 0:(z_level-1), replace = T, size = n_sample)

x_true = rbinom(n = n_sample,
                prob = exp(1 * u + 0.5 * z_true)/(1 + exp(1 * u + 0.5 * z_true)),
                size = 1)
y_true = rbinom(n = n_sample,
                prob = exp(3 * u)/(1 + exp(3 * u)),
                size = 1)

cor(z_true, y_true)
cor(z_true, x_true)

data_true = data.frame(z = z_true, x = x_true, y = y_true)

count_true = data_true %>% group_by(z, x, y) %>%
      summarise(n = n()) %>%
      as.data.frame()

print(xtable(count_true, caption = "Generation of Observed Data"), include.rownames=FALSE)

arm_size = data_true %>% group_by(z) %>%
        summarise(n = n()) %>%
        as.data.frame()

data_true = merge(count_true, arm_size, by = "z")

pyx.z_true = data_true[,4] / data_true[,5]

ace.bounds(pyx.z_true)
# note that z=3, x=0, y = 1, the count is zero
#count_true = rbind(count_true[1:13,], c(3,0,1,0), count_true[14:15, ])

#ace_true = (sum(count_true[(1:zlevel) * 4,4])-sum(count_true[(0:(zlevel-1)) * 4 + 2, 4]))/n_sample
#---------------------------------------------------#
# Bayesian Inference
#---------------------------------------------------#

# Specification of prior for p(x,y|z)
prior.z = rep(1, 4) # for each arm

## post
post.z = rep(prior.z, z_level) + count_true[, 4]

num.sims <- 10000  # Number of posterior simulations for each arm

theta.sims = t(replicate(num.sims, dirichlet(post.z[1:4])))
colnames(theta.sims) = c("X0,Y0|Z0","X0,Y1|Z0","X1,Y0|Z0","X1,Y1|Z0")
for(i in 1: (z_level-1)){
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

### posterior.theta.sims.iv contains simulations from the posterior
### distribution of p(x,y|z) under our model

n <- sum(all.ivs.ok)
n # number of sims remaining after removing those
# violating IV inequalities


#ace = apply(posterior.theta.sims.iv, 1, function(x) sum(x[(1:z_level) * 4])-sum(x[(0:(z_level-1)) * 4 + 2]))/4

#plot(density(ace))

ace.bound.mat = matrix(ncol = 2, nrow = num.sims)
for(k in 1: num.sims){
        ace.bound.mat[k,] = ace.bounds(posterior.theta.sims.iv[k,])
}
png("figure/no-effec-tri-plot.png")
do.tri.plot(ace.bound.mat, title.txt = "ACE Bounds")
dev.off()

h <- hexbin(ace.bound.mat)

png("figure/no-effec-hexbin-plot.png")
plot(h, xlab = "Lower Bound on ACE", ylab = "Upper Bound on ACE")
dev.off()

# ace.bound.df = as.data.frame(ace.bound.mat)
# colnames(ace.bound.df) = c("lower", "upper")
# p = ggplot(ace.bound.df, aes(lower, upper))
# p + stat_bin2d(bins = 40) +
#       scale_fill_gradient(low = "grey", high = "black")
#


# use Thomas' codes
# change the matrix a little bit, swith position. make pyx.z to pxy.z
T_posterior.theta.sims.iv = posterior.theta.sims.iv
tmp = posterior.theta.sims.iv[,(c(1: z_level) - 1) * 4 + 3]
T_posterior.theta.sims.iv[,(c(1: z_level) - 1) * 4 + 3] = posterior.theta.sims.iv[,(c(1: z_level) - 1) * 4 + 2]
T_posterior.theta.sims.iv[,(c(1: z_level) - 1) * 4 + 2] = tmp

T_ace.bound.mat = matrix(ncol = 2, nrow = num.sims)
for(k in 1: num.sims){

        T_ace.bound.mat[k,] = ace.bounds.tr(T_posterior.theta.sims.iv[k,])
}

do.tri.plot(T_ace.bound.mat, title.txt = "ACE Bounds")
T_h <- hexbin(T_ace.bound.mat)
plot(T_h)

#--------------------------------------------------------------------------#
# X has effect on Y
#--------------------------------------------------------------------------#
n_sample = 5000
z_level = 3

set.seed(1)
u = rnorm(n = n_sample, 0, 1)
z_true = sample(x = 0:(z_level-1), replace = T, size = n_sample)

x_true = rbinom(n = n_sample,
                prob = exp(1 * u + 0.5 * z_true)/(1 + exp(1 * u + 0.5 * z_true)),
                size = 1)
y_true = rbinom(n = n_sample,
                prob = exp(3 * u + 5*x_true)/(1 + exp(3 * u + 5*x_true)),
                size = 1)

cor(z_true, y_true)
cor(z_true, x_true)

data_true = data.frame(z = z_true, x = x_true, y = y_true)

count_true = data_true %>% group_by(z, x, y) %>%
      summarise(n = n()) %>%
      as.data.frame()

print(xtable(count_true, caption = "Generation of Observed Data"), include.rownames=FALSE)

arm_size = data_true %>% group_by(z) %>%
      summarise(n = n()) %>%
      as.data.frame()

data_true = merge(count_true, arm_size, by = "z")

# true ACE
y_miss = rbinom(n = n_sample,
                prob = exp(3 * u + 5 *  (1-x_true))/(1 + exp(3 * u + 5 *  (1-x_true))),
                size = 1)

ace_true = (sum(y_true[x_true == 1]) + sum(y_miss[x_true==0]))/n_sample - (sum(y_true[x_true == 0]) + sum(y_miss[x_true==1]))/n_sample

#pyx.z_true = data_true[,4] / data_true[,5]

#ace.bounds(pyx.z_true)


# Specification of prior for p(x,y|z)
prior.z = rep(1, 4) # for each arm

## post
post.z = rep(prior.z, z_level) + count_true[, 4]

num.sims <- 10000  # Number of posterior simulations for each arm

theta.sims = t(replicate(num.sims, dirichlet(post.z[1:4])))
colnames(theta.sims) = c("X0,Y0|Z0","X0,Y1|Z0","X1,Y0|Z0","X1,Y1|Z0")
for(i in 1: (z_level-1)){
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

### posterior.theta.sims.iv contains simulations from the posterior
### distribution of p(x,y|z) under our model

n <- sum(all.ivs.ok)
n # number of sims remaining after removing those
# violating IV inequalities


#ace = apply(posterior.theta.sims.iv, 1, function(x) sum(x[(1:z_level) * 4])-sum(x[(0:(z_level-1)) * 4 + 2]))/4

#plot(density(ace))

ace.bound.mat = matrix(ncol = 2, nrow = num.sims)
for(k in 1: num.sims){
      ace.bound.mat[k,] = ace.bounds(posterior.theta.sims.iv[k,])
}
png("figure/effec-tri-plot.png")
do.tri.plot(ace.bound.mat, title.txt = "ACE Bounds", line = ace_true)
dev.off()

h <- hexbin(ace.bound.mat)

png("figure/effec-hexbin-plot.png")
plot(h, xlab = "Lower Bound on ACE", ylab = "Upper Bound on ACE")
dev.off()

