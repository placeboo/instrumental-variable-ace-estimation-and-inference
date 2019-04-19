rm(list = ls())

source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")
source("check-iv-z-arbitrary.R")

#load("bayesian-simulation/sim1000.Rdata")
load("bayesian-simulation/sim1000_0913.Rdata")
#load(file = "bayesian-simulation/true_ace_and_ace-bound.Rdata")
load(file = "bayesian-simulation/general_true_ace_and_ace-bound.Rdata")
ace.bound.mat = t(mat[c(1,2), ])


png("bayesian-simulation/general-posterior95-true-ace.png")
do.tri.plot(ace.bound.mat, title.txt = "ACE Bounds", line = c(ace_true,ace_true))
dev.off()


png("bayesian-simulation/general-posterior95-true-ace-bound.png")
do.tri.plot(ace.bound.mat, title.txt = "ACE Bounds", line = ace_bound_true)
dev.off()

t.test(ace.bound.mat[,1], mu = ace_bound_true[1])
t.test(ace.bound.mat[,2], mu = ace_bound_true[2])

# P(ACE bound \subset ACE simulate bounds)
length(which(ace.bound.mat[,1] <= ace_bound_true[1] & ace.bound.mat[,2] >= ace_bound_true[2]))/nrow(ace.bound.mat) # 0.965

plot(density(mat[5,]))
length(which(mat[5,] >= 0.95))
