rm(list = ls())

source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")
source("check-iv-z-arbitrary.R")

filename = "transparent95_bounds_1004"
filename = "transparent95_bounds"


load(paste("bayesian-simulation/",filename, ".Rdata", sep = ""))
# load(file = "bayesian-simulation/true_ace_and_ace-bound.Rdata")
load(file = "bayesian-simulation/general_true_ace_and_ace-bound.Rdata")
ace.bound.mat = t(bound.mat)

na.index = which(ace.bound.mat[,1] == "NA")
num.na = length(na.index)

mat = apply(ace.bound.mat[-na.index, ], 2, as.numeric)

png(paste("bayesian-simulation/", filename, "true-ace.png", sep = ""))
do.tri.plot(mat, title.txt = "ACE Bounds", line = c(ace_true,ace_true))
dev.off()

# P(ACE ib ACE simulate bounds)
length(which(mat[,1] <= ace_true & mat[,2] >= ace_true))/nrow(mat)
# 0.7579909

# png("bayesian-simulation/general-transparent95-true-ace-bound.png")
# do.tri.plot(mat, title.txt = "ACE Bounds", line = ace_bound_true)
# dev.off()

t.test(as.numeric(mat[,1]), mu = ace_bound_true[1])
t.test(as.numeric(mat[,2]), mu = ace_bound_true[2])

