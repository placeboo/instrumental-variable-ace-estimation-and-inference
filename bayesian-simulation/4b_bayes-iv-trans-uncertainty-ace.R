rm(list = ls())
library(ggplot2)
library(MCMCpack)

sample_size.vec = c(0, 100, 200, 400, 800, 1600)
load(file = "general_true_ace_and_ace-bound.Rdata")

size = c()
ace = c()


for(i in 1: length(sample_size.vec)){
  size_i = sample_size.vec[i]
 
  load(paste("bounds_sample_size_", size_i, "priorz_1",".Rdata", sep = ""))
  size = c(size, rep(size_i, nrow(ace.bound.mat)))
  ace = c(ace, apply(ace.bound.mat, 1, function(x) runif(1, min = x[1], max = x[2])))
}

ace.df = data.frame(ace=ace, size = as.factor(size))

ace.p = ggplot(data = ace.df, aes(ace)) + geom_density() + scale_x_continuous(breaks = seq(-1, 1, 0.2)) + geom_vline(aes(xintercept = ace_true), color = "red", linetype = "dashed", size = 1)
ace.p + facet_wrap(~size, ncol =2)
ggsave(filename = "bayes-iv-tran-ace-prior_z1.pdf", width = 20, height = 15, units = "cm")
