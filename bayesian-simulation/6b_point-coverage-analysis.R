rm(list = ls())

library(ggplot2)
library(tidyr)
library(HDInterval)

load(file = "ace_point_coverage_2019-01-31.RData")

# compare the lower quantile of L and upper quantile of U
priorz = 1
sample_size.vec = c(0, 100, 200, 400, 800, 1600)
ci = 0.95
quant.mat = matrix(NA, ncol = 2, nrow = length(sample_size.vec))
for (i in 1 : length(sample_size.vec)){

        size_i = sample_size.vec[i]
        load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))

        quant.mat[i, ] = c(hdi(ace.bound.mat[,1], ci)[1], hdi(ace.bound.mat[,2], ci)[2])
        #quant.mat[i, ] = c(quantile(ace.bound.mat[,1], 0.05), quantile(ace.bound.mat[,2], 0.95))
}

ci.mat = cbind(ci.mat, quant.mat)
colnames(ci.mat) = c("L_star","U_star","L0025", "U0975")

# visualization
# sample size 0

# load(paste("bounds_sample_size_", 0,"priorz_", priorz,".Rdata", sep = ""))
# ace.bound.df = as.data.frame(ace.bound.mat)
# names(ace.bound.df) = c("L", "U")
# ace.bound.long = ace.bound.df %>%
#         gather(bounds, x, L, U)
#
# ace.bound.long %>% ggplot(aes(x = x, color = bounds)) +
#         geom_density() +
#         geom_vline(aes(xintercept = quantile(ace.bound.mat[,1], 0.025)), size = 1) +
#         geom_vline(aes(xintercept = quantile(ace.bound.mat[,2], 0.975)), size = 1) +
#         geom_vline(aes(xintercept = ci.mat[mat[,1] %in% 0, 2]), color = "red", linetype = "dashed", size = 1) +
#         geom_vline(aes(xintercept = ci.mat[mat[,1] %in% 0, 3]), color = "red", linetype = "dashed", size = 1) +
#         theme_bw()

drawDensity = function(sample.size, priorz = 1){
        load(paste("bounds_sample_size_", sample.size,"priorz_", priorz,".Rdata", sep = ""))
        ace.bound.df = as.data.frame(ace.bound.mat)
        names(ace.bound.df) = c("L", "U")
        ace.bound.long = ace.bound.df %>%
                gather(bounds, x, L, U)
        ace.bound.long %>% ggplot(aes(x = x, color = bounds)) +
                geom_density() +
                geom_vline(aes(xintercept = hdi(ace.bound.mat[,1], 0.95)[1]), size = 1) +
                geom_vline(aes(xintercept = hdi(ace.bound.mat[,2], 0.95)[2]), size = 1) +
                geom_vline(aes(xintercept = ci.mat[sample_size.vec %in% sample.size, 2]), color = "red", linetype = "dashed", size = 1) +
                geom_vline(aes(xintercept = ci.mat[sample_size.vec %in% sample.size, 3]), color = "red", linetype = "dashed", size = 1) +
                theme_bw()
}

drawDensity(0)
drawDensity(800)
drawDensity(200)
drawDensity(100)
# display ci according different sample size
load(file = "general_true_ace_and_ace-bound.Rdata")

ci.df = as.data.frame(ci.mat)
ci.df$sample_size = sample_size.vec
ci.df$y = 1: length(sample_size.vec)
ci.long = ci.df %>% gather(key = CI, x, L_star:U_star)

ci.long %>% ggplot(aes(x = x, y = y, group = sample_size)) +
        geom_point(size = 4, color = "black", shape = "|") +
        geom_line(size = 2) +
        xlab("Interval") + ylab("Sample Size") +
        scale_y_continuous(breaks = c(1:6),
                           labels = sample_size.vec) +
        geom_vline(aes(xintercept = ace_bound_true[1]), color = "red", linetype = "dashed", size = 1) +
        geom_vline(aes(xintercept = ace_bound_true[2]), color = "red", linetype = "dashed", size = 1) +
        theme_bw()

