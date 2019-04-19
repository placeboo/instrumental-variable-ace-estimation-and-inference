library(dplyr)
library(ggplot2)
library(tidyr)
rm(list = ls())
load(file = "general_true_ace_and_ace-bound.Rdata")

priorz = 1
sample_size.vec = c(0, 100, 200, 400, 800, 1600)

mat = matrix(0, ncol = 7, nrow = length(sample_size.vec))
mat[,1] = c(1: length(sample_size.vec))
colnames(mat) = c("L_0025", "L_005","L_max",
                  "U_min","U_095","U_0975",
                  "N")

for(i in 1: length(sample_size.vec)){
        size_i = sample_size.vec[i]
        load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))
        mat[i, c(1:3)] = quantile(ace.bound.mat[,1], c(0.025, 0.05, 1))
        mat[i, c(4:6)] = quantile(ace.bound.mat[,2], c(0, 0.95, 0.975))
        mat[i, 7] = size_i

}
mat = as.data.frame(mat)
mat$y = 1: length(sample_size.vec)
# prepara dataframe for plot

mat.wide = mat %>%
        dplyr::select(L_005, U_095, y, N) %>%
        filter(N != 0) %>%
        mutate(length = U_095 - L_005)
mat.long = mat.wide %>% gather(key = ci, x, L_005, U_095)

mat.long %>% ggplot(aes(x = x, y = y, group = N)) +
        geom_point(size = 4, color = "black", shape = "|") +
        geom_line(size = 2) +
        xlab("Interval") + ylab("Sample Size") +
        scale_y_continuous(breaks = c(2,3,4,5,6),
                           labels = sample_size.vec[-1]) +
        geom_vline(aes(xintercept = ace_bound_true[1]), color = "red", linetype = "dashed", size = 1) +
        geom_vline(aes(xintercept = ace_bound_true[2]), color = "red", linetype = "dashed", size = 1) +
        theme_bw()
