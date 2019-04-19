rm(list = ls())
library(xtable)
library(ggplot2)
library(tidyr)

priorz = 0.25
sample_size.vec = c(0, 100, 200, 400, 800, 1600)
load(file = "general_true_ace_and_ace-bound.Rdata")

load(paste("bounds_sample_size_vary_prior_z", priorz,"_seed_17.Rdata", sep = ""))
bounds.mat = apply(bounds.mat, 2, as.numeric)
rownames(bounds.mat) = sample_size.vec
print(xtable(bounds.mat, caption = "Transparent 95% credible interval for ACE with different sample size", digits = 3))


x = c()
lower = c()
upper = c()
size = c()
for(i in 1: length(sample_size.vec)){
      size_i = sample_size.vec[i]
      load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))

      x = c(x, runif(n=nrow(ace.bound.mat), min = i-0.4, max = i+0.4))
      size = c(size, rep(size_i, nrow(ace.bound.mat)))
      lower = c(lower,ace.bound.mat[,1])
      upper = c(upper, ace.bound.mat[,2])
}

sim.mat = cbind(x, lower, upper, size)
sim.df= as.data.frame(rbind(sim.mat[,c(1,2,4)], sim.mat[,c(1,3,4)]))
names(sim.df) = c("x", "bounds", "size")
sim.df$bounds_lab = rep(c("lower", "upper"), each = nrow(sim.df)/2)

p = ggplot(data = sim.df, aes(x = x, y = bounds, color = bounds_lab)) + ylim(c(-1,1))
p.points = p + geom_point(size = 0.5)
       #geom_segment(aes(x = 0.6, xend = 1.4, y = bounds.mat[1,1], yend = bounds.mat[1,1]), size = 1.5, color = "black") + geom_segment(aes(x = 0.6, xend = 1.4, y = bounds.mat[1,2], yend = bounds.mat[1,2]), size = 1.5, color = "black") +
pp.points = p.points + geom_segment(aes(x = 1.6, xend = 2.4, y = bounds.mat[2,1], yend = bounds.mat[2,1]), size = 1.5, color = "black") + geom_segment(aes(x = 1.6, xend = 2.4, y = bounds.mat[2,2], yend = bounds.mat[2,2]), size = 1.5, color = "black") +
      geom_segment(aes(x = 2.6, xend = 3.4, y = bounds.mat[3,1], yend = bounds.mat[3,1]), size = 1.5, color = "black") + geom_segment(aes(x = 2.6, xend = 3.4, y = bounds.mat[3,2], yend = bounds.mat[3,2]), size = 1.5, color = "black") +
      geom_segment(aes(x = 3.6, xend = 4.4, y = bounds.mat[4,1], yend = bounds.mat[4,1]), size = 1.5, color = "black") + geom_segment(aes(x = 3.6, xend = 4.4, y = bounds.mat[4,2], yend = bounds.mat[4,2]), size = 1.5, color = "black") +
      geom_segment(aes(x = 4.6, xend = 5.4, y = bounds.mat[5,1], yend = bounds.mat[5,1]), size = 1.5, color = "black") + geom_segment(aes(x = 4.6, xend = 5.4, y = bounds.mat[5,2], yend = bounds.mat[5,2]), size = 1.5, color = "black") + 
  geom_segment(aes(x = 5.6, xend = 6.4, y = bounds.mat[6,1], yend = bounds.mat[6,1]), size = 1.5, color = "black") + geom_segment(aes(x = 5.6, xend = 6.4, y = bounds.mat[6,2], yend = bounds.mat[6,2]), size = 1.5, color = "black") 

pp.points + xlab("Sample Size") + scale_x_continuous(breaks = c(1,2,3,4,5,6),
                                                    labels = sample_size.vec)+
        geom_hline(aes(yintercept = ace_bound_true[1]), color = "red", linetype = "dashed", size = 1) +
        geom_hline(aes(yintercept = ace_bound_true[2]), color = "red", linetype = "dashed", size = 1) + 
      theme(legend.position="none")
ggsave(filename = "variability_sample_size_z025.png")
#ggsave(filename = "variability_sample_size.pdf")

# hexbin
# reshape the sim.df
sim_wide.df = sim.df %>% spread(bounds_lab, bounds)
hex.p = ggplot(sim_wide.df, aes(lower, upper)) + stat_binhex(colour="gray",na.rm=TRUE) + scale_fill_gradientn(colours=c("gray","black"),name = "Frequency",na.value=NA)+
        scale_x_continuous(limits = c(-1,1), breaks = seq(-1, 1, 0.25)) +
        scale_y_continuous(limits = c(-1,1), breaks = seq(-1, 1, 0.25))
hex.p + facet_wrap(~size, ncol =2) + 
        geom_vline(aes(xintercept = ace_bound_true[1]), color = "red", linetype = "dashed", size = 1) +
        geom_hline(aes(yintercept = ace_bound_true[2]), color = "red", linetype = "dashed", size = 1) + theme_bw()
ggsave(filename = "variability_sample_size_hexbin_z025_seed17.pdf", width = 20, height = 20, units = "cm")
