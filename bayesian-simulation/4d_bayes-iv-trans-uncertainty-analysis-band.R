rm(list = ls())

library(ggplot2)
library(tidyr)
library(dplyr)

#-------------------------------------------------------------------------------#
# In the following, we consider different methods to interpret the "confidence interval"
# 1. change (L_i, U_i) to (U_i-L_i, (U_i + L_i)/2)
# 2. Inference for ACE, pointwise caverage
# 3. Inference for ACE bounds, strong consistency
#-------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------#
# 1. change (L_i, U_i) to (U_i-L_i, (U_i + L_i)/2)
#-------------------------------------------------------------------------------#

priorz = 1
sample_size.vec = c(0, 100, 200, 400, 800, 1600)
load(file = "general_true_ace_and_ace-bound.Rdata")

mat = matrix(0, ncol = 7, nrow = length(sample_size.vec))
mat[,1] = c(1: length(sample_size.vec))
colnames(mat) = c("y", "L", "L_width", "M_min","M_max", "U_width", "U")

for(i in 1: length(sample_size.vec)){
        size_i = sample_size.vec[i]
        load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))
        
        mid.vec = (ace.bound.mat[,1] + ace.bound.mat[,2])/2
        width.vec = ace.bound.mat[,2] - ace.bound.mat[,1]
        
        
        # m_min = min(mid_quant)
        # m_max = max(mid_quant)
        m_min = quantile(mid.vec, 0.025)
        m_max = quantile(mid.vec, 0.925)
        
        w = median(width.vec)
        # find 95% upper bounds of the bound width
        w_quant = quantile(width.vec, 0.95)
        
        mat[i, "M_min"] = m_min
        mat[i, "M_max"] = m_max
           
        mat[i, "L_width"] = m_min - 0.5 * w
        mat[i, "U_width"] = m_max + 0.5 * w
        
        mat[i, "L"] = m_min - 0.5 * w_quant
        mat[i, "U"] = m_max + 0.5 * w_quant
}
mat = as.data.frame(mat)
mat$sample_size = sample_size.vec
mat_long = mat %>% gather(type, x, L:U, factor_key = T)
mat_long$sample_size = factor(mat_long$sample_size, levels = sample_size.vec)

uncertainty = rep("o.w.", nrow(mat_long))
uncertainty[mat_long$type == "L"] = "sampling uncertainty"
uncertainty[mat_long$type == "M_min"]= "sampling uncertainty"
uncertainty[mat_long$type == "U_width"] = "sampling uncertainty"
mat_long$uncertainty = uncertainty
#plot
ggplot(data = mat_long, aes(x = x, y = y, group = sample_size, color = uncertainty )) + geom_point(size = 4, color = "black", shape = "|") + geom_line(size = 2) + xlab("Interval") + ylab("Sample Size") + scale_y_continuous(breaks = c(1,2,3,4,5,6),labels = sample_size.vec) + theme_bw()
ggsave(filename = "bayes_iv_uncertainty_line_plot.pdf", width = 20, height = 10, units = "cm")


#-------------------------------------------------------------------------------#
# 2. Inference for ACE, pointwise caverage
#-------------------------------------------------------------------------------#
alpha = function(bounds, c){
        # input of c in [-4,4]
        # return alpha
        # it decreases in c
        lower = bounds[,1]
        upper = bounds[,2]
        
        tmp = min(pnorm(c) - pnorm(-c - (mean(upper) - mean(lower))/sd(upper)),
                  pnorm(c + (mean(upper) - mean(lower))/sd(lower)) - pnorm(-c))
        return(1 - tmp)
}


binarySearch = function(bounds, l=-5, r=5, x=0.05){
        while(l <= r){
                mid = (l + r)/2
                y = alpha(bounds, mid)
                #print(y)
                if(abs(y-x) < 10^-8){
                        return(mid)
                }else if(y < x){
                        r = mid
                }else{
                        l = mid
                }
                        
        }
        return("No result")
}



mat2 = matrix(0, ncol = 6, nrow = length(sample_size.vec))
mat2[,1] = c(1: length(sample_size.vec))
colnames(mat2) = c("y","mean_beta_l","mean_beta_u", "sd_beta_l",  "sd_beta_u", "c")

for(i in 1: length(sample_size.vec)){
        size_i = sample_size.vec[i]
        load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))
        mat2[i,-1 ] = c(apply(ace.bound.mat, 2, mean), 
                      apply(ace.bound.mat, 2, sd),
                      binarySearch(ace.bound.mat))
}

mat2 = as.data.frame(mat2)
var_l = mat2$mean_beta_l - mat2$c * mat2$sd_beta_l
var_u = mat2$mean_beta_u + mat2$c * mat2$sd_beta_u
mat2$var_l = var_l
mat2$var_u = var_u
mat2$sample_size = sample_size.vec
mat2_long = mat2 %>% gather(type, x, mean_beta_l, mean_beta_u,var_l, var_u, factor_key = T)
mat2_long$sample_size = factor(mat2_long$sample_size, levels = sample_size.vec)

uncertainty = rep("o.w.", nrow(mat2_long))
uncertainty[mat2_long$type == "var_l"] = "sampling uncertainty"
uncertainty[mat2_long$type == "mean_beta_u"]= "sampling uncertainty"
mat2_long$uncertainty = uncertainty

#plot
ggplot(mat2_long, aes(x = x, y = y, group = sample_size, color = uncertainty)) + geom_point(size = 4, color = "black", shape = "|") + geom_line(size = 2) + xlab("Interval") + ylab("Sample Size") + scale_y_continuous(breaks = c(1,2,3,4,5,6),labels = sample_size.vec) + theme_bw() 
ggsave(filename = "bayes_iv_uncertainty_pointwise.pdf", width = 20, height = 10, units = "cm")

#-------------------------------------------------------------------------------#
# 3. Inference for ACE bounds, strong consistency
#-------------------------------------------------------------------------------#
mat3 = matrix(0, ncol = 5, nrow = length(sample_size.vec))
mat3[,1] = c(1: length(sample_size.vec))
colnames(mat3) = c("y","mean_beta_l","mean_beta_u", "var_l",  "var_u")

for(i in 1: length(sample_size.vec)){
        size_i = sample_size.vec[i]
        load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))
        mean = apply(ace.bound.mat, 2, mean)
        sd = apply(ace.bound.mat, 2, sd)
        mat3[i,-1] = c(mean, mean[1] - qnorm(0.975)* sd[1], mean[2]+qnorm(0.975)* sd[2])
}

mat3 = as.data.frame(mat3)
mat3$sample_size = sample_size.vec
mat3_long = mat3 %>% gather(type, x, mean_beta_l:var_u, factor_key = T)
mat3_long$sample_size = factor(mat3_long$sample_size, levels = sample_size.vec)

uncertainty = rep("o.w.", nrow(mat3_long))
uncertainty[mat3_long$type == "var_l"] = "sampling uncertainty"
uncertainty[mat3_long$type == "mean_beta_u"]= "sampling uncertainty"
mat3_long$uncertainty = uncertainty

#plot
ggplot(mat3_long, aes(x = x, y = y, group = sample_size, color = uncertainty)) + geom_point(size = 4, color = "black", shape = "|") + geom_line(size = 2) + xlab("Interval") + ylab("Sample Size") + scale_y_continuous(breaks = c(1,2,3,4,5,6),labels = sample_size.vec) + theme_bw() 
ggsave(filename = "bayes_iv_uncertainty_ace_bounds_strong_consistency.pdf", width = 20, height = 10, units = "cm")


# make a new dataframe to combine the three results
mat_new = mat_long %>% select(x, y, sample_size, uncertainty)
mat2_new = mat2_long %>% select(x, y, sample_size, uncertainty)
mat3_new = mat3_long %>% select(x, y, sample_size, uncertainty)

mat_tol_long = rbind(mat_new, mat2_new, mat3_new)
mat_tol_long$methods = c(rep("transform", nrow(mat_new)),
                             rep("pointwise", nrow(mat2_new)),
                             rep("IR strong consistency", nrow(mat3_new)))

ggplot(mat_tol_long, aes(x = x, y = y, group = sample_size, color = uncertainty)) + 
        geom_point(size = 4, color = "black", shape = "|") + geom_line(size = 2) +
        facet_wrap(~methods, nrow = 3) + 
        scale_x_continuous(breaks = round(seq(min(mat_tol_long$x), max(mat_tol_long$x), by = 0.1),1)) + 
        theme_bw() + 
        theme(strip.text = element_text(face="bold", size=12,lineheight=5.0),
              strip.background = element_rect(fill="lightblue", colour="black",
                                              size=1),
              legend.text = element_text(size = 10, face = "bold"),
              legend.title = element_text(size = 12, face = "bold"),
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12, face = "bold"))+
        xlab("Interval") + ylab("Sample Size") + scale_y_continuous(breaks = c(1,2,3,4,5,6),labels = sample_size.vec) 
ggsave(filename = "bayes_iv_uncertainty_ace_overall.pdf", width = 20, height = 20, units = "cm")



#-------------------------------------------------------------------------------#
# 4. Inference for ACE, joint quantitle
#-------------------------------------------------------------------------------#
mat4 = matrix(0, ncol = 2, nrow = length(sample_size.vec))
mat4[,1] = c(1: length(sample_size.vec))
colnames(mat4) = c("Lower", "Upper")

for(i in 1: length(sample_size.vec)){
        size_i = sample_size.vec[i]
        load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))
        sd = apply(ace.bound.mat, 2, sd)
        mat4[i,] = tmp = quantile(ace.bound.mat, c(0.025, 0.975))
}
mat4 = as.data.frame(mat4)
mat4$sample_size = sample_size.vec
mat4_long$sample_size = factor(mat4_long$sample_size, levels = sample_size.vec)