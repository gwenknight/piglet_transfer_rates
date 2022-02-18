### Iterative lhs run 

library(data.table)
library(janitor)
library(foreach)
library(doParallel)
library(tidyverse)
library(lhs)
library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)
library(scales)
library(reshape2)
library(dplyr)
library(Rfast)
library(prodlim)
library(ggforce)
library(ggridges)
library(patchwork)
theme_set(theme_bw(base_size = 11))

source("code/iterative_lhs.R")
source("code/piglet_mrsa_functions.R")

### INITIAL CONDITON 1
Initial.Values = c(mu = 0.01,
                   gamma = 0.00000001,
                   f = 0.000001, 
                   grow = 0.17, 
                   rel_fit = 0.99)
# check not -Inf
# run_sim_logPosterior(Initial.Values) # -1122

Scen = "scn1"

#run LHS on this + / - 100%
m1 <- lhs_build_run(Initial.Values, limit = 0.5, paste0(Scen,"/iv1"), nsamples = 5)

# Run LHS on this + / - 100% 
new_iv = m1$max_ll_para
names(new_iv) <- names(Initial.Values)
m2 <- lhs_build_run(new_iv, 0.5, paste0(Scen,"/iv1_max"), nsamples = 1e4)

# Run LHS on this + / - 100% 
new_iv = m2$max_ll_para
names(new_iv) <- names(Initial.Values)
m3 <- lhs_build_run(new_iv, 0.5, paste0(Scen,"/iv1_max2"), nsamples = 1e4)

## Zoom in on max
all_worked = rbind(m1$worked,m2$worked,m3$worked)

max_para <- as.numeric(all_worked[which.max(all_worked[,1]),2:ncol(all_worked)])
names(max_para) <- names(Initial.Values)
# + / - 50%
mz_1 <- lhs_build_run(max_para, 2, paste0(Scen,"/ivz1"), nsamples = 1e4)
# + / - 30% 
new_iv = mz_1$max_ll_para
names(new_iv) <- names(Initial.Values)
mz_2 <- lhs_build_run(new_iv, 3, paste0(Scen,"/ivz2"), nsamples = 1e4)
# + / - 10% 
new_iv = mz_2$max_ll_para
names(new_iv) <- names(Initial.Values)
mz_3 <- lhs_build_run(new_iv, 10, paste0(Scen,"/ivz3"), nsamples = 1e4)

# Combine all of above
all_worked_with_zoom = rbind(all_worked, 
                             mz_1$worked,mz_2$worked,mz_3$worked)

#Take top 1000 likelihoods and look at parameter sets
para_gn <- all_worked_with_zoom %>% filter(!ll == -Inf)
top_ll <- para_gn[order(para_gn$ll,decreasing = TRUE)[1:1000],]

ggplot(top_ll, aes(x=value,group = name)) + geom_histogram() + facet_wrap(~name, scales = "free")

