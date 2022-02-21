### Run lhs on big range for all lhs

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

# Scen1
Initial.Values = c(mu = 0.01,
                   gamma = 0.00000001,
                   f = 0.000001, 
                   grow = 0.17, 
                   rel_fit = 0.99)
limit1 = cbind(c(rep(0,2),rep(-0.5,1), rep(0,2)),c(rep(0.5,2),rep(0.5,1),3,1.5))
m1 <- lhs_build_run(Initial.Values, limit = limit1, "lhs_all/sc1", nsamples = 5)

# Scen2
Initial.Values = c(mu_phage = 0.01, mu_plasmid = 0.01, 
                   gamma_phage = 0.00000001,gamma_plasmid = 0.00000001,
                   f_phage = 0.000001, f_plasmid = 0.00000001,
                   grow = 0.17, 
                   rel_fit = 0.99)
limit2 = cbind(c(rep(0,4),rep(-1,2), rep(0,2)),c(rep(1,4),rep(1,2),3,1.7))
m2 <- lhs_build_run(Initial.Values, limit = limit2, "lhs_all/sc2", nsamples = 5)

# Scen3
Initial.Values = c(mu = 0.01,
                   gamma = 0.00000001,
                   f2 = 0.000001, f5 = 0.000001, f6 = 0.000001, 
                   f7 = 0.000001, f8 = 0.000001, f10 = 0.000001, 
                   grow = 0.17, 
                   rel_fit = 0.99)
limit3 = cbind(c(rep(0,2),rep(-0.7,6), rep(0,2)), c(rep(0.7,2),rep(0.7,6),3,1.5))
m3 <- lhs_build_run(Initial.Values, limit = limit3, "lhs_all/sc3", nsamples = 5)

# Scen4
Initial.Values = c(mu2 = 0.02, mu5 = 0.00005, mu6 = 0.00005,
                   mu7 = 0.00005, mu8 = 0.00005, mu10 = 0.00005,
                   gamma2 = 0.00005, gamma5 = 0.00005,
                   gamma6 = 0.00005, gamma7 = 0.00005,
                   gamma8 = 0.0000508, gamma10 = 0.00005,
                   f2 = -0.00005, f5 = -0.00005, f6 = 0.00005,
                   f7 = -0.00005, f8 = 0.00005, f10 = 0.00005,
                   grow = 0.12, rel_fit = 1.08)
limit4 = cbind(c(rep(0,12),rep(-1,6), rep(0,2)), c(rep(1,12),rep(1,6),3,1.5))
m4 <- lhs_build_run(Initial.Values, limit = limit4, "lhs_all/sc4", nsamples = 5)

## No loss 
# Scen4
Initial.Values = c(mu2 = 0, mu5 = 0, mu6 = 0,
                   mu7 = 0, mu8 = 0, mu10 = 0,
                   gamma2 = 0.00005, gamma5 = 0.00005,
                   gamma6 = 0.00005, gamma7 = 0.00005,
                   gamma8 = 0.0000508, gamma10 = 0.00005,
                   f2 = -0.00005, f5 = -0.00005, f6 = 0.00005,
                   f7 = -0.00005, f8 = 0.00005, f10 = 0.00005,
                   grow = 0.12, rel_fit = 1.08)
limit4_nl = cbind(c(rep(0,12),rep(-1,6), rep(0,2)), c(rep(0,6),rep(1,6),rep(1,6),3,1.5))
m4_nl <- lhs_build_run(Initial.Values, limit = limit4_nl, "lhs_all/sc4_nl", nsamples = 5)