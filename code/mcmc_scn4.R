#### Scenario 4: all different

### Look at output all same

library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)
library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(Rfast)
library(patchwork)
library(prodlim)
source("code/piglet_mrsa_functions.R")
theme_set(theme_bw(base_size = 11))
set.seed(42) # to get reproducible results

library(fmcmc) #https://uscbiostats.github.io/fmcmc/
library(coda)
library(adaptMCMC)
### Data
data <- read.csv("data/data_to_fit.csv")[,-1]
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
### Likelihood
dist_like <- read.csv("data/seen_predicted.csv")[,-1]
pigg_elements <- read.csv("data/pigg_elements.csv")[,-1]
totalsp <- read.csv("data/totals_bug.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>%
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

ini <- initial_piglet_setup(384)

# parameters
Initial.Values = c(mu2 = 2.58492135046923e-05, mu5 = 1.23694298645124e-07, mu6 = 0.0091543430464811, 
                   mu7 = 5.63000577873516e-06, mu8 = 0.0100363784007989, mu10 = 0.00644058325472071, 
                   gamma2 = 1.03260990391265e-05, gamma5 = 1.29510851107398e-06, 
                   gamma6 = 7.43707130531284e-12, gamma7 = 9.90849379073332e-07, 
                   gamma8 = 1.18084441913925e-08, gamma10 = 1.25138645913259e-12, 
                   f2 = -0.222792921010004, f5 = -0.311241161192935, f6 = 0.635171096823821, 
                   f7 = -0.172511628201815, f8 = 0.281506548980894, f10 = 0.781415083339214, 
                   grow = 0.125465546487709, rel_fit = 1.08138114848714)

# breaks: was dividing by zero 
# Initial.Values = c(mu2 = 0.0163, mu5 = 0.0147, mu6 = 0.0028, 
#                    mu7 = 0.0113, mu8 = 0.0503, mu10 = 0.0144, 
#                    gamma2 = 0.0088, gamma5 = 0.0294, 
#                    gamma6 = 0.0227, gamma7 = 0.0102, 
#                    gamma8 = 0.0198, gamma10 = 0.0027, 
#                    f2 = -0.43, f5 = -0.23, f6 = -0.24, 
#                    f7 = 0.289, f8 = -0.36, f10 = -0.19, 
#                    grow = 2.5, rel_fit = 0.23)

# Check non-zero likelihood for start
#out <- piglet_mrsa_movement(tsteps, Initial.Values, ini$bacteria, ini$difference_list)
run_sim_logPosterior(Initial.Values)


# # # out <- fmcmc::MCMC(
# # #   initial = Initial.Values,
# # #   fun     = run_sim_logPosterior,
# # #   nsteps  = 1e3,
# # #   kernel  = kernel_normal(scale = .0000001) 
# # # )
# # # plot(out[,1:3])
# # 
# # 
# # out <- fmcmc::MCMC(
# #   initial = out,
# #   fun     = run_sim_logPosterior,
# #   nsteps  = 5e3,
# #   kernel  = kernel_normal(scale = .0000001)
# # )
# 
# 
# out <- fmcmc::MCMC(
#   initial = out,
#   fun     = run_sim_logPosterior,
#   nsteps  = 5e3,
#   kernel  = kernel_normal(scale = .0000001)
# )

matrix_corr <- read.csv("fits/scn4_cov_lhs3.csv")[,-1]
sd <- read.csv("fits/scn4_sd_lhs3.csv")[,-1]
# try adapt
khaario <- kernel_adapt(Sd = sd,freq = 1, warmup = 500, ub = c(rep(0.1,12),rep(0.5,6),3,1.5),
                        lb = c(rep(0,12),rep(-0.5,6), rep(0,2)))

set.seed(12) 

out_haario_1 <- fmcmc::MCMC(
  initial   = Initial.Values,                       
  fun       = run_sim_logPosterior, 
  nsteps    = 5000,    
  kernel    = khaario, # We passed the predefined kernel
  thin      = 1,       # No thining here
  nchains   = 1L,      # A single chain
  multicore = FALSE    # Running in serial
)
plot(out_haario_1[,1:5])

khaario <- kernel_adapt(Sd = .00000001, freq = 1, warmup = 500, lb = 0)

out_haario_2 <- fmcmc::MCMC(
  initial   = out_haario_1,                       
  fun       = run_sim_logPosterior, 
  nsteps    = 3000,    
  kernel    = khaario, # We passed the predefined kernel
  thin      = 1,       # No thining here
  nchains   = 1L,      # A single chain
  multicore = FALSE    # Running in serial
)
plot(out_haario_2[,1:5])
write.csv(out_haario_2, "fits/3101_scn2_3e3.csv")


###### GO BACK to mfidd
library("fitR")
source("code/mcmcmh.r")
matrix_corr <- read_csv("fits/scn4_cov_lhs3.csv")[,-1]
colnames(matrix_corr) <- names(Initial.Values)
rownames(matrix_corr) <- names(Initial.Values)
matrix_corr <- as.matrix(matrix_corr)

mcmc.epi3 <- mcmcMH(target = run_sim_logPosterior,
                       init.theta = Initial.Values,
                    limits = list(lower = c(rep(0,12),rep(-0.5,6), rep(0,2)), 
                                  upper = c(rep(0.1,12),rep(0.5,6),3,1.5)),
                       proposal.sd = sd,
                    covmat = matrix_corr,
                       n.iterations = 1000)

