#### Scenario 1: All same loss / gain / fitness costs

###### Find odd behaviour

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

library(fmcmc) #https://uscbiostats.github.io/fmcmc/
library(adaptMCMC)

theme_set(theme_bw(base_size = 11))

#### CODE
source("code/piglet_mrsa_functions.R")

set.seed(42) # to get reproducible results

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
Initial.Values = rbind(c(mu = 0.2,
                   gamma = 0.4,
                   f = 0.001, 
                   grow = 2, 
                   rel_fit = 0.5), 
                   c(mu = 0.1,
                   gamma = 0.5,
                   f = -0.00001, 
                   grow = 1, 
                   rel_fit = 0.9), 
                   c(mu = 0.4,
                   gamma = 0.2,
                   f = -0.01, 
                   grow = 0.8, 
                   rel_fit = 1.2))

Initial.Values = c(mu = 0.2,
                         gamma = 0.4,
                         f = 0.001, 
                         grow = 2, 
                         rel_fit = 0.5)
# Check non-zero likelihood for start
#out <- piglet_mrsa_movement(tsteps, Initial.Values, ini$bacteria, ini$difference_list)
#run_sim_logPosterior(Initial.Values)

### read in 
#Initial.Values <- as.matrix(read.csv("fits/scn1_2022-02-05_00-26-30_GMT_trace.csv")[,-1])


### Try 
#khaario_scn1_1 <- kernel_adapt(Sd = .0000000000000001,freq = 1, warmup = 500, ub = c(rep(0.5,2),rep(0.5,1),3,1.5),
                             #lb = c(rep(0,2),rep(-0.5,1), rep(0,2)))

khaario_scn1_1 <- kernel_adapt(freq = 1, warmup = 500, ub = c(rep(0.5,2),rep(0.5,1),3,1.5),
                               lb = c(rep(0,2),rep(-0.5,1), rep(0,2)))

out_final1 <- fmcmc::MCMC(
  initial   = Initial.Values,                       # Automatically takes the last 2 points
  fun       = run_sim_logPosterior, 
  nsteps    = 1e3,                       # Increasing the sample size: about 6hrs for 5e3
  kernel    = khaario_scn1_1,
  thin      = 1
)

# Save output
filename = gsub(c(" "), "_", format(as.POSIXct(Sys.time()), tz = "Europe/London", usetz = TRUE))
filename = gsub(":", "-", filename)

write.csv(out_final1, here::here("fits/",paste0("scn1_a_",filename,"_","trace",".csv")))

# if multiple chains
#save(out_final1, file= here::here("fits/",paste0("scn1_a_",filename,"_","trace",".csv")))
#model = load("FILEPATH")

# Next 
khaario_scn1_2 <- khaario_scn1_1
out_final2 <- fmcmc::MCMC(
  initial   = out_final1,                       # Automatically takes the last 2 points
  fun       = run_sim_logPosterior_para, 
  nsteps    = 1e4,                       # Increasing the sample size: about 6hrs for 5e3
  kernel    = khaario_scn1_2, 
  thin      = 1
)

# Save output
filename = gsub(c(" "), "_", format(as.POSIXct(Sys.time()), tz = "Europe/London", usetz = TRUE))
filename = gsub(":", "-", filename)

save(out_final2, file= here::here("fits/",paste0("scn1_a_",filename,"_","trace",".csv")))

write.csv(out_final2, here::here("fits/",paste0("scn1_a_",filename,"_","trace",".csv")))

# #### Look at
library("lattice")  ## for xyplot
library("fitR")
plotESSBurn(out_final1) # Burn in = 1000?
plotESSBurn(out_final2) # Burn in = 100?

mcmc.trace.burned1 <- burnAndThin(out_final1, burn = 0)
mcmc.trace.burned2 <- burnAndThin(out_final2, burn = 1000)
mcmc.trace.burned3 <- burnAndThin(out_final3, burn = 1000)

acceptanceRate1 <- 1 - rejectionRate(mcmc.trace.burned1)
acceptanceRate2 <- 1 - rejectionRate(mcmc.trace.burned2)
acceptanceRate1
acceptanceRate2

plot(mcmc.trace.burned1)
plot(mcmc.trace.burned2)
plot(mcmc.trace.burned3)

autocorr.plot(mcmc.trace.burned1)
autocorr.plot(mcmc.trace.burned2)


xyplot(mcmc.trace.burned1)
xyplot(mcmc.trace.burned2)
effectiveSize(mcmc.trace.burned1) # aiming for 200 - 1000, <100 bad
effectiveSize(mcmc.trace.burned2) # First chain better

### Fits?
length(unique(mcmc.trace.burned1[,"mu"]))
plotPosteriorDensity(mcmc.trace.burned1)
plotPosteriorDensity(mcmc.trace.burned2)

levelplot(mcmc.trace.burned1, col.regions = heat.colors(100))
levelplot(mcmc.trace.burned2, col.regions = heat.colors(100))


library(foreach)
library(doParallel)
nc = detectCores()-3
cl = makeCluster(nc)
registerDoParallel(cl)

parasets <- as.data.frame(mcmc.trace.burned[!duplicated(mcmc.trace.burned[,"mu"]), ])
parasets$fit <- 0
for(i in 548:length(parasets[,1])){
  para <- parasets[i,1:(ncol(parasets)-1)]
  #parasets[i,"fit"] <- run_sim_logPosterior(para)
  write.csv(run_sim_logPosterior(para), paste0("fits/1fits_",i,".csv"))
}
stopCluster()


###
source("code/function_sample_posterior.R")
samples <- as.data.frame(posterior_sample_trace(as.matrix(mcmc.trace.burned), break_size = 0.01, plotit = TRUE, plotname = "1_1102")$samples)
nc = detectCores()-3
cl = makeCluster(nc)
registerDoParallel(cl)

samples$fit <- 0
for(i in 548:length(parasets[,1])){
  para <- parasets[i,1:(ncol(parasets)-1)]
  #parasets[i,"fit"] <- run_sim_logPosterior(para)
  write.csv(run_sim_logPosterior(para), paste0("fits/1fits_",i,".csv"))
}
stopCluster()
