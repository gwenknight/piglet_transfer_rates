##### Fitting
library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)


### Data 
data <- read.csv("data/data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))

### Likelihood 
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

### Functions
source("code/functions.R")
source("code/mcmcmh.r") # had to change line172

## Phage same, plasmids same 
init.theta = c(mu2 = 0.09, mu5 = 0.096, mu6 = 0.106, mu7 = 0.032, mu8 = 0.053, 
             mu10 = 0.037, gamma2 = 0.015, gamma5 = 0.011, gamma6 = 0.05, 
             gamma7 = 0.011, gamma8 = 0.027, gamma10 = 0.059, f2 = 0.018, 
             f5 = 0.007, f6 = 0.002, f7 = 0.015, f8 = 0.012, f10 = 0.011, 
             grow = 0.1)

lower.p <- init.theta
lower.p[] <- 0

mcmc.epi3_79 <- mcmcMH(target = run_sim_logPosterior,
                       limits=list(lower = lower.p),
                       init.theta = c(mu2 = 0.09, mu5 = 0.096, mu6 = 0.106, mu7 = 0.032, mu8 = 0.053, 
                                      mu10 = 0.037, gamma2 = 0.015, gamma5 = 0.011, gamma6 = 0.05, 
                                      gamma7 = 0.011, gamma8 = 0.027, gamma10 = 0.059, f2 = 0.018, 
                                      f5 = 0.007, f6 = 0.002, f7 = 0.015, f8 = 0.012, f10 = 0.011, 
                                      grow = 0.1),
                       proposal.sd = c(rep(0.005,19)),
                       n.iterations = 100,
                       adapt.size.start = 100,
                       adapt.shape.start = 500,
                       adapt.size.cooling=0.999, 
                       verbose = TRUE)

# Save output
filename = gsub(c(" "), "_", format(as.POSIXct(Sys.time()), tz = "Europe/London", usetz = TRUE))
filename = gsub(":", "-", filename)

write.csv(mcmc.epi3_79$trace, here::here("fits/everything/",paste0(filename,"_","trace",".csv")))
write.csv(mcmc.epi3_79$acceptance.rate, here::here("fits/everything/",paste0(filename,"_","acceptance_rates",".csv")))
write.csv(mcmc.epi3_79$covmat.empirical, here::here("fits/everything/",paste0(filename,"_","covmat",".csv")))


# # Look at output
# mcmc.trace <- mcmc(mcmc.epi3_79$trace)
# summary(mcmc.trace)
# acceptanceRate <- 1 - rejectionRate(mcmc.trace)
# acceptanceRate
# plot(mcmc.trace)
