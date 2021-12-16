##### Fitting
library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)
library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(Rfast)


### Data 
data <- read.csv("data/data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))

### Likelihood 
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

### Functions
source("code/piglet_mrsa_functions.R")
source("code/mcmcmh.r") # had to change line172

## Phage same, plasmids same 
# Morning of 15/12 - only 1,000. 
# init.theta = c(mu2 = 0.09, mu5 = 0.096, mu6 = 0.106, mu7 = 0.032, mu8 = 0.053, 
#              mu10 = 0.037, gamma2 = 0.015, gamma5 = 0.011, gamma6 = 0.05, 
#              gamma7 = 0.011, gamma8 = 0.027, gamma10 = 0.059, f2 = 0.018, 
#              f5 = 0.007, f6 = 0.002, f7 = 0.015, f8 = 0.012, f10 = 0.011, 
#              grow = 0.1)

# Last trace from previous
init.theta = c(mu2 = 0.188878541032485, mu5 = 0.191811127529857, mu6 = 1.06342477430845,
  mu7 = 0.650240103957843, mu8 = 1.13092685382576, mu10 = 1.36327733057528,
  gamma2 = 0.299144953920434, gamma5 = 1.10429707130617, gamma6 = 1.15332183481689,
  gamma7 = 0.586259749844624, gamma8 = 1.42781770425473, gamma10 = 0.669795153200574,
  f2 = 0.0752560231516298, f5 = 0.149273263311754, f6 = 0.104081458666405,
  f7 = 0.00298936941316733, f8 = 0.29765944219356, f10 = 0.30328253933773,
  grow = 1.72361721520587)
# # From 3,000 git early 16th Dec
# c(mu2 = 0.24236127162781, mu5 = 0.291186237125602, mu6 = 1.39972101446994, 
#   mu7 = 0.475071870965611, mu8 = 2.11453825138162, mu10 = 1.27364080641482, 
#   gamma2 = 0.923926481290458, gamma5 = 0.288443809560827, gamma6 = 2.49297498091131, 
#   gamma7 = 1.50374708470661, gamma8 = 0.655867138827866, gamma10 = 1.58363153572681, 
#   f2 = 0.105630705638364, f5 = 0.0773147543120973, f6 = 0.0812436946234052, 
#   f7 = 0.0847738803512515, f8 = 0.0304607459475252, f10 = 0.548471037221678, 
#   grow = 2.35978774108334)
# From 10,000 here 16th Dec
c(mu2 = 0.182572898539728, mu5 = 0.20403048887804, mu6 = 2.71946067053579, 
  mu7 = 0.647146010768494, mu8 = 3.15087842581532, mu10 = 2.43570506680262, 
  gamma2 = 0.250529800725404, gamma5 = 0.187362977110421, gamma6 = 3.96217620667374, 
  gamma7 = 0.882019151324511, gamma8 = 0.428272435768007, gamma10 = 2.39253967826668, 
  f2 = 0.00982410671178733, f5 = 0.187112702210853, f6 = 0.0966775914439774, 
  f7 = 0.0292466113278493, f8 = 0.0762295307859667, f10 = 0.452847036866703, 
  grow = 2.16900158786653)

lower.p <- init.theta
lower.p[] <- 0

mcmc.epi3_79 <- mcmcMH(target = run_sim_logPosterior,
                       limits=list(lower = lower.p),
                       init.theta = c(mu2 = 0.182572898539728, mu5 = 0.20403048887804, mu6 = 2.71946067053579, 
                                      mu7 = 0.647146010768494, mu8 = 3.15087842581532, mu10 = 2.43570506680262, 
                                      gamma2 = 0.250529800725404, gamma5 = 0.187362977110421, gamma6 = 3.96217620667374, 
                                      gamma7 = 0.882019151324511, gamma8 = 0.428272435768007, gamma10 = 2.39253967826668, 
                                      f2 = 0.00982410671178733, f5 = 0.187112702210853, f6 = 0.0966775914439774, 
                                      f7 = 0.0292466113278493, f8 = 0.0762295307859667, f10 = 0.452847036866703, 
                                      grow = 2.16900158786653),
                       proposal.sd = c(rep(0.005,19)),
                       n.iterations = 2000,
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
