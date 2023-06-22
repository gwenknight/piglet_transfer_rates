
# script to run mcmc

library(BayesianTools)
library(Rfast)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggforce)
library(patchwork)
library(truncnorm)
#library(lognorm)

source(here::here("model", "model.R"))
source(here::here("model", "model_functions.R"))
source(here::here("mcmc", "mcmc.R"))
source(here::here("mcmc", "createPrior_gk.R"))

dist_like = read.csv(here::here("data", "seen_predicted.csv"))[,-1]
data = read.csv(here::here("data", "data_to_fit.csv"))[,-1]
data_6 = data %>%
  filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
totals = read.csv(here::here("data", "totals_bug.csv"))[,-1] %>%
  select(time,value, parent, lim) %>%
  pivot_wider(names_from = lim) %>%
  mutate(weight = 1/(max - min))

init = initial_piglet_setup()
bacteria = init$bacteria
difference_list = init$difference_list

profs = apply(bacteria[,-11],1, function(x) paste0(x, collapse = ""))
profile_end1.1 = paste0(c(1,1,1,1,1,1,1,1,0,0,1), collapse = "") # pig - same as at start
profile_end2.1 = paste0(c(0,1,0,0,1,0,0,0,1,0,2), collapse = "") # human - gains phi6 / phi2 / loses p4 (2/5/10)
profile_end2.2 = paste0(c(0,1,0,0,1,0,1,0,1,0,2), collapse = "") # human
profiles_needed_end = c(which(profs == profile_end1.1),
                        which(profs == profile_end2.1),
                        which(profs == profile_end2.2))


## Define parameters based on scenario ###############

scenario = 1
parameters_in = define_parameters(scenario)


## Create priors based on chosen scenario ###############
priors = define_priors(parameters_in)

# Error later to do with names of parameters - priors$density needs parameters_in to have names
priors$density(parameters_in)
priors$density(c(0.0001,0.0002,0.0002,0.02,0)) # doesn't work without names

priors$sampler()
priors$sampler(1)  # generates samples with names now (changed for n = 5 case)
priors$sampler(10)  # generates samples with names now (changed for n = 5 case)

noname_parameters_in <- unname(parameters_in)
priors_nn = define_priors(noname_parameters_in)
priors_nn$density(noname_parameters_in) # but priors built with names in 

## Run the sampler #####################

bayesianSetup = createBayesianSetup(run_sim_logPosterior, priors) #, names = names(parameters_in)) # I think names is for the settings? not the parameters? 

#settings for the sampler, see documentation for more info
settings = list(iterations = 100)

#run the mcmc
out = runMCMC(bayesianSetup = bayesianSetup, settings = settings, sampler = "DE")

## TODO:
# - optimise likelihood function, currently too slow
# - replace names with indices to access parameter vector

summary(out)
plot(out)
getSample(out)
