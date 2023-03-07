
# script to do a simple model run

library(Rfast)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggforce)
library(patchwork)
library(truncnorm)

source(here::here("Model", "model.R"))
source(here::here("Model", "model_functions.R"))

tsteps = 384
scenario = 3

init = initial_piglet_setup()
parameters_in = define_parameters(scenario)

bacteria = init$bacteria
# bacteria[256,"freq"]=0 #removes parent 1
# bacteria[1793,"freq"]=0 #removes parent 2
difference_list = init$difference_list

results = piglet_mrsa_movement(tsteps, parameters_in, bacteria, difference_list)

ggplot(results$totl_predict) +
  geom_line(aes(time, total, colour = as.factor(parent))) +
  theme_bw()

ggplot(results$prev_predict) +
  geom_line(aes(time, value, colour = variable)) +
  facet_grid(~parent) +
  theme_bw()

#plot_circles(ini$bacteria, results$all_results)

