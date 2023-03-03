
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

Initial.Values4 = c(mu2 = 0, mu5 = 0, mu6 = 0, 
                    mu7 = 0.65, mu8 = 0, mu10 = 0.1, 
                    gamma2 = 0.01, gamma5 = 0.01, 
                    gamma6 = 0.000000000001, gamma7 = 0.0005, 
                    gamma8 = 0.000000000001, gamma10 = 0.00000001, 
                    f2 = -0.4, f5 = -0.4, f6 = 0.5, 
                    f7 = -0.5, f8 = 0.4, f10 = 0.4, 
                    grow = 0.08, rel_fit = 0.95)
Initial.Values3 = c(mu = 0.01,
                    gamma = 0.00000001,
                    f2 = 0.000001, f5 = 0.000001, f6 = 0.000001, 
                    f7 = 0.000001, f8 = 0.000001, f10 = 0.000001, 
                    grow = 0.17, 
                    rel_fit = 0.99)
Initial.Values2 = c(mu_phage = 0.01, mu_plasmid = 0.01, 
                    gamma_phage = 0.00000001,gamma_plasmid = 0.00000001,
                    f_phage = 0.000001, f_plasmid = 0.00000001,
                    grow = 0.17, 
                    rel_fit = 0.99)
Initial.Values1 = c(mu = 0.1,
                    gamma = 0.000000000001,
                    f = 0.1, 
                    grow = 0.08, 
                    rel_fit = 0.95)

tsteps = 384
ini = initial_piglet_setup(tsteps)

parameters_in = Initial.Values3
bacteria = ini$bacteria
# bacteria[256,"freq"]=0 #removes parent 1
# bacteria[1793,"freq"]=0 #removes parent 2
difference_list = ini$difference_list

results = piglet_mrsa_movement(tsteps, parameters_in, bacteria, difference_list)

ggplot(results$totl_predict) +
  geom_line(aes(time, total, colour = as.factor(parent))) +
  theme_bw()

ggplot(results$prev_predict) +
  geom_line(aes(time, value, colour = variable)) +
  facet_grid(~parent) +
  theme_bw()

plot_circles(ini$bacteria, results$all_results)

