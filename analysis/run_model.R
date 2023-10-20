
# script to do a simple model run

library(Rfast)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggforce)
library(patchwork)
library(truncnorm)
library(prodlim)

source(here::here("Model", "model.R"))
source(here::here("Model", "model_functions.R"))
theme_set(theme_bw(base_size = 11))

tsteps = 384
scenario = 1


###### Initial conditions
init = initial_piglet_setup()
bacteria = init$bacteria

## Grab final profile numbers from this setup
profile_end1.1 <- c(1,1,1,1,1,1,1,1,0,0,1) # pig - same as at start
profile_end2.1 <- c(0,1,0,0,1,0,0,0,1,0,2) # human - gains phi6 / phi2 / loses p4 (2/5/10)
profile_end2.2 <- c(0,1,0,0,1,0,1,0,1,0,2) # human
profiles_needed_end <- c(row.match(profile_end1.1, as.data.frame(bacteria)%>%select(-freq), nomatch = NA),
                         row.match(profile_end2.1, as.data.frame(bacteria)%>%select(-freq), nomatch = NA),
                         row.match(profile_end2.2, as.data.frame(bacteria)%>%select(-freq), nomatch = NA))

# bacteria[256,"freq"]=0 #removes parent 1
# bacteria[1793,"freq"]=0 #removes parent 2
difference_list = init$difference_list

##### Initial parameters 
parameters_in = define_parameters(scenario)

#### Run model 
results = piglet_mrsa_movement(10*tsteps, parameters_in, bacteria, difference_list)

# ggplot(results$totl_predict) +
#   geom_line(aes(time, total, colour = as.factor(parent))) +
#   theme_bw()
# 
# ggplot(results$prev_predict) +
#   geom_line(aes(time, value, colour = variable)) +
#   facet_grid(~parent) +
#   theme_bw()

#plot_circles(ini$bacteria, results$all_results)

###### compare to data 
### Data 
pigg_elements <- read.csv("data/pigg_elements.csv")[,-1]
totalsp <- read.csv("data/totals_bug.csv")[,-1]
### Reformat output
model_outputp <- results$prev_predict %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
model_outputp$prev <- round(model_outputp$value,2)
model_outputp <- rename(model_outputp, parent_strain = parent)
model_outputp$name <- recode(model_outputp$variable, V2 = "phi6",V5 = "phi2",V6 = "p1",V7 = "p2",V8 = "p3",V10 ="p4")
total_end <- unlist(results$all_results %>% filter(time == tsteps) %>% group_by(parent) %>% summarise(total = sum(value)) %>%select(total))
# Plot
g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")),
             aes(x=time, y = sum_prop, group = interaction(name, pig))) +
  geom_line(aes(col = name, linetype = factor(pig)),lwd = 1.5, alpha = 0.6) +
  geom_line(data =  model_outputp, 
            aes(x = time, y = prev, group = interaction(parent_strain,name)),linewidth = 1) +
  geom_point(aes(col = name),size = 1.5) +
  facet_wrap(name~parent_strain, ncol = 2) +
  scale_color_discrete("MGE") + 
  scale_linetype_discrete("Piglet") + 
  scale_y_continuous("Proportion of population", lim = c(0,1.1)) 

g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) +
  geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") +
  scale_y_log10("Total bacterial count") +
  scale_colour_discrete("Parent") + 
  geom_line(data = results$totl_predict, aes(x=time, y = total, group = parent, col = factor(parent)))

g3 <- ggplot(results$all_results %>% filter(variable %in% profiles_needed_end), aes(x=time, y = value)) + geom_line(aes(col = variable)) + geom_vline(xintercept = tsteps) + 
  geom_hline(yintercept = c(0.8, 0.2) * total_end) + scale_color_manual(values = c("red","green","blue"), breaks = c(profiles_needed_end), "Profile") + 
  scale_y_continuous("Number of bacteria")

g1 / (g2 + g3) + plot_layout(heights = c(3,1))
#ggsave("plots/example.jpeg")
