###### Convert data to fitting data

### read in distributions from back_calc_prev.R
seen_predicted <- read.csv("seen_predicted.csv")[,-1]
### read in data on elements
pigg_elements <- read.csv("pigg_elements.csv")[,-1]

### Filter distributions 
conversion_data = seen_predicted %>% #filter(prob > 0.05) %>% # Want to only have >1% probability
  select(prev, n_colonies, prob) %>% mutate(prev_in_pig = prev/100,
                                            prev_in_sample = n_colonies / 20) %>%
  select(prev_in_pig, prev_in_sample, prob) # This is the probability see y/20 with this element if model has x 

### Piglets
## Need to average over the 4 piglets 
## P (see prev in sample | prev in model) = P(see prev in samples from 4 piglets | prev in model pig)

element_data = pigg_elements %>% mutate(prev_in_sample = round(sum_prop,2)) %>% select(-sum_prop) # This has the prevalence in each pig
    
ggplot(element_data, aes(x=pig, y = prev_in_sample)) + geom_point(aes(col = factor(time))) + 
  facet_wrap(parent_strain ~ name)

element_data_prob = pigg_elements %>% mutate(prev_in_sample = round(sum_prop,2)) %>% select(-sum_prop) %>% # This has the prevalence in each pig
  group_by(parent_strain, time, name) %>% count(prev_in_sample) %>% mutate(prob_over_piglets = n / 4)

ggplot(element_data_prob, aes(x=prev_in_sample, y = prob_over_piglets)) + geom_point(aes(col = factor(time))) + 
  geom_line(aes(col = factor(time))) + 
  scale_y_continuous(lim = c(0,1)) + 
  facet_wrap(parent_strain ~ name)

#### Fit distributions to this??? otherwise lose the variation across the 4 piglets? 

### or just take mean... 
element_data_prob_mean = element_data_prob %>% mutate(mean_pig = prob_over_piglets * prev_in_sample) %>% group_by(parent_strain, time, name) %>% summarise(mean = sum(mean_pig))
ggplot(element_data_prob, aes(x=time, y = mean)) + geom_point(aes(col = factor(time))) + 
  scale_y_continuous(lim = c(0,1)) + 
  facet_wrap(parent_strain ~ name)

elemtn <- element_data_prob_mean %>% select(parent_strain, time, name, mean) %>%
  rename(n_colonies_prev = mean, parent = parent_strain) %>%
  filter(time > 0)

write.csv(elemtn,"data_to_fit.csv")



##### DO weighted mean to give likelihood
### No can't as the data is not in the piglet - it is colony counts. So just use the above mean! 

# 
# data_fit = left_join(element_data, conversion_data, by = "prev_in_sample") %>%
#   #  select(-prev_in_pig) %>%
#   group_by(parent_strain, time, name, prev_in_pig)
# 
# # Should add up to 1
# tail(data_fit %>% group_by(pig, prev_in_sample, parent_strain, time, name) %>% summarise(sum(prob)))
# 
# ## data_fit has 
# # prev_in_pig = what the data says
# # prev_in_sample = if see this prevalence = prob it is true
# 
# 
# # Check what adding together for phi6
# theme_set(theme_bw(base_size = 6))
# ggplot(data_fit %>% filter(name == "phi6"), aes(x=prev_in_pig, y = factor(time), group = time, 
#                                                 height = prob)) +
#   geom_ridgeline(aes(fill = time)) +
#   facet_wrap(name~parent_strain + pig, ncol = 2)
# 
# data_fit_we = data_fit  %>% ungroup() %>% select(-prev_in_pig) %>% # remove the data points
#   group_by(parent_strain, time, name, prev_in_sample) %>% 
#   summarise(weighted_prob = sum(prob)/4, .groups = "drop") #%>% # equal weight to each piglet: divide by 4
# 
# data_fit_we %>% group_by(parent_strain, time, name) %>% summarise(sum(weighted_prob))
# 
# data_fit_we %>% filter(parent_strain == 1, time == 384, name == "p2")
# data_fit %>% filter(parent_strain == 1, time == 384, name == "p2") %>% print(n=Inf)
# 
# 
# 
# ggplot(data_fit_we %>% filter(name == "phi6"), aes(x=prev_in_pig, y = factor(time), group = time, 
#                                                    height = weighted_prob)) +
#   geom_ridgeline(aes(fill = time)) +
#   facet_wrap(name~parent_strain, ncol = 2)
# ggsave("plots/phi6_fit_over_time.pdf")
# 
# ggplot(data_fit_we %>% filter(name == "p1"), aes(x=prev_in_sample, y = factor(time), group = time, 
#                                                  height = weighted_prob)) +
#   geom_ridgeline(aes(fill = time)) +
#   facet_wrap(name~parent_strain, ncol = 2)
# ggsave("plots/p1_fit_over_time.pdf")
# 
# ggplot(data_fit_we %>% filter(name == "p2"), aes(x=prev_in_pig, y = factor(time), group = time, 
#                                                  height = weighted_prob)) +
#   geom_ridgeline(aes(fill = time)) +
#   facet_wrap(name~parent_strain, ncol = 2)
# ggsave("plots/p2_fit_over_time.pdf")
# 
# ggplot(data_fit_we, aes(x=prev_in_pig, y = factor(time), group = time, 
#                         height = weighted_prob)) +
#   geom_ridgeline(aes(fill = name)) +
#   facet_wrap(name~parent_strain, ncol = 2) + 
#   scale_y_discrete("Time") + 
#   scale_x_continuous("Prevalence in piglet") 
# ggsave("plots/fit_over_time.pdf")
# 
# write.csv(data_fit_we,"data_to_fit.csv")
# 
# 
# #### OLD 
# ### Filter 
# # conversion_data = seen_predicted %>% #filter(prob > 0.05) %>% # Want to only have >1% probability
# #   select(prev, n_colonies, prob) %>% mutate(prev_in_pig = prev/100,
# #                                             prev_in_sample = n_colonies / 20) %>%
# #   select(prev_in_pig, prev_in_sample, prob)
# 
# conversion_data = seen_predicted %>% filter(prob > 0.01)
# 
# element_data = pigg_elements %>% mutate(prev_in_pig = round(sum_prop,2)) %>% select(-sum_prop)
# 
# data_fit = left_join(element_data, conversion_data, by = "prev_in_pig") %>%
#   #  select(-prev_in_sample) %>%
#   group_by(parent_strain, time, name, prev_in_pig)
# 
# # Should add up to 1
# tail(data_fit %>% group_by(pig, prev_in_sample, parent_strain, time, name) %>% summarise(sum(prob)))
# 
# ## data_fit has 
# # prev_in_pig = what the data says
# # prev_in_sample = if see this prevalence = prob it is true
# 
# 
# # Check what adding together for phi6
# theme_set(theme_bw(base_size = 6))
# ggplot(data_fit %>% filter(name == "phi6"), aes(x=prev_in_pig, y = factor(time), group = time, 
#                                                 height = prob)) +
#   geom_ridgeline(aes(fill = time)) +
#   facet_wrap(name~parent_strain + pig, ncol = 2)
# 
# data_fit_we = data_fit  %>% ungroup() %>% select(-prev_in_pig) %>% # remove the data points
#   group_by(parent_strain, time, name, prev_in_sample) %>% 
#   summarise(weighted_prob = sum(prob)/4, .groups = "drop") #%>% # equal weight to each piglet: divide by 4
# 
# data_fit_we %>% group_by(parent_strain, time, name) %>% summarise(sum(weighted_prob))
# 
# data_fit_we %>% filter(parent_strain == 1, time == 384, name == "p2")
# data_fit %>% filter(parent_strain == 1, time == 384, name == "p2") %>% print(n=Inf)
# 
# 
# 
# ggplot(data_fit_we %>% filter(name == "phi6"), aes(x=prev_in_pig, y = factor(time), group = time, 
#                                                    height = weighted_prob)) +
#   geom_ridgeline(aes(fill = time)) +
#   facet_wrap(name~parent_strain, ncol = 2)
# ggsave("plots/phi6_fit_over_time.pdf")
# 
# ggplot(data_fit_we %>% filter(name == "p1"), aes(x=prev_in_sample, y = factor(time), group = time, 
#                                                  height = weighted_prob)) +
#   geom_ridgeline(aes(fill = time)) +
#   facet_wrap(name~parent_strain, ncol = 2)
# ggsave("plots/p1_fit_over_time.pdf")
# 
# ggplot(data_fit_we %>% filter(name == "p2"), aes(x=prev_in_pig, y = factor(time), group = time, 
#                                                  height = weighted_prob)) +
#   geom_ridgeline(aes(fill = time)) +
#   facet_wrap(name~parent_strain, ncol = 2)
# ggsave("plots/p2_fit_over_time.pdf")
# 
# ggplot(data_fit_we, aes(x=prev_in_pig, y = factor(time), group = time, 
#                         height = weighted_prob)) +
#   geom_ridgeline(aes(fill = name)) +
#   facet_wrap(name~parent_strain, ncol = 2) + 
#   scale_y_discrete("Time") + 
#   scale_x_continuous("Prevalence in piglet") 
# ggsave("plots/fit_over_time.pdf")
# 
# write.csv(data_fit_we,"data_to_fit.csv")