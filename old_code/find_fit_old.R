##### Example / play code

### Find fits
library(zoo)
library(patchwork)

### Data 
data <- read.csv("data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))

### Likelihood 
dist_like <- read.csv("seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

# Just focus on 6 elements that actually move 
#c("phi6","phi2","p1","p2","p3","p4")
theta = c(# gain
  mu2 = 0.008,mu5 = 0.008,
  mu6 = 0.006,mu7 = 0.008,mu8 = 0.000195,mu10 = 0.0008,
  # loss
  gamma2 = 0.00004,gamma5 = 0.000009,
  gamma6 = 0.00039,gamma7 = 0.00019,gamma8 = 0.0003,gamma10 = 0.23,
  # fitness cost
  f2 = 0.6,f5 = 0.06,
  f6 = 0.68,f7 = 0.6,f8 = 0.5,f10 = 0.5,
  grow = 0.09)

theta = c(# gain
  mu2 = 0.008,mu5 = 0.008,
  mu6 = 0.006,mu7 = 0.008,mu8 = 0.000195,mu10 = 0.0008,
  # loss
  gamma2 = 0.00000004,gamma5 = 0.00000009,
  gamma6 = 0.0000039,gamma7 = 0.0000019,gamma8 = 0.000003,gamma10 = 0.000000023,
  # fitness
  f2 = 0.002,f5 = 0.006,
  f6 = 0.000068,f7 = 0.25,f8 = 0.005,f10 = 0.0063,
  grow = 0.142)


tsteps = 384

## Run for these parameters
out <- run_sim(tsteps, theta)
max(out$P_all$time)
out$error

if(max(out$P_all$time)==tsteps){ # if get to end 
  #### element prevalence from model 
  model_outputp <- out$prev_predict %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
  model_outputp$prev <- round(model_outputp$prev,2)
  
  # Check what distribution of n_colonies at this prevalence in the model pig
  distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, name, prob_all, n_colonies_prev)
  # e.g. to check 
  #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
  
  # lookup the probability from this distribution for the data
  likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","name","n_colonies_prev")) %>% summarise(sum(log(prob_all)))
  
  #### total bugs output from model 
  model_outputt <- out$totl_predict
  
  likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
    as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))
  
  #### Compare to data 
  compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals
  
  compare_dat
  # compare_dat <- as.numeric(left_join(model_outputp, data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup() %>% 
  #                             summarise(sum(log(weighted_prob_all))))  + # log likelihood for elements
  #   as.numeric(left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
  #                as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))) # log likelhood for totals
}else{compare_dat <- as.numeric(-Inf)}
# return log likelihood
compare_dat
likelihood_lookup_elements
likelihood_lookup_totals

# left_join(model_outputp, data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup() %>% print(n=Inf)
# 
# left_join(model_outputp %>% filter(time == 384), data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup() %>% summarise(sum(log(weighted_prob)))
# 
# left_join(model_outputp %>% filter(time == 384), data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup()
# 
# 
# left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max)))
# 
# 

model_outputp <- rename(model_outputp, parent_strain = parent)
pigg_elements <- read.csv("pigg_elements.csv")[,-1]

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")), 
             aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) + 
  geom_point(aes(col = name),size = 1.5) + 
  facet_wrap(name~parent_strain, ncol = 2)  


totalsp <- read.csv("totals_bug.csv")[,-1]
g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) + 
  scale_y_log10() + 
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

g1 / g2
ggsave("plots/first_fit.pdf")

oo <- out$P_all
oo$parent <- 1
pp <- out$Q_all
pp$parent <- 2
both <- rbind(oo, pp)
ggplot(both, aes(x=time, y = freq, group = label)) + geom_line(aes(col = label)) + scale_y_log10() + 
  facet_wrap(~parent)
mt <- max(oo$time) 
100 * sum(oo %>% filter(time == max(oo$time)) %>% select(freq))/(sum(oo %>% filter(time == max(oo$time)) %>% select(freq)) + sum(pp %>% filter(time == max(oo$time)) %>% select(freq)))

total_p <- both %>% group_by(time, parent) %>% summarise(total = sum(freq)) %>% ungroup() %>%
  group_by(time) %>% mutate(freq_q = ifelse(parent == 2, total, 0), sum = sum(total), freqq = max(freq_q), prop_q = 100 * freqq/sum)

ggplot(total_p, aes(x=time, y = prop_q)) + geom_line() + scale_y_continuous(lim = c(0,100)) + geom_hline(yintercept = c(60,40), lty = "dashed")

#### prev predict when time < 384
P_all <- out$P_all; Q_all <- out$Q_all
colnames(P_all) <- c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4","freq","fitness","prob_survive","time","label")
colnames(Q_all) <- c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4","freq","fitness","prob_survive","time","label")
pnew <- P_all %>% pivot_longer(cols = c("SCCmec":"p4")) %>% 
  filter(time %in% c(4,48,72,288,384)) %>% group_by(time) %>% 
  mutate(total = sum(freq)/10) %>% group_by(time, name) %>% # /10 for total as 10 elements
  mutate(nbugs_with = value * freq) %>% 
  summarise(nbugs_with_total = sum(nbugs_with), total = min(total), # min here but could be max or anything - just need to move total into summary
            prev = ifelse(nbugs_with_total == 0, 0, nbugs_with_total / total),.groups = "drop") %>%
  select(time,name,prev) %>% mutate(parent = 1)

qnew <-  Q_all %>% pivot_longer(cols = c("SCCmec":"p4")) %>% 
  filter(time %in% c(4,48,72,288,384)) %>% group_by(time) %>% 
  mutate(total = sum(freq)/10) %>% group_by(time, name) %>% # /10 for total as 10 elements
  mutate(nbugs_with = value * freq) %>% 
  summarise(nbugs_with_total = sum(nbugs_with), total = min(total), # min here but could be max or anything - just need to move total into summary
            prev = ifelse(nbugs_with_total == 0, 0, nbugs_with_total / total),.groups = "drop") %>%
  select(time,name,prev) %>% mutate(parent = 2)

prev_predict <- rbind(pnew,qnew)
model_outputp <- prev_predict %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
model_outputp <- rename(model_outputp, parent_strain = parent)
pigg_elements <- read.csv("pigg_elements.csv")[,-1]

# Totals 
ptotals <- P_all %>% select(time,freq) %>%
  filter(time %in% c(4,48,72,288,384)) %>% group_by(time) %>% 
  summarise(total = sum(freq),.groups = "drop") %>% mutate(parent = 1)
qtotals <- Q_all %>% select(time,freq) %>%
  filter(time %in% c(4,48,72,288,384)) %>% group_by(time) %>% 
  summarise(total = sum(freq),.groups = "drop") %>% mutate(parent = 2)
model_outputt <- rbind(ptotals, qtotals)

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")), 
             aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) + 
  geom_point(aes(col = name),size = 1.5) + 
  facet_wrap(name~parent_strain, ncol = 2)  


totalsp <- read.csv("totals_bug.csv")[,-1]
g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) + 
  scale_y_log10() + 
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

g1 / g2


# OLD 
# 
# 
# # This gives the distribution for the data = likelihood
# 
# data <- rename(data, parent = parent_strain)
# #data <- rename(data, prev = prev_in_pig)
# data <- rename(data, prev = mean)
# data_6_filter <- data %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
# 
# # for elements over time plot
# pigg_elements <- read.csv("pigg_elements.csv")[,-1]
# 
# ### Interpolate 
# # data_6 <- data_6_filter %>%
# #   group_by(parent, time, name) %>%
# #   complete(prev = seq(0, 1, len = 101)) %>%
# #   mutate(weighted_prob_all = na.approx(weighted_prob)) %>%
# #   select(-weighted_prob)
# # data_6$prev <- round(data_6$prev,2)
# 
# # Also the total number of bugs needs to fit 
# totals <- read.csv("totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
#   pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))
# 
# 
# #c("phi6","phi2","p1","p2","p3","p4")
# theta = c(# gain
#   mu2 = 0.008,mu5 = 0.008,
#   mu6 = 0.0079,mu7 = 0.043,mu8 = 0.000095,mu10 = 0.0008,
#   # loss
#   gamma2 = 0.004,gamma5 = 0.0009,
#   gamma6 = 0.049,gamma7 = 0.089,gamma8 = 0.65,gamma10 = 0.4,
#   # fitness
#   f2 = 0.2,f5 = 0.06,
#   f6 = 0.48,f7 = 0.06,f8 = 0.06,f10 = 0.38,
#   grow = 0.102)
# 
# tsteps = 384
# 
# 
# 
# ## Run for these parameters
# out <- run_sim(tsteps, theta)
# max(out$P_all$time)
# 
# if(max(out$P_all$time)==tsteps){
#   model_outputp <- out$prev_predict %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
#   model_outputp$prev <- round(model_outputp$prev,2)
#   model_outputt <- out$totl_predict
#   
#   ## Compare to data 
#   compare_dat <- as.numeric(left_join(model_outputp, data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup() %>% 
#                               summarise(sum(log(weighted_prob_all))))  + # log likelihood for elements
#     as.numeric(left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
#                  as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))) # log likelhood for totals
# }else{compare_dat <- as.numeric(-Inf)}
# # return log likelihood
# compare_dat
# 
# # left_join(model_outputp, data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup() %>% print(n=Inf)
# # 
# # left_join(model_outputp %>% filter(time == 384), data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup() %>% summarise(sum(log(weighted_prob)))
# # 
# # left_join(model_outputp %>% filter(time == 384), data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup()
# # 
# # 
# # left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max)))
# # 
# # 
# 
# model_outputp <- rename(model_outputp, parent_strain = parent)
# 
# g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")), 
#              aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
#   geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
#   geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) + 
#   geom_point(aes(col = name),size = 1.5) + 
#   facet_wrap(name~parent_strain, ncol = 2)  
# 
# 
# totals <- read.csv("totals_bug.csv")[,-1]
# g2 <- ggplot(totals, aes(x=time, y = value, group = interaction(time,name,parent))) +
#   geom_point(aes(col = factor(parent))) + 
#   scale_y_log10() + 
#   geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))
# 
# g1 / g2
# ggsave("plots/first_fit.pdf")