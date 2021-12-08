##### Example / play code

### Find fits
library(zoo)
library(patchwork)
library(tidyverse)
library(mvtnorm)
theme_set(theme_bw(base_size = 11))
source("code/functions.R")

### Data 
data <- read.csv("data/data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))

### Likelihood 
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

# Just focus on 6 elements that actually move 
#c("phi6","phi2","p1","p2","p3","p4")
theta = c(# gain
  mu2 = 0.008,mu5 = 0.008,
  mu6 = 0.006,mu7 = 0.008,mu8 = 0.000195,mu10 = 0.0008,
  # loss
  gamma2 = 0.00004,gamma5 = 0.000009,
  gamma6 = 0.00039,gamma7 = 0.00019,gamma8 = 0.0003,gamma10 = 0.3,
  # fitness cost
  f2 = 0.6,f5 = 0.06,
  f6 = 0.68,f7 = 0.6,f8 = 0.5,f10 = 0.7,
  grow = 0.09)

theta0 = c(# gain
  mu2 = 0.08,mu5 = 0.08,
  mu6 = 0.00006,mu7 = 0.008,mu8 = 0.0006,mu10 = 0.0008,
  # loss
  gamma2 = 0.00000000004,gamma5 = 0.00000000004,
  gamma6 = 0.008,gamma7 = 0.00000000004,gamma8 = 0.00000000004,gamma10 = 0.1,
  # fitness
  f2 = 0.00002,f5 = 0.006,
  f6 = 0.000068,f7 = 0.006,f8 = 0.005,f10 = 0.0063,
  grow = 0.142)

#c("phi6","phi2","p1","p2","p3","p4")
# worked 8th Sept to 303
# theta = c(# gain
#   mu2 = 0.08,mu5 = 0.008,
#   mu6 = 0.002,mu7 = 0.008,mu8 = 0.00595,mu10 = 0.08,
#   # loss
#   gamma2 = 0.00004,gamma5 = 0.000009,
#   gamma6 = 0.00039,gamma7 = 0.019,gamma8 = 0.0003,gamma10 = 0.003,
#   # fitness
#   f2 = 0.02,f5 = 0.05,
#   f6 = 0.26,f7 = 0.8,f8 = 0.1,f10 = 0.25,
#   grow = 0.142)

## Gets to 384
theta = c(# gain
  mu2 = 0.08,mu5 = 0.008,
  mu6 = 0.002,mu7 = 0.008,mu8 = 0.00595,mu10 = 0.08,
  # loss
  gamma2 = 0.04,gamma5 = 0.0009,
  gamma6 = 0.0039,gamma7 = 0.019,gamma8 = 0.00003,gamma10 = 0.0003,
  # fitness
  f2 = 0.03,f5 = 0.05,
  f6 = 0.26,f7 = 0.9,f8 = 0.1,f10 = 0.8,
  grow = 0.242)

names <- names(theta)
theta_fit = c(0.0886372585286972,
          0.0868323051597757,0.0874064038510944,0.0641586697096416,0.0583446399961726,0.0591852470784067,0.014938995871502,
          0.00489269572535201,0.0205840879838839,0.0283551835318363,0.0106046896970885,0.0339511374181242,0.0124465877907019,
          0.00979985126076986,0.0095708660929162,0.0100261619714865,0.0129826217189951,0.00992632192187117,0.103416758890225)
names(theta_fit) <- names

## All same
#para <- 0.1215405 0.04811499 0.04409806 0.09738901
gain_all <- 0.1215405
gamma_all <- 0.04811499
fitness_all <- 0.04409806


theta = c(# gain
  mu2 = gain_all, mu5 = gain_all,
  mu6 = gain_all,mu7 = gain_all,mu8 = gain_all,mu10 = gain_all,
  # loss
  gamma2 = gamma_all,gamma5 = gamma_all,
  gamma6 = gamma_all,gamma7 = gamma_all,gamma8 = gamma_all,gamma10 = gamma_all,
  # fitness
  f2 = fitness_all,f5 = fitness_all,
  f6 = fitness_all,f7 = fitness_all,f8 = fitness_all,f10 = fitness_all,
  grow = 0.09738901)


theta_fit1 <- c(0.0902391342190606,0.0955462697865904,0.106214814356656,0.0320580334507821,0.053400872560795,
                0.036871987212835,0.0145792645731623,0.0113884777998453,0.0502397535174183,0.011160982802245,
                0.0267229541391116,0.0585937074381367,0.0179654291647543,0.00740612615527331,0.00206276494039848,
                0.0149934380222259,0.0117866755813373,0.0108041440759182,0.0999725982538078)
names(theta_fit1)<-names(theta)

tsteps = 384

### theta same for all 
theta = c(mu = 0.1, gamma = 0.1, f = 0.000000001, grow = 0.2)

## Run for these parameters
out <- run_sim(tsteps, theta)
max(out$P_all$time)
out$error

out$P_all %>% filter(time == 25) %>% arrange(desc(freq))
out$Q_all %>% filter(time == 25) %>% arrange(desc(freq))

if(max(out$P_all$time)==tsteps){ # if get to end 
  #### element prevalence from model 
  model_outputp <- out$prev_predict %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
  model_outputp$prev <- round(model_outputp$prev,2)
  
  ## check likelihood works
  #mock_data <- read_csv("data/mock_perfect_data.csv")[,-1]
  #model_outputp <- mock_data
  
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
compare_dat #### - 718 for mock data
likelihood_lookup_elements
likelihood_lookup_totals


model_outputp <- rename(model_outputp, parent_strain = parent)
pigg_elements <- read.csv("data/pigg_elements.csv")[,-1]

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")), 
             aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) + 
  geom_point(aes(col = name),size = 1.5) + 
  facet_wrap(name~parent_strain, ncol = 2)  


totalsp <- read.csv("data/totals_bug.csv")[,-1]
g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) + 
  geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") + 
  scale_y_log10() + 
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

g1 / g2
#ggsave("plots/first_fit.pdf")

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
  filter(time %in% seq(0,385,20)) %>% group_by(time) %>% 
  mutate(total = sum(freq)/10) %>% group_by(time, name) %>% # /10 for total as 10 elements
  mutate(nbugs_with = value * freq) %>% 
  summarise(nbugs_with_total = sum(nbugs_with), total = min(total), # min here but could be max or anything - just need to move total into summary
            prev = ifelse(nbugs_with_total == 0, 0, nbugs_with_total / total),.groups = "drop") %>%
  select(time,name,prev) %>% mutate(parent = 1)

qnew <-  Q_all %>% pivot_longer(cols = c("SCCmec":"p4")) %>% 
  filter(time %in% seq(0,385,20)) %>% group_by(time) %>% 
  mutate(total = sum(freq)/10) %>% group_by(time, name) %>% # /10 for total as 10 elements
  mutate(nbugs_with = value * freq) %>% 
  summarise(nbugs_with_total = sum(nbugs_with), total = min(total), # min here but could be max or anything - just need to move total into summary
            prev = ifelse(nbugs_with_total == 0, 0, nbugs_with_total / total),.groups = "drop") %>%
  select(time,name,prev) %>% mutate(parent = 2)

prev_predict <- rbind(pnew,qnew)
model_outputp <- prev_predict %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
model_outputp <- rename(model_outputp, parent_strain = parent)
pigg_elements <- read.csv("data/pigg_elements.csv")[,-1]

# Totals 
ptotals <- P_all %>% select(time,freq) %>%
  filter(time %in% seq(0,385,20)) %>% group_by(time) %>% 
  summarise(total = sum(freq),.groups = "drop") %>% mutate(parent = 1)
qtotals <- Q_all %>% select(time,freq) %>%
  filter(time %in% seq(0,385,20)) %>% group_by(time) %>% 
  summarise(total = sum(freq),.groups = "drop") %>% mutate(parent = 2)
model_outputt <- rbind(ptotals, qtotals)

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")), 
             aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) + 
  geom_point(aes(col = name),size = 1.5) + 
  facet_wrap(name~parent_strain, ncol = 2)  


totalsp <- read.csv("data/totals_bug.csv")[,-1]
g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) + 
  scale_y_log10() + 
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

g1 / g2
