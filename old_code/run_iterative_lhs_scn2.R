### Iterative lhs run 
library(data.table)
library(janitor)
library(foreach)
library(doParallel)
library(tidyverse)
library(lhs)
library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)
library(scales)
library(reshape2)
library(dplyr)
library(Rfast)
library(prodlim)
library(ggforce)
library(ggridges)
library(patchwork)
theme_set(theme_bw(base_size = 11))

source("code/iterative_lhs.R")
source("code/piglet_mrsa_functions.R")

### INITIAL CONDITON 1
## Best from full exploration 
m2work <- read_csv("fits/lhs_all/sc2/worked.csv")[,-1]
max(m2work$ll)
ma <- m2work %>% filter(ll > -600, rel_fit < 1)
m <- as.numeric(ma[1,])

Initial.Values = c(mu_phage = m[2], mu_plasmid = m[3], 
                   gamma_phage = m[4],gamma_plasmid = m[5],
                   f_phage = m[6], f_plasmid = m[7],
                   grow = m[8], 
                   rel_fit = m[9])
# check not -Inf
# run_sim_logPosterior(Initial.Values) # -1122

Scen = "scn2"

#run LHS on this + / - 100%
m1 <- lhs_build_run(Initial.Values, limit = 0.5, paste0(Scen,"/iv1"), nsamples = 1e3, ranges = cbind(c(rep(0,4),rep(-1,2), rep(0,2)),c(rep(1,4),rep(1,2),3,1.7)))

# Run LHS on this + / - 100% 
new_iv = m1$max_ll_para
names(new_iv) <- names(Initial.Values)
m2 <- lhs_build_run(new_iv, 0.5, paste0(Scen,"/iv1_max"), nsamples = 1e3, ranges = cbind(c(rep(0,4),rep(-1,2), rep(0,2)),c(rep(1,4),rep(1,2),3,1.7)))

# Run LHS on this + / - 100% 
new_iv = m2$max_ll_para
names(new_iv) <- names(Initial.Values)
m3 <- lhs_build_run(new_iv, 0.5, paste0(Scen,"/iv1_max2"), nsamples = 1e3, ranges = cbind(c(rep(0,4),rep(-1,2), rep(0,2)),c(rep(1,4),rep(1,2),3,1.7)))

## Zoom in on max
all_worked = rbind(m1$worked,m2$worked,m3$worked)

max_para <- as.numeric(all_worked[which.max(all_worked[,1]),2:ncol(all_worked)])
names(max_para) <- names(Initial.Values)
# + / - 50%
mz_1 <- lhs_build_run(max_para, 2, paste0(Scen,"/ivz1"), nsamples = 1e3, ranges = cbind(c(rep(0,4),rep(-1,2), rep(0,2)),c(rep(1,4),rep(1,2),3,1.7)))
# + / - 30% 
new_iv = mz_1$max_ll_para
names(new_iv) <- names(Initial.Values)
mz_2 <- lhs_build_run(new_iv, 3, paste0(Scen,"/ivz2"), nsamples = 1e3, ranges = cbind(c(rep(0,4),rep(-1,2), rep(0,2)),c(rep(1,4),rep(1,2),3,1.7)))
# + / - 10% 
new_iv = mz_2$max_ll_para
names(new_iv) <- names(Initial.Values)
mz_3 <- lhs_build_run(new_iv, 10, paste0(Scen,"/ivz3"), nsamples = 1e3, ranges = cbind(c(rep(0,4),rep(-1,2), rep(0,2)),c(rep(1,4),rep(1,2),3,1.7)))

# Combine all of above
all_worked_with_zoom = as.data.frame(rbind(all_worked, 
                             mz_1$worked,mz_2$worked,mz_3$worked))

#Take top 1000 likelihoods and look at parameter sets
names(all_worked_with_zoom) <- c("ll", names(Initial.Values))

para_gn <- all_worked_with_zoom %>% filter(!ll == -Inf) 
para_gn <- para_gn %>% mutate(paraset = seq(1,nrow(para_gn)))
top_ll <- para_gn[order(para_gn$ll,decreasing = TRUE)[1:min(nrow(para_gn),1000)],] %>%
  pivot_longer(cols=names(Initial.Values))

ggplot(top_ll, aes(x=value,group = name)) + geom_density(alpha = 0.2, aes(fill = name)) + facet_wrap(~name,ncol = 2, scales = "free")
setwd(here::here())
filename = gsub(c(" "), "_", format(as.POSIXct(Sys.time()), tz = "Europe/London", usetz = TRUE))
filename = gsub(":", "-", filename)
ggsave(paste0("fits/",Scen,"lhs_density",filename,".pdf"))

write.csv(all_worked_with_zoom,paste0("fits/",Scen,"all_worked",filename,".csv"))
write.csv(top_ll,paste0("fits/",Scen,"top",filename,".csv"))

ggplot(top_ll, aes(x = ll, y = value, group = name)) + geom_point(aes(col= factor(paraset))) + facet_wrap(~name, scales = "free") + theme(legend.position = "none")

#### Output
top <- para_gn %>% filter(ll > -540) # SORT then look at top 10 and how they are distributed - present this? 


best <- para_gn %>% filter(ll == max(para_gn$ll))
out <- piglet_mrsa_movement(tsteps, best[-c(1,ncol(best))], ini$bacteria, ini$difference_list)


### Likelihood
if(!is.null(out$prev_predict)){
  theta_in <- best[-c(1,ncol(best))]
  #### element prevalence from model (a)
  #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
  model_outputp <- out$prev_predict %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
  model_outputp$prev <- round(model_outputp$value,2)
  
  # Check what distribution of n_colonies at this prevalence in the model pig
  distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, variable, prob_all, n_colonies_prev)
  # e.g. to check 
  #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
  
  # lookup the probability from this distribution for the data
  likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","variable","n_colonies_prev")) %>% 
    mutate(exp_miss = dnorm(prob_all, 1, 0.3)) %>% # experimental measure - 10% around a prob of 1
    summarise(sum(log(exp_miss))) #summarise(sum(log(prob_all)))
  
  #### total bugs output from model (b)
  model_outputt <- out$totl_predict
  
  likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% 
    #mutate(val_in = as.numeric(between(total,min,max))) %>% # instead of 1 / 0 make distance 
    ungroup() %>% 
    mutate(val_in = dnorm(log10(total), mean = (log10(max) - log10(min))/2 + log10(min), sd = (log10(max) - log10(min))/20)) %>% #mean = (max - min)/2 + min, sd = (max - min)/10000)) %>% # instead of 1 / 0 make distance 
    as.data.frame() %>% mutate(likelihood = val_in) %>% summarise(sum(log(likelihood)))
  
  #### ensure certain profiles present (c)
  total_end <- unlist(out$all_results %>% filter(time == tsteps) %>% group_by(parent) %>% summarise(total = sum(value)) %>%select(total))
  ## Need to be present at > 80% and > 20% for parent 1 and parent 2 respectively 
  profile_end <- out$all_results %>% filter(time == tsteps, variable %in% profiles_needed_end) 
  profile_end$variable <- as.numeric(profile_end$variable)
  if(nrow(profile_end) == 3){prof_end <- profile_end$value}else{
    prf_end <- as.data.frame(cbind(profiles_needed_end,c(0,0,0)));colnames(prf_end) <- c("variable","end") 
    p <- left_join(prf_end, profile_end, by = "variable")
    p[which(is.na(p$value)), "value"] <- 0
    prof_end <- p$value
  }
  # If total_end = 0 then prof_end will be 0 too 
  if(total_end[1]>0){ prof_end[1] <-  prof_end[1]/total_end[1]}
  if(total_end[2]>0){ prof_end[2:3] <-  prof_end[2:3]/total_end[2]}
  # Could change sd etc if not close enough 
  likelihood_profile_end <- sum(log(dnorm(prof_end,mean = c(0.8,0.2,0.2), sd = 0.1)))
  
  #### Add in priors
  # Set up for all - if para not there then 0 
  if(length(theta_in) < 32 && length(theta_in) > 10){
    prior.mu = dunif(as.numeric(theta_in["mu2"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["mu5"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["mu6"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["mu7"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["mu8"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["mu10"]), min = 0, max = 1, log = TRUE)
    prior.gamma = dunif(as.numeric(theta_in["gamma2"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["gamma5"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["gamma6"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["gamma7"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["gamma8"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(theta_in["gamma10"]), min = 0, max = 1, log = TRUE)
    prior.f = dnorm(as.numeric(theta_in["f2"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f5"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f6"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f7"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f8"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f10"]), mean = 0, sd = 0.1, log = TRUE) 
    prior.grow = dunif(as.numeric(theta_in[["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(theta_in["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.f + prior.grow + prior.relfit
  }
  
  # If fixed input - same rates for all elements
  if(length(theta_in) == 5){
    prior.mu = dunif(as.numeric(theta_in["mu"]), min = 0, max = 1, log = TRUE) 
    prior.gamma = dunif(as.numeric(theta_in["gamma"]), min = 0, max = 1, log = TRUE) 
    prior.f = dnorm(as.numeric(theta_in["f"]), mean = 0, sd = 0.1, log = TRUE) 
    prior.grow = dunif(as.numeric(theta_in[["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(theta_in["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.f + prior.grow + prior.relfit
  }
  
  # If fixed input - same rates for all elements and no fitness cost 
  if(length(theta_in) == 4){
    prior.mu = dunif(as.numeric(theta_in["mu"]), min = 0, max = 1, log = TRUE) 
    prior.gamma = dunif(as.numeric(theta_in["gamma"]), min = 0, max = 1, log = TRUE) 
    prior.grow = dunif(as.numeric(theta_in[["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(theta_in["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.grow + prior.relfit
  }
  
  
  # If fixed input - same rates for phage vs plasmids
  if(length(theta_in) == 8){
    prior.mu = dunif(as.numeric(theta_in["mu_phage"]), min = 0, max = 1, log = TRUE)  + dunif(as.numeric(theta_in["mu_plasmid"]), min = 0, max = 1, log = TRUE) 
    prior.gamma = dunif(as.numeric(theta_in["gamma_phage"]), min = 0, max = 1, log = TRUE) + dunif(as.numeric(theta_in["gamma_plasmid"]), min = 0, max = 1, log = TRUE) 
    prior.f = dnorm(as.numeric(theta_in["f_phage"]), mean = 0, sd = 0.1, log = TRUE) + dnorm(as.numeric(theta_in["f_plasmid"]), mean = 0, sd = 0.1, log = TRUE)
    prior.grow = dunif(as.numeric(theta_in[["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(theta_in["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.f + prior.grow + prior.relfit
  }
  
  
  # If fixed input - same loss/gain rates different fitness
  if(length(theta_in) == 10){
    prior.mu = dunif(as.numeric(theta_in["mu"]), min = 0, max = 1, log = TRUE) 
    prior.gamma = dunif(as.numeric(theta_in["gamma"]), min = 0, max = 1, log = TRUE) 
    prior.f = dnorm(as.numeric(theta_in["f2"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f5"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f6"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f7"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f8"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(theta_in["f10"]), mean = 0, sd = 0.1, log = TRUE)
    prior.grow = dunif(as.numeric(theta_in[["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(theta_in["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.f + prior.grow + prior.relfit
  }
  
  #### Compare to data 
  compare_dat <- log.prior + likelihood_lookup_elements + likelihood_lookup_totals + 10 * likelihood_profile_end # add in a 10* weight for profile_end as otherwise only a small contribution relatively
}else{compare_dat <- -Inf}

# return log likelihood
compare_dat #### - 718 for mock data
likelihood_lookup_elements
likelihood_lookup_totals
likelihood_profile_end


model_outputp <- rename(model_outputp, parent_strain = parent)
model_outputp$name <- recode(model_outputp$variable, V2 = "phi6",V5 = "phi2",V6 = "p1",V7 = "p2",V8 = "p3",V10 ="p4")

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")),
             aes(x=time, y = sum_prop, group = interaction(name, pig))) +
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) +
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) +
  geom_point(aes(col = name),size = 1.5) +
  facet_wrap(name~parent_strain, ncol = 2) + 
  ggtitle(paste0("likelihood = ", round(compare_dat,3)))

g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) +
  geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") +
  scale_y_log10() +
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

g3 <- ggplot(out$all_results %>% filter(variable %in% profiles_needed_end), aes(x=time, y = value)) + geom_line(aes(col = variable)) + geom_vline(xintercept = tsteps) + 
  geom_hline(yintercept = c(0.05) * max(out$all_results$value)) + scale_color_manual(values = c("red","green","blue"), breaks = c(profiles_needed_end))

g1 / (g2 + g3)
ggsave(paste0("fits/",Scen,"best.pdf"))
