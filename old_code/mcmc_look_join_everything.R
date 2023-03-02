#### mcmc_everything
# Join chains together

f1 <- "2021-12-16_23-38-55_GMT"
f2 <- "2021-12-16_19-47-33_GMT"
f3 <- "2021-12-16_21-11-08_GMT"
#f4 <- "2021-12-16_08-14-05_GMT"
f5 <- "2021-12-16_03-12-04_GMT"


mcmc.epi_every <- c()
t1 <- read.csv(here::here("fits/everything",paste0(f1,"_","trace",".csv")))[,-1]
t2 <- read.csv(here::here("fits/everything",paste0(f2,"_","trace",".csv")))[,-1]
t3 <- read.csv(here::here("fits/everything",paste0(f3,"_","trace",".csv")))[1:5000,-1] # gets stuck
#t4 <- read.csv(here::here("fits/everything",paste0(f4,"_","trace",".csv")))[,-1]
t5 <- read.csv(here::here("fits/everything",paste0(f5,"_","trace",".csv")))[,-1]

t12 <- rbind(t1,t2)
t123 <- rbind(t12, t3)
#t1234 <- rbind(t123, t4)
mcmc.epi_every$trace <- rbind(t123, t5)


mcmc.epi_every$acceptance.rate <- rbind(c(read.csv(here::here("fits/everything",paste0(f1,"_","acceptance_rates",".csv"))),
                                          read.csv(here::here("fits/everything",paste0(f2,"_","acceptance_rates",".csv"))),
                                          read.csv(here::here("fits/everything",paste0(f3,"_","acceptance_rates",".csv")))))
mcmc.epi_every$covmat.empirical <- rbind(c(read.csv(here::here("fits/everything",paste0(f1,"_","covmat",".csv"))),
                                           read.csv(here::here("fits/everything",paste0(f2,"_","covmat",".csv"))),
                                           read.csv(here::here("fits/everything",paste0(f3,"_","covmat",".csv")))))




# # # Look at output
mcmc.trace <- mcmc(mcmc.epi_every$trace)
summary(mcmc.trace)
acceptanceRate <- 1 - rejectionRate(mcmc.trace)
acceptanceRate
plot(mcmc.trace)
effectiveSize(mcmc.trace)

mcmc.epi_every$covmat.empirical

mcmc.trace.burned <- burnAndThin(mcmc.trace, burn = 2000)
plot(mcmc.trace.burned)

autocorr.plot(mcmc.trace.burned)
#
mcmc.trace.burned.thinned <- burnAndThin(mcmc.trace.burned, thin = 5)
autocorr.plot(mcmc.trace.burned.thinned)

plotESSBurn(mcmc.trace)

tail(mcmc.trace)
# #### Fit 15th Dec
parameters_every <- mcmc.trace[nrow(mcmc.trace),1:(ncol(mcmc.trace)-1)]

## Run for these parameters
out <- piglet_mrsa_movement(tsteps, parameters_every)

### Data
data <- read.csv("data/data_to_fit.csv")[,-1]
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
### Likelihood
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>%
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))


if(!is.null(out$prev_predict)){
  
  #### element prevalence from model 
  #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
  model_outputp <- out$prev_predict %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
  model_outputp$prev <- round(model_outputp$value,2)
  
  # Check what distribution of n_colonies at this prevalence in the model pig
  distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, variable, prob_all, n_colonies_prev)
  # e.g. to check 
  #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
  
  # lookup the probability from this distribution for the data
  likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","variable","n_colonies_prev")) %>% summarise(sum(log(prob_all)))
  
  #### total bugs output from model 
  model_outputt <- out$totl_predict
  
  likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
    as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))
  
  #### Compare to data 
  compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals
}else{compare_dat <- -Inf}
# return log likelihood
as.numeric(compare_dat)

# return log likelihood
compare_dat #### - 718 for mock data
likelihood_lookup_elements
likelihood_lookup_totals


model_outputp <- rename(model_outputp, parent_strain = parent)
model_outputp$name <- recode(model_outputp$variable, V2 = "phi6",V5 = "phi2",V6 = "p1",V7 = "p2",V8 = "p3",V10 ="p4")

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
ggsave(here::here("fits/everything/",paste0(filename,"_","fit",".pdf")))

dput(mcmc.trace[nrow(mcmc.trace),1:(ncol(mcmc.trace)-1)])