##### LHS search 

library(data.table)
library(janitor)
library(foreach)
library(doParallel)
library(tidyverse)
library(lhs)
#numCores <- parallel::detectCores() - 2
#numCores

## Functions
source("code/functions.R")


## Theta that works: 
init.theta = c(mu2 = 0.188878541032485, mu5 = 0.191811127529857, mu6 = 1.06342477430845,
               mu7 = 0.650240103957843, mu8 = 1.13092685382576, mu10 = 1.36327733057528,
               gamma2 = 0.299144953920434, gamma5 = 1.10429707130617, gamma6 = 1.15332183481689,
               gamma7 = 0.586259749844624, gamma8 = 1.42781770425473, gamma10 = 0.669795153200574,
               f2 = 0.0752560231516298, f5 = 0.149273263311754, f6 = 0.104081458666405,
               f7 = 0.00298936941316733, f8 = 0.29765944219356, f10 = 0.30328253933773,
               grow = 1.72361721520587)

mu_mean = mean(init.theta[c("mu2","mu5","mu6","mu7","mu8","mu10")])
gamma_mean = mean(init.theta[c("gamma2","gamma5","gamma6","gamma7","gamma8","gamma10")])
f_mean = mean(init.theta[c("f2","f5","f6","f7","f8","f10")])

#Parameter ranges for LHS;

param_ranges <- as.data.frame(cbind(c(rep(0.0000001,10),
                                      rep(0.000000001,10),
                                      rep(0.0001,10),
                                      0),
                                    c(rep(0.0001,10),
                                      rep(0.001,10),
                                      rep(0.5,10),
                                      0.15)))
colnames(param_ranges) <- c("min","max")

# ##### NO FITNESS COST
# param_ranges <- as.data.frame(cbind(c(rep(0.0000000001,10), # loss
#                                       rep(0.000000001,10), # gain
#                                       rep(0,10), # fitness cost
#                                       0),
#                                     c(rep(0.000001,10),
#                                       rep(0.001,10),
#                                       rep(0,10),
#                                       0.15)))
# colnames(param_ranges) <- c("min","max")
# 
# ##### P4 exploration
# param_ranges <- as.data.frame(cbind(c(c(rep(0.0000000001,9),0.001), # loss
#                                       rep(0.000000001,10), # gain
#                                       c(rep(0,9),0), # fitness cost
#                                       0),
#                                     c(c(rep(0.00001,9),0.3),
#                                       rep(0.001,10),
#                                       c(rep(0,9),1),
#                                       0.15)))
# colnames(param_ranges) <- c("min","max")
# # For p4
# # High loss, equal fitness
# Y[1,10] <- 0.26
# Y[1,30] <- 0
# # Low lost, high fitness cost
# Y[2,10] <- 0.0001
# Y[2,30] <- 0.5


# From the ?lhs example page
# transform a Latin hypercube
nsamples <- 50000
nparameters <- 31
X <- randomLHS(nsamples, nparameters) # first = number of samples, 1000 maybe? second = number of parameters, here 2
Y <- matrix(0, nrow=nsamples, ncol=nparameters)
# Assume parameters uniformly arranged over the range (could assume normal etc if have evidence...)
for(ii in 1:nparameters){
  Y[,ii] <- qunif(X[,ii], min = param_ranges[ii,1], max = param_ranges[ii,2])
}
write.csv(Y, "fits/paraset.csv")

#store_para <- c() # which samples are good?




#### Parallel

#registerDoParallel(numCores)

#foreach (ii=1:nsamples) %dopar% {
for(ii in 1:nsamples){
  
  print(c("Samples number: ", ii))
  
  # gamma_here <- Y[ii,1:10]
  # mu_here <- Y[ii,11:20]
  # growth_here = Y[ii,21:30]
  # grate_here = Y[ii,31]
  
  parameters <- c(mu2 = Y[ii,2],mu5 = Y[ii,5],mu6 = Y[ii,6],
                  mu7 = Y[ii,7],mu8 = Y[ii,8],mu10 = Y[ii,10],
                  gamma2 = Y[ii,12],gamma5 = Y[ii,15],gamma6 = Y[ii,16],
                  gamma7 = Y[ii,17],gamma8 = Y[ii,18],gamma10 = Y[ii,20],
                  f2 = Y[ii,22],f5 = Y[ii,25],f6 = Y[ii,26],
                  f7 = Y[ii,27],f8 = Y[ii,28],f10 = Y[ii,30],
                  grow = Y[ii,31])
  
  ## Time step = 1 hr
  tsteps = 16 * 24
  
  ## Run
  print(Y[ii,])
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  out <- run_sim(tsteps, parameters)
  
  ## Look at output
  P_all <- out$P_all
  Q_all <- out$Q_all
  P_all$parent = "S0385"
  Q_all$parent = "H398"
  P_all$para <- ii
  Q_all$para <- ii
  # Both <- rbind(P_all, Q_all)
  # ggplot(Q_all %>% ungroup() %>% group_by(time) %>% summarise(p4_prev=sum(v10*freq)/sum(freq)),
  #        aes(x=time, y = p4_prev)) + geom_point() + geom_line()
  # ggplot(Both, aes(x=time, y = freq, group = label)) + geom_line(aes(col = label)) + facet_wrap(~parent,ncol = 1)
  error <- c(ii, out$error)
  
  # If a fit then reaches max time so store:
  if(max(P_all$time) == tsteps){
    write.csv(P_all,paste0("fits/P_all_",ii,".csv")) 
    write.csv(Q_all,paste0("fits/Q_all_",ii,".csv")) 
  }else{write.csv(error,paste0("not_fit/",ii,".csv"))}
  
  
}

#stopImplicitCluster()

##### 210806: running at about 10% acceptance! 
name_runs <- "_210806"
### Read in 
errors <-list.files(path="not_fit/", full.names = TRUE) %>%
  lapply(fread) %>% t() %>% 
  bind_rows
errors <- as.data.frame(errors)
errors2 <- as.data.frame(cbind(errors[seq(1,dim(errors)[1],2),"x"],errors[seq(2,dim(errors)[1],2),"x"]))
colnames(errors2) <- c("samp_no","error")
table(errors2$error)
write.csv(errors2, paste("fits",name_runs,"/not_fits_errors.csv"))

# 210806: most errors: fewer than 3 profiles at day 2 /  too fast / P less than 40% / too slow 

### Explore space
Y <- read.csv(paste0("fits",name_runs,"/paraset.csv"))[,-1]
Y <- as.data.frame(Y)
Y <- Y[,c(2,5:8,10,12,15:18,20,25:28,30:31)]
Y$error <- 0
Y[errors2$samp_no,"error"] <- errors2$error
Ylong <- Y %>% pivot_longer(cols = c("V2":"V31"))

Ylong$name <- factor(Ylong$name, levels = c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10",
                                            "V11","V12","V13","V14","V15","V16","V17","V18","V19","V20",
                                            "V21","V22","V23","V24","V25","V26","V27","V28","V29","V30","V31"))

ggplot(Ylong, 
       aes(x = value, y = error)) + 
  geom_point() + 
  facet_wrap(~name, ncol = 10, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(Ylong %>% filter(name %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10",
                                    "V11","V12","V13","V14","V15","V16","V17","V18","V19","V20")), 
       aes(x = value, y = error)) + 
  geom_point() + 
  facet_wrap(~name, ncol = 10, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(Ylong %>% filter(error == 0), aes(x=value, group = name)) + geom_point()
### Correlations between values that reached time 
myd <- Y%>% filter(error == 0) %>% select(-c("error"))
cormat <- round(cor(myd),2) 
melted_cormat <- reshape2::melt(cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + geom_tile()

## Do the LHS fit runs actually come near the data? 
files <- list.files(path=paste0("fits",name_runs), pattern = "_all", full.names = TRUE)
fits <-  files %>%
  lapply(fread) %>% bind_rows(.id = "ID")

fits_long <- fits %>% filter(time %in% c(4,48,96,288,384)) %>% select(SCCmec:freq, time, parent, ID) %>% pivot_longer(cols = SCCmec : p4) %>%
  mutate(n_bugs = freq * value) %>% group_by(time, parent, ID) %>% mutate(total = sum(freq)/10) %>%
  ungroup() %>% mutate(prev = round(n_bugs / total,5))

theme_set(theme_bw(base_size = 8))
ggplot(fits_long, aes(x=time, y = prev, group = interaction(parent, name, ID))) + 
  geom_point(aes(col = ID)) + 
  geom_line(aes(col = ID)) + 
  facet_wrap(name~parent, ncol = 2) + 
  theme(legend.position = "none")
ggsave("plots/lhs_fits_210806.pdf")

### how many above 5%? 

