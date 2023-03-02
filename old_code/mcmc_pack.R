#### MCMC package
library(mcmc)
set.seed(42) # to get reproducible results

Initial.Values = c(mu = 0.1, gamma = 0.1, grow = 0.1)

out <- metrop(run_sim_logPosterior, Initial.Values, 10)
names(out)

# doesn't do adaptive? 

library(fmcmc) #https://uscbiostats.github.io/fmcmc/
library(coda)

ans <- MCMC(run_sim_logPosterior, initial = Initial.Values, nsteps = 1000, 
            kernel  = kernel_adapt())
plot(ans)
summary(ans)
ans <- MCMC(run_sim_logPosterior, initial = Initial.Values, nsteps = 1000, 
            kernel  = kernel_adapt())
plot(ans[1:4000])
summary(ans)


ans <- MCMC(
  run_sim_logPosterior_para, 
  initial = Initial.Values,
  nsteps  = 1000,
  kernel  = kernel_adapt(),
  nchains = 3,                           # Multiple chains
  conv_checker = convergence_gelman(200), # Checking for conv. every 200 steps
  multicore = TRUE
)

ans <- MCMC(
  run_sim_logPosterior, 
  initial = rbind(Initial.Values,Initial.Values*2, Initial.Values*3),
  nsteps  = 3000,
  kernel  = kernel_adapt(),
  nchains = 3,
  conv_checker = convergence_gelman(200) # Checking for conv. every 200 steps
)

Initial.Values = c(mu2 = 0.265289370464156, mu5 = 0.130932829264206, mu6 = 2.58577755844879, 
                   mu7 = 0.481693463586813, mu8 = 3.65375125367999, mu10 = 1.18123827559817, 
                   gamma2 = 0.592694899220884, gamma5 = 1.11892206012185, gamma6 = 3.75164109872491, 
                   gamma7 = 0.71398191791436, gamma8 = 0.240601848692034, gamma10 = 2.51037647414252, 
                   f2 = 0.0164854283546256, f5 = 0.213784972241758, f6 = 0.0382536099278368, 
                   f7 = 0.0400716377754492, f8 = 0.203513266049325, f10 = 0.466694316467684, 
                   grow = 1.97708717931669, rel_fit = 0.8)
out <- piglet_mrsa_movement(tsteps, Initial.Values/3, ini$bacteria, ini$difference_list)

ans <- MCMC(
  run_sim_logPosterior, 
  initial = rbind(Initial.Values,Initial.Values/2, Initial.Values/3),
  nsteps  = 3000,
  kernel  = kernel_adapt(
    lb = c(rep(0, 12), rep(-1,6), 0.5, 0),
    ub = c(rep(10, 6),rep(1,6), rep(0.1,6), 2.935, 1), # Gain maximum 1 as relative proportion of contacts
  ),
  nchains = 3,
  conv_checker = convergence_gelman(200) # Checking for conv. every 200 steps
)


### From 3000: Gelman-Rubin's R: 7.0103.
list(c(mu2 = 0.205586807324394, mu5 = 0.200943634306506, mu6 = 1.43214227025678, 
       mu7 = 0.509097403589266, mu8 = 1.19672115106187, mu10 = 1.54153649256017, 
       gamma2 = 1.2193937320621, gamma5 = 2.09100294518291, gamma6 = 1.49881135307676, 
       gamma7 = 1.0924492968656, gamma8 = 0.561307180279803, gamma10 = 2.07930061884609, 
       f2 = 0.0853074236054307, f5 = 0.0626762521572307, f6 = 0.190218563660229, 
       f7 = 0.0762184989292291, f8 = 0.0240887795738617, f10 = 0.555270549669397, 
       grow = 2.43013935939726), c(mu2 = 0.20909389060422, mu5 = 0.167631491075283, 
                                   mu6 = 1.30635229807877, mu7 = 0.568472554219395, mu8 = 2.92986609082682, 
                                   mu10 = 2.61678208760249, gamma2 = 0.456217619472711, gamma5 = 1.11631808784579, 
                                   gamma6 = 1.73076185497618, gamma7 = 1.55214373492393, gamma8 = 0.197025469827435, 
                                   gamma10 = 2.52220388856468, f2 = 0.188700636827048, f5 = 0.0971748097659286, 
                                   f6 = 0.0151031924610764, f7 = 0.0588604733332894, f8 = 0.0972174912491175, 
                                   f10 = 0.430369572284168, grow = 2.14864355478075), c(mu2 = 0.221085662945999, 
                                                                                        mu5 = 0.202244255060284, mu6 = 0.695947475302244, mu7 = 0.377348092689187, 
                                                                                        mu8 = 2.2668521458179, mu10 = 2.19856860961439, gamma2 = 0.140640191963348, 
                                                                                        gamma5 = 1.2916744139894, gamma6 = 1.21408085305423, gamma7 = 2.07119013067996, 
                                                                                        gamma8 = 0.216378227186089, gamma10 = 0.482836135660895, f2 = 0.00259591365427901, 
                                                                                        f5 = 0.00470119968868288, f6 = 0.284403735738665, f7 = 0.028780730214476, 
                                                                                        f8 = 0.0514617203691673, f10 = 0.552989217984867, grow = 2.33682122132295
                                   ))


ans_2nd3000 <- MCMC(
  run_sim_logPosterior, 
  initial = ans,
  nsteps  = 3000,
  kernel  = kernel_adapt(
    lb = 0,
    ub = 3.0
  ),
  nchains = 3,
  conv_checker = convergence_gelman(200) # Checking for conv. every 200 steps
)



###### Jan 19th
Initial.Values = c(mu2 = 0.0000265289370464156, mu5 = 0.000000130932829264206, mu6 = 0.01, 
                   mu7 = 0.00000481693463586813, mu8 = 0.01, mu10 = 0.005, 
                   gamma2 = 0.00001, gamma5 = 0.000001, gamma6 = 0.00000000001, 
                   gamma7 = 0.000001, gamma8 = 0.00000001, gamma10 = 0.000000000001, 
                   f2 = -0.3, f5 = -0.3, f6 = 0.9, 
                   f7 = -0.2, f8 = 0.3, f10 = 0.8, 
                   grow = 0.15, rel_fit = 1)

ans <- MCMC(run_sim_logPosterior, 
            initial = Initial.Values, 
            nsteps = 10000, 
            kernel  = kernel_adapt(
              lb = c(rep(1e-12, 12), rep(-1,6), 0.01, 0),
              ub = c(rep(0.01, 6),rep(1,6), rep(1,6), 2.935, 1), # Gain maximum 1 as relative proportion of contacts
            ))
plot(ans)
summary(ans)



ans2 <- MCMC(run_sim_logPosterior, 
            initial = ans, 
            nsteps = 10000, 
            kernel  = kernel_adapt(
              lb = c(rep(0, 12), rep(-1,6), 0.5, 0),
              ub = c(rep(10, 6),rep(1,6), rep(0.1,6), 2.935, 1), # Gain maximum 1 as relative proportion of contacts
            ))
plot(ans2)
summary(ans2)


ans <- MCMC(
  run_sim_logPosterior, 
  initial = rbind(Initial.Values,Initial.Values/2, Initial.Values/3),
  nsteps  = 3000,
  kernel  = kernel_adapt(
    lb = c(rep(0, 12), rep(-1,6), 0.5, 0),
    ub = c(rep(10, 6),rep(1,6), rep(0.1,6), 2.935, 1), # Gain maximum 1 as relative proportion of contacts
  ),
  nchains = 3,
  conv_checker = convergence_gelman(200) # Checking for conv. every 200 steps
)



##### friday 21st
library(adaptMCMC)
samp <- MCMC(run_sim_logPosterior, 
             n=100, 
             init=Initial.Values, 
             scale=c(rep(0.00000000000001,12), rep(0.000000000001,8)),
             adapt=TRUE, 
             acc.rate=0.234)

samp.coda <- convert.to.coda(samp)
class(samp.coda)
## ----------------------
## use functions of package 'coda'
require(coda)
plot(samp.coda)
cumuplot(samp.coda)



samp_lots <- MCMC(run_sim_logPosterior, 
             n=10000, 
             init=Initial.Values, 
             scale=c(rep(0.0000000000000001,12), rep(0.00000000000001,8)),
             adapt=TRUE, 
             acc.rate=0.234)

samp.coda_lots <- convert.to.coda(samp_lots)
write.csv(samp.coda_lots,"fits/220125_10000_samp_all.csv")
plot(samp.coda_lots)
cumuplot(samp.coda_lots)


