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
