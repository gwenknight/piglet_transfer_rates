
# script for mcmc functions

define_priors = function(parameters_in){
  
  # Set up for all - if para not there then 0 
  if(length(parameters_in) < 32 && length(parameters_in) > 10){
    density = function(parameters_in){
      prior.mu = 
        dlnorm(as.numeric(parameters_in["mu2"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["mu5"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["mu6"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["mu7"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["mu8"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["mu10"]), meanlog = log(1e-3), sdlog = 2, log = TRUE)
      prior.gamma = dlnorm(as.numeric(parameters_in["gamma2"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["gamma5"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["gamma6"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["gamma7"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["gamma8"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["gamma10"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) 
      prior.f = dnorm(as.numeric(parameters_in["f2"]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in["f5"]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in["f6"]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in["f7"]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in["f8"]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in["f10"]), mean = -0.1, sd = 0.1, log = TRUE) 
      prior.grow = dunif(as.numeric(parameters_in["grow"]), min = 0, max = 3, log = TRUE)
      prior.relfit = dnorm(as.numeric(parameters_in["rel_fit"]), mean = 0.8, sd = 0.1, log = TRUE) 
      return(prior.mu+prior.gamma+prior.f+prior.grow+prior.relfit)
    }
    
    sampler = function(n=1){
      cc <- cbind(
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2),
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rnorm(n, mean = -0.1, sd = 0.1), 
        rnorm(n, mean = -0.1, sd = 0.1), 
        rnorm(n, mean = -0.1, sd = 0.1), 
        rnorm(n, mean = -0.1, sd = 0.1), 
        rnorm(n, mean = -0.1, sd = 0.1), 
        rnorm(n, mean = -0.1, sd = 0.1), 
        runif(n, min = 0, max = 3),
        rnorm(n, mean = 1, sd = 0.1) 
      )
      colnames(cc) <- names(parameters_in)
      return(cc)
    }
    
    lower = c(rep(0, 12), rep(-1, 6), 0, 0)
    upper = c(rep(1, 12), rep(1, 6), 3, 2)
    
  } else if(length(parameters_in) == 5){
    # If fixed input - same rates for all elements
    density = function(parameters_in){
      prior.mu = dlnorm(as.numeric(parameters_in["mu"]), meanlog = log(1e-3), sdlog = 2, log = TRUE)
      prior.gamma = dlnorm(as.numeric(parameters_in["gamma"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) 
      prior.f = dnorm(as.numeric(parameters_in["f"]), mean = -0.1, sd = 0.1, log = TRUE) 
      prior.grow = dunif(as.numeric(parameters_in["grow"]), min = 0, max = 3, log = TRUE)
      prior.relfit = dnorm(as.numeric(parameters_in["rel_fit"]), mean = 0.8, sd = 0.1, log = TRUE) 
      return(prior.mu+prior.gamma+prior.f+prior.grow+prior.relfit)
    }
    
    sampler = function(n=1){
      cc <- cbind(
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rnorm(n, mean = -0.1, sd = 0.1), 
        runif(n, min = 0, max = 3),
        rnorm(n, mean = 1, sd = 0.1) 
      )
      colnames(cc) <- names(parameters_in)
      return(cc)
    }
    
    lower = c(0, 0, -1, 0, 0)
    upper = c(1, 1, 1, 3, 2)
    
  } else if(length(parameters_in) == 4){
    # If fixed input - same rates for all elements and no fitness cost 
    density = function(parameters_in){
      prior.mu = dlnorm(as.numeric(parameters_in["mu"]), meanlog = log(1e-3), sdlog = 2, log = TRUE)
      prior.gamma = dlnorm(as.numeric(parameters_in["gamma"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) 
      prior.grow = dunif(as.numeric(parameters_in["grow"]), min = 0, max = 3, log = TRUE)
      prior.relfit = dnorm(as.numeric(parameters_in["rel_fit"]), mean = 0.8, sd = 0.1, log = TRUE) 
      return(prior.mu+prior.gamma+prior.grow+prior.relfit)
    }
    
    sampler = function(n=1){
      cc <- cbind(
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        runif(n, min = 0, max = 3),
        rnorm(n, mean = 1, sd = 0.1) 
      )
      colnames(cc) <- names(parameters_in)
      return(cc)
    }
    
    lower = c(0, 0, 0, 0)
    upper = c(1, 1, 3, 2)
    
  } else if(length(parameters_in) == 8){
    # If fixed input - same rates for phage vs plasmids
    density = function(parameters_in){
      prior.mu = dlnorm(as.numeric(parameters_in["mu_phage"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) + 
        dlnorm(as.numeric(parameters_in["mu_plasmid"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) 
      prior.gamma = dlnorm(as.numeric(parameters_in["gamma_phage"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) +
        dlnorm(as.numeric(parameters_in["gamma_plasmid"]), meanlog = log(1e-3), sdlog = 2, log = TRUE) 
      prior.f = dnorm(as.numeric(parameters_in["f_phage"]), mean = -0.1, sd = 0.1, log = TRUE) +
        dnorm(as.numeric(parameters_in["f_plasmid"]), mean = -0.1, sd = 0.1, log = TRUE)
      prior.grow = dunif(as.numeric(parameters_in["grow"]), min = 0, max = 3, log = TRUE)
      prior.relfit = dnorm(as.numeric(parameters_in["rel_fit"]), mean = 0.8, sd = 0.1, log = TRUE)
      return(prior.mu+prior.gamma+prior.f+prior.grow+prior.relfit)
    }
    
    sampler = function(n=1){
      cc <- cbind(
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rnorm(n, mean = -0.1, sd = 0.1), 
        rnorm(n, mean = -0.1, sd = 0.1),
        runif(n, min = 0, max = 3),
        rnorm(n, mean = 1, sd = 0.1)
      )
      colnames(cc) <- names(parameters_in)
      return(cc)
    }
    
    lower = c(0, 0, 0, 0, -1, -1, 0, 0)
    upper = c(1, 1, 1, 1, 1, 1, 3, 2)
    
  } else if(length(parameters_in) == 10){
    # If fixed input - same loss/gain rates different fitness
    density = function(parameters_in){
      prior.mu = dlnorm(as.numeric(parameters_in[1]), meanlog = log(1e-3), sdlog = 2, log = TRUE)
      prior.gamma = dlnorm(as.numeric(parameters_in[2]), meanlog = log(1e-3), sdlog = 2, log = TRUE) 
      prior.f = dnorm(as.numeric(parameters_in[3]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in[4]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in[5]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in[6]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in[7]), mean = -0.1, sd = 0.1, log = TRUE) + 
        dnorm(as.numeric(parameters_in[8]), mean = -0.1, sd = 0.1, log = TRUE)
      prior.grow = dunif(as.numeric(parameters_in[9]), min = 0, max = 3, log = TRUE)
      prior.relfit = dnorm(as.numeric(parameters_in[10]), mean = 0.8, sd = 0.1, log = TRUE)
      return(prior.mu+prior.gamma+prior.f+prior.grow+prior.relfit)
    }
    
    sampler = function(n=1){
      cc <- cbind(
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rlnorm(n, meanlog = log(1e-3), sdlog = 2), 
        rnorm(n, mean = -0.1, sd = 0.1),
        rnorm(n, mean = -0.1, sd = 0.1),
        rnorm(n, mean = -0.1, sd = 0.1),
        rnorm(n, mean = -0.1, sd = 0.1),
        rnorm(n, mean = -0.1, sd = 0.1),
        rnorm(n, mean = -0.1, sd = 0.1),
        runif(n, min = 0, max = 3),
        rnorm(n, mean = 1, sd = 0.1)
      )
      colnames(cc) <- names(parameters_in)
      return(cc)
    }
    
    lower = c(0, 0, -1, -1, -1, -1, -1, -1, 0, 0)
    upper = c(1, 1, 1, 1, 1, 1, 1, 1, 3, 2)
    
  } else stop("Invalid number of parameters!")
  
  priors = createPrior_gk(density = density,
                          sampler = sampler,
                          lower = lower, 
                          upper = upper,
                          best = parameters_in)
  
  return(priors)
  
}



run_sim_logPosterior = function(parameters_in){
  
  tsteps = 384
  
  ## Run for these parameters
  out = piglet_mrsa_movement(tsteps, parameters_in, bacteria, difference_list)
  
  ### Likelihood
  if(!is.null(out$prev_predict)){
    
    #### element prevalence from model (a)
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    model_outputp = out$prev_predict %>%
      filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
    model_outputp$prev = round(model_outputp$value,2)
    
    # Check what distribution of n_colonies at this prevalence in the model pig
    distributs = left_join(model_outputp, dist_like, by = "prev") %>%
      select(parent, time, variable, prob_all, n_colonies_prev)
    # e.g. to check 
    #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
    
    # lookup the probability from this distribution for the data
    likelihood_lookup_elements = left_join(data_6, distributs,
                                           by = c("parent","time","variable","n_colonies_prev")) %>% 
      mutate(exp_miss = dnorm(prob_all, 1, 0.3)) %>% # experimental measure - 10% around a prob of 1
      pull(exp_miss) %>%
      log %>%
      sum
    
    #### total bugs output from model (b)
    model_outputt = out$totl_predict
    
    likelihood_lookup_totals = left_join(model_outputt, totals, by = c("time", "parent")) %>%
      mutate(val_in = dnorm(log10(total),
                            mean = (log10(max) - log10(min))/2 + log10(min),
                            sd = (log10(max) - log10(min))/20)) %>% #mean = (max - min)/2 + min, sd = (max - min)/10000)) %>% # instead of 1 / 0 make distance 
      pull(val_in) %>%
      log %>%
      sum
    
    #### ensure certain profiles present (c)
    total_end = out$all_results %>%
      filter(time == tsteps) %>%
      group_by(parent) %>%
      summarise(total = sum(value)) %>%
      pull(total)
    
    # Sometimes one parent lost: then only one value in total_end
    if(length(total_end)==1){parent_at_end = as.numeric(out$all_results %>% filter(time == tsteps) %>% select(parent) %>% unique())}
    
    ## Need to be present at > 80% and > 20% for parent 1 and parent 2 respectively 
    profile_end = out$all_results %>%
      filter(time == tsteps,
             variable %in% profiles_needed_end) %>%
      mutate(variable = as.numeric(variable))
    
    if(nrow(profile_end) == 3){
      prof_end = profile_end$value
    } else {
      prf_end = as.data.frame(cbind(profiles_needed_end,
                                    c(0,0,0)))
      colnames(prf_end) = c("variable","end") 
      p = left_join(prf_end, profile_end, by = "variable")
      p[which(is.na(p$value)), "value"] <- 0
      prof_end = p$value
    }
    
    # If total_end = 0 then prof_end will be 0 too 
    if(length(total_end) > 1){
      if(total_end[1]>0) prof_end[1] = prof_end[1]/total_end[1]
      if(total_end[2]>0) prof_end[2:3] = prof_end[2:3]/total_end[2]
    } else {
      if(parent_at_end == 1){prof_end[1] = prof_end[1]/total_end[1]; prof_end[2:3] = 0}
      if(parent_at_end == 2){prof_end[1] = 0; prof_end[2:3] = prof_end[2:3]/total_end[2]}
    }
    # Could change sd etc if not close enough 
    #likelihood_profile_end <- sum(log(dnorm(prof_end,mean = c(0.8,0.2,0.2), sd = 0.5)))
    ## Strain 1 has to be mostly the parent strain! 
    ## Strain 2: two profiles are highly dominant: sum of these has to be > 0.5
    likelihood_profile_end = log(dnorm(prof_end[1], mean = 0.975, sd = 0.0288)) +
      log(dtruncnorm(sum(prof_end[2:3]), a = 0, b = 1, mean = 0.75, sd = 0.17))
    
    #### Compare to data 
    compare_dat = likelihood_lookup_elements + likelihood_lookup_totals + likelihood_profile_end # add in a 10* weight for profile_end as otherwise only a small contribution relatively
    
  } else {
    compare_dat = -Inf
  }
  
  # return log likelihood
  compare_dat
}
