
# script for model

piglet_mrsa_movement = function(tsteps, parameters_in, bacteria, difference_list){
  # tsteps = number of time steps
  # parameters_in = parameters needed to run model
  # bacteria = initial distribution and all possible profiles
  # difference_list = possible transitions
  
  #change to assign values to missing parameters based on how many not included
  #start from most complete case
  #eg if length is 5 then need to assign the fixed values to a whole bunch of non-defined params
  
  # Assign parameter values
  if(length(parameters_in) == 32){ # Not included atm
    
    # If have rates for all 10 elements
    rate_loss = as.numeric(parameters_in[1:10])
    rate_gain = as.numeric(parameters_in[11:20])
    fitness_costs = as.numeric(parameters_in[21:30])
    growth_rate = as.numeric(parameters_in[31])
    rel_fit_human = as.numeric(parameters_in[32])
    
  } else if(length(parameters_in) < 32 && length(parameters_in) > 10){
    
    # If just looking at the elements that move
    rate_loss = as.numeric(c(0,parameters_in["mu2"],0,0,parameters_in["mu5"],parameters_in["mu6"],parameters_in["mu7"],parameters_in["mu8"],0,parameters_in["mu10"]))
    rate_gain = as.numeric(c(0,parameters_in["gamma2"],0,0,parameters_in["gamma5"],parameters_in["gamma6"],parameters_in["gamma7"],parameters_in["gamma8"],0,parameters_in["gamma10"]))
    fitness_costs = as.numeric(c(0,parameters_in["f2"],0,0,parameters_in["f5"],parameters_in["f6"],parameters_in["f7"],parameters_in["f8"],0,parameters_in["f10"]))
    growth_rate = as.numeric(parameters_in["grow"])
    rel_fit_human = as.numeric(parameters_in["rel_fit"])
    
  } else if(length(parameters_in) == 5){
    
    # If fixed input - same rates for all elements
    if(is.null(names(parameters_in))){stop("No names for parameters_in")}
    
    rate_loss = as.numeric(c(0,parameters_in["mu"],0,0,parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],0,parameters_in["mu"]))
    rate_gain = as.numeric(c(0,parameters_in["gamma"],0,0,parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],0,parameters_in["gamma"]))
    fitness_costs = as.numeric(c(0,parameters_in["f"],0,0,parameters_in["f"],parameters_in["f"],parameters_in["f"],parameters_in["f"],0,parameters_in["f"]))
    growth_rate = as.numeric(parameters_in["grow"])
    rel_fit_human = as.numeric(parameters_in["rel_fit"])
    
  } else if(length(parameters_in) == 4){
    
    # If fixed input - same rates for all elements and no fitness cost 
    if(is.null(names(parameters_in))){stop("No names for parameters_in")}
    
    rate_loss = as.numeric(c(0,parameters_in["mu"],0,0,parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],0,parameters_in["mu"]))
    rate_gain = as.numeric(c(0,parameters_in["gamma"],0,0,parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],0,parameters_in["gamma"]))
    fitness_costs = as.numeric(c(0,0,0,0,0,0,0,0,0,0))
    growth_rate = as.numeric(parameters_in["grow"])
    rel_fit_human = 1
    
  } else if(length(parameters_in) == 8){
    
    # If fixed input - same rates for phage vs plasmids
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    
    rate_loss = as.numeric(c(0,parameters_in["mu_phage"],0,0,parameters_in["mu_phage"],parameters_in["mu_plasmid"],parameters_in["mu_plasmid"],parameters_in["mu_plasmid"],0,parameters_in["mu_plasmid"]))
    rate_gain = as.numeric(c(0,parameters_in["gamma_phage"],0,0,parameters_in["gamma_phage"],parameters_in["gamma_plasmid"],parameters_in["gamma_plasmid"],parameters_in["gamma_plasmid"],0,parameters_in["gamma_plasmid"]))
    fitness_costs = as.numeric(c(0,parameters_in["f_phage"],0,0,parameters_in["f_phage"],parameters_in["f_plasmid"],parameters_in["f_plasmid"],parameters_in["f_plasmid"],0,parameters_in["f_plasmid"]))
    growth_rate = as.numeric(parameters_in["grow"])
    rel_fit_human = as.numeric(parameters_in["rel_fit"])
    
  } else if(length(parameters_in) == 10){
    
    # If fixed input - same loss/gain rates different fitness
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    
    rate_loss = as.numeric(c(0,parameters_in[1],0,0,parameters_in[1],parameters_in[1],parameters_in[1],parameters_in[1],0,parameters_in[1]))
    rate_gain = as.numeric(c(0,parameters_in[2],0,0,parameters_in[2],parameters_in[2],parameters_in[2],parameters_in[2],0,parameters_in[2]))
    fitness_costs = as.numeric(c(0,parameters_in[3],0,0,parameters_in[4],parameters_in[5],parameters_in[6],parameters_in[7],0,parameters_in[8]))
    growth_rate = as.numeric(parameters_in[9])
    rel_fit_human = as.numeric(parameters_in[10])
    
  } else stop("Length of parameters_in does not match any possible option.")
  
  # Errors
  if(length(rate_gain) < 10){stop("Not enough gain parameters")}
  if(length(rate_loss) < 10){stop("Not enough loss parameters")}
  if(length(fitness_costs) < 10){stop("Not enough fitness parameters")}
  
  #parameters for growth
  #nb that's just a logistic deterministic function in the model at the moment
  Nmax = 0.5 * 1e7
  # Death rate: This captures immune system killing = 10% die initially 
  death_rate = 0.1
  
  # fitness cost - firstly without strain issues, then with human strain cost
  # Can be a sum of positive and negative values: fitness_costs_all could be negative (benefit)
  fitness_cost_all = rowSums(t(t(bacteria[,1:(ncol(bacteria)-2)]) * fitness_costs))
  fitness_cost_all = fitness_cost_all + c(rep(0, nrow(bacteria)/2),rep(1-rel_fit_human, nrow(bacteria)/2))
  
  #summary matrix to store number of bacteria in each strain at each timepoint
  all_results = matrix(0, nrow = tsteps, ncol = nrow(bacteria))
  all_results[1,] = bacteria[,"freq"]
  
  #summary matrix to store MGE prevalence at each timepoint in each parent
  all_mge_prev = matrix(0, nrow = 2*tsteps, ncol = (ncol(bacteria)-2))
  
  # Final check on input parameters  - upper bound on fitness
  if(any(fitness_cost_all > 20) || growth_rate > 3){
    print(c("Input error (fitness or growth)",fitness_cost_all))
    prev_predict = c()
    totl_predict = c()
  } else {
    
    for(t in 2:tsteps){
      
      # death randomly remove bacteria using rmultinom even before calculating MGE prev
      bacteria[,"freq"] = bacteria[,"freq"]-rmultinom(1, round(death_rate*sum(bacteria[,"freq"])), bacteria[,"freq"])

      #copy over the bacteria matrix to store new bacteria numbers as we go along
      new_bacteria = bacteria
      new_bacteria[,"freq"] = 0
      
      #we can already do some calculations here
      #total bacteria in environment currently:
      tot_bacteria = sum(bacteria[,"freq"])
      
      #MGE prevalence currently:
      #(remember first 10 columns are the MGE profiles)
      MGE_prevalence = bacteria[,c(1:10)]*bacteria[,"freq"]
      MGE_prevalence_all = Rfast::colsums(MGE_prevalence)/tot_bacteria
      
      # MGE prevalence in each parent strain
      if(sum(bacteria[1:(nrow(bacteria)/2),"freq"]) > 0){
        MGE_prevalence_1 = Rfast::colsums(MGE_prevalence[1:(nrow(bacteria)/2),])/sum(bacteria[1:(nrow(bacteria)/2),"freq"])
      } else { MGE_prevalence_1 = rep(0,10) }  # just in parent 1
      if(sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) > 0){
        MGE_prevalence_2 = Rfast::colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"])
      } else { MGE_prevalence_2 = rep(0,10) } # just in parent 2
      
      #Store MGE prev
      all_mge_prev[(t-1),] = MGE_prevalence_1
      all_mge_prev[(tsteps) + (t-1),] = MGE_prevalence_2
      
      #for each strain present
      present_strains = which(bacteria[,"freq"]>0)
      if(length(present_strains)>1) present_strains = sample(present_strains) # do sample to give random order otherwise prioritise parent 1 growth at stationary phase
      
      for(i in present_strains){
        
        #extract strain profile
        bacteria_profile = bacteria[i,-c(ncol(bacteria)-1, ncol(bacteria))]
        
        #loss proba is either the loss probability or 0 (if the MGE is already absent in that strain)
        lose_probas = pmin(rate_loss, bacteria_profile)
        
        #gain proba is either a density dependent proba (rate*n_recipient*n_donors/all_bacteria)
        #   or 0 (if MGE is already present in that strain)
        #print(c(rate_gain, bacteria[i,"freq"],"MGE",MGE_prevalence_all))
        gain_probas = 1-exp(-rate_gain*bacteria[i,"freq"]*MGE_prevalence_all) 
        gain_probas = pmin(gain_probas, 1-bacteria_profile)
        
        #sum probas (since each MGE can either be gained or lost only)
        #now, these probas represent the probability for each MGE status to change during this timestep
        #ie an MGE already present can only be lost, and an MGE currently absent can only be gained
        probas_sum = lose_probas + gain_probas
        
        #recover the difference matrix for our strain i
        differences = difference_list[[i]]
        
        #this aligns the probabilities with the other strains they correspond to
        #eg if there's a 0.5 proba to lose MGE 1, and that losing MGE 1 will change
        # the profile of strain i to the profile of strain 10, this will translate 
        # to a 0.5 proba of bacteria of strain i becoming strain 10 during this timestep
        #inversing the matrices twice is the fastest option to do this (I think?)
        probas = abs(t(t(differences[,-ncol(differences)]) * probas_sum))
        
        #collapse to vector using rowsums (only 1 value per row will be greater than 0)
        probas = Rfast::rowsums(probas)
        #okay, this is needed to set a probability of the strain NOT changing profile
        #it's not ideal, because sometimes the probabilities add up to more than 1, hence the max() function
        #something to look into...
        probas[Rfast::rowsums(differences[,1:10])==0] = max(0, 1 - sum(probas))
        #print(c("probas",probas,"gain",gain_probas))
        
        #use multinomial sampling to decide what the bacteria from strain i now become,
        # then add that amount to the updated matrix of bacteria numbers
        #here's where the "id" column in "differences" is useful: to align the indexing
        # between "differences" (which only contains valid transitions for strain i) and
        # "new_bacteria" (which contains all 2048 possible strains)
        probas[is.na(probas)] = 0 # got na errors in rmultinom - fix with this for now
        if(sum(probas) > 0){#print(parameters_in); break}
          new_bacteria[differences[,"id"],"freq"] = new_bacteria[differences[,"id"],"freq"] +
            rmultinom(1, bacteria[i, "freq"], probas) # due to discrete time step, don't want to end up with negative bugs as more than Nmax at some point
        } else {
          new_bacteria[i,"freq"] = bacteria[i,"freq"]
        }
        
        #if some bacteria remain in their original strain i, they now grow
        #currently just a deterministic logistic calculation
        # Need to add in fitness cost of elements: assume additive atm 
        # the use of min with new_bacteria is nice here, to essentially "stop" further growth once Nmax is reached
        new_bacteria[i,"freq"] = new_bacteria[i,"freq"] +
          min((Nmax - sum(new_bacteria[,"freq"])), round(bacteria[i,"freq"] * (max(0,(1 - fitness_cost_all[i]))) * growth_rate * (1 - tot_bacteria/Nmax)))
        
        #print(c(sum(bacteria[,"freq"]),round(bacteria[i,"freq"] * (1 - fitness_cost_all[i]) * growth_rate * (1 - sum(bacteria[,"freq"])/Nmax))))
      }
      
      #update main bacteria matrix
      bacteria = new_bacteria
      #store numbers for each strain at that timepoint
      all_results[t,] = bacteria[,"freq"]
      
    }
    all_results_out = all_results # store to check later
    
    tot_bacteria = sum(bacteria[,"freq"])
    
    # Last MGE prevalence
    MGE_prevalence = bacteria[,c(1:10)]*bacteria[,"freq"]
    MGE_prevalence_all = Rfast::colsums(MGE_prevalence)/tot_bacteria
    if(sum(bacteria[1:(nrow(bacteria)/2),"freq"]) > 0){
      MGE_prevalence_1 = Rfast::colsums(MGE_prevalence[1:(nrow(bacteria)/2),])/sum(bacteria[1:(nrow(bacteria)/2),"freq"])
    } else { MGE_prevalence_1 = rep(0,10) }  # just in parent 1
    if(sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) > 0){
      MGE_prevalence_2 = Rfast::colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"])
    } else { MGE_prevalence_2 = rep(0,10) } # just in parent 2
    
    #Store MGE prev
    all_mge_prev[(tsteps),] = MGE_prevalence_1
    all_mge_prev[(2*tsteps),] = MGE_prevalence_2
    
    # are these steps really necessary, as well as the conversion to dataframes?
    # sticking to matrix would be best, and indexing can be recovered based on tsteps and baseline row number
    
    #clean up 
    all_results = as.data.frame(all_results_out)
    colnames(all_results) = c(1:2048)
    all_results$time = c(1:nrow(all_results))
    all_results = reshape2::melt(all_results, id.vars = "time")
    all_results$parent = c(rep(1, tsteps * nrow(bacteria)/2), rep(2, tsteps * nrow(bacteria)/2))
    
    # Totals - want prevalence of each parent strains
    totl_predict = all_results %>%
      select(time,value, parent) %>%
      filter(time %in% c(4,48,96,288,384)) %>%
      group_by(time, parent) %>% 
      summarise(total = sum(value),.groups = "drop") 
    
    # Prevalence of MGE at same time as data
    #clean up MGE
    all_mge_prev = as.data.frame(all_mge_prev)
    all_mge_prev$time = c(1:tsteps)
    all_mge_prev$parent = c(rep(1, tsteps), rep(2, tsteps))
    all_mge_prev = reshape2::melt(all_mge_prev, id.vars = c("parent","time"))
    
    prev_predict = all_mge_prev %>% 
      filter(time %in% c(4,48,96,288,384))
  }
  
  # Output
  return(list(all_results = all_results,
              prev_predict = prev_predict, totl_predict = totl_predict))
}
