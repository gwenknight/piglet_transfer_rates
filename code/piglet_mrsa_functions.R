#### code from Quentin - using transition matrix 

### Initial population set up 

initial_piglet_setup <- function(tsteps){
  #matrix to store all possible MGE combinations (ie all possible strains)
  #MGE1 presence / MGE2 presence / MGE3 presence / ... / freq / parent
  bacteria = cbind(rbind(as.matrix(expand.grid(rep(list(0:1), 10))),
                         as.matrix(expand.grid(rep(list(0:1), 10)))),0)
  
  # Initial conditions
  Pinit = t(c(1,1,1,1,1,1,1,1,0,0,0)) # pig adapted
  Qinit = t(c(0,0,0,0,0,0,0,0,1,1,0)) # human adapted
  pw <- which(apply(bacteria, 1, function(x) return(all(x == Pinit))))
  qw <- which(apply(bacteria, 1, function(x) return(all(x == Qinit))))
  
  # min and max pw / qw as same row exists in parent 1 and parent 2
  bacteria[min(pw),ncol(bacteria)] <- 7 * 10^2 # start with all at P1
  bacteria[max(qw),ncol(bacteria)] <- 7 * 10^2 # start with all at Q1
  
  # Add in parent label
  bacteria = cbind(bacteria, c(rep(1, nrow(bacteria)/2), rep(2, nrow(bacteria)/2)))
  
  #set column names for easier reference later
  colnames(bacteria)[c(ncol(bacteria)-1, ncol(bacteria))] = c("freq", "parent")
  
  #first, we work out valid transitions for each strain
  # since we are assuming that only one gain/loss event can happen per timestep
  difference_list = list()
  
  for(i in 1:nrow(bacteria)){
    
    #valid transitions must have the same parent type
    same_parents = which(bacteria[,"parent"] == bacteria[i,"parent"])
    
    #calculate the difference between each strain and every other strain
    differences = t(t(bacteria[,-c(ncol(bacteria)-1, ncol(bacteria))]) -
                      bacteria[i,-c(ncol(bacteria)-1, ncol(bacteria))])
    
    #only keep transitions with at most 1 difference in MGE profile with the strain we're currently looking at
    possible_transitions = which(rowsums(abs(differences)) < 2,)
    
    #add an "id" variable to use later for indexing
    differences = cbind(differences, "id" = c(1:nrow(differences)))
    
    #final possible transitions are those with same parent AND at most 1 difference in MGE profile
    possible_transitions = intersect(same_parents, possible_transitions)
    
    #only keep the section of the difference matrix that's for valid transitions
    differences = differences[possible_transitions,]
    
    #store in list
    difference_list[[i]] = differences
    
  }
  
  return(list(difference_list = difference_list, bacteria = bacteria)) 
}


### Function to move MGE
piglet_mrsa_movement_old <- function(tsteps, parameters_in){
  # tsteps = number of time steps
  # parameters_in = parameters needed to run model
  # bacteria = initial distribution and all possible profiles
  # difference_list = possible transitions
  
  # Assign parameter values
  # If have rates for all 10 elements
  if(length(parameters_in) == 32){
    rate_loss <- as.numeric(parameters_in[1:10])
    rate_gain <- as.numeric(parameters_in[11:20] )
    fitness_costs <- as.numeric(parameters_in[21:30])
    growth_rate <- as.numeric(parameters_in[31])
    rel_fit_human <- as.numeric(parameters_in[32])
  }
  
  # If just looking at the elements that move
  if(length(parameters_in) < 32 && length(parameters_in) > 8){
    rate_loss <- as.numeric(c(0,parameters_in["mu2"],0,0,parameters_in["mu5"],parameters_in["mu6"],parameters_in["mu7"],parameters_in["mu8"],0,parameters_in["mu10"]))
    rate_gain <- as.numeric(c(0,parameters_in["gamma2"],0,0,parameters_in["gamma5"],parameters_in["gamma6"],parameters_in["gamma7"],parameters_in["gamma8"],0,parameters_in["gamma10"]))
    fitness_costs <- as.numeric(c(0,parameters_in["f2"],0,0,parameters_in["f5"],parameters_in["f6"],parameters_in["f7"],parameters_in["f8"],0,parameters_in["f10"]))
    growth_rate <- as.numeric(parameters_in["grow"])
    rel_fit_human <- as.numeric(parameters_in["rel_fit"])
  }
  
  # If fixed input - same rates for all elements
  if(length(parameters_in) == 5){
    if(is.null(names(parameters_in))){stop("No names for parameters_in")}
    rate_loss <- as.numeric(c(0,parameters_in["mu"],0,0,parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],0,parameters_in["mu"]))
    rate_gain <- as.numeric(c(0,parameters_in["gamma"],0,0,parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],0,parameters_in["gamma"]))
    fitness_costs <- as.numeric(c(0,parameters_in["f"],0,0,parameters_in["f"],parameters_in["f"],parameters_in["f"],parameters_in["f"],0,parameters_in["f"]))
    growth_rate <- as.numeric(parameters_in["grow"])
    rel_fit_human <- as.numeric(parameters_in["rel_fit"])
  }
  
  # If fixed input - same rates for phage vs plasmids
  if(length(parameters_in) == 8){
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    rate_loss <- as.numeric(c(0,parameters_in["mu_phage"],0,0,parameters_in["mu_phage"],parameters_in["mu_plasmid"],parameters_in["mu_plasmid"],parameters_in["mu_plasmid"],0,parameters_in["mu_plasmid"]))
    rate_gain <- as.numeric(c(0,parameters_in["gamma_phage"],0,0,parameters_in["gamma_phage"],parameters_in["gamma_plasmid"],parameters_in["gamma_plasmid"],parameters_in["gamma_plasmid"],0,parameters_in["gamma_plasmid"]))
    fitness_costs <- as.numeric(c(0,parameters_in["f_phage"],0,0,parameters_in["f_phage"],parameters_in["f_plasmid"],parameters_in["f_plasmid"],parameters_in["f_plasmid"],0,parameters_in["f_plasmid"]))
    growth_rate <- as.numeric(parameters_in["grow"])
    rel_fit_human <- as.numeric(parameters_in["rel_fit"])
  }
  
  # Errors
  if(length(rate_gain) < 10){stop("Not enough gain parameters")}
  if(length(rate_loss) < 10){stop("Not enough loss parameters")}
  if(length(fitness_costs) < 10){stop("Not enough fitness parameters")}
  
  #matrix to store all possible MGE combinations (ie all possible strains)
  #MGE1 presence / MGE2 presence / MGE3 presence / ... / freq / parent
  bacteria = cbind(rbind(as.matrix(expand.grid(rep(list(0:1), 10))),
                         as.matrix(expand.grid(rep(list(0:1), 10)))),0)
  
  # Initial conditions
  Pinit = t(c(0,1,0,0,1,1,1,1,0,0,0))
  Qinit = t(c(0,0,0,0,0,0,0,0,0,1,0))
  pw <- which(apply(bacteria, 1, function(x) return(all(x == Pinit))))
  qw <- which(apply(bacteria, 1, function(x) return(all(x == Qinit))))
  
  # min and max pw / qw as same row exists in parent 1 and parent 2
  bacteria[min(pw),ncol(bacteria)] <- 7 * 10^2 # start with all at P1
  bacteria[max(qw),ncol(bacteria)] <- 7 * 10^2 # start with all at Q1
  
  # Add in parent label
  bacteria = cbind(bacteria, c(rep(1, nrow(bacteria)/2), rep(2, nrow(bacteria)/2)))
  
  #set column names for easier reference later
  colnames(bacteria)[c(ncol(bacteria)-1, ncol(bacteria))] = c("freq", "parent")
  
  #first, we work out valid transitions for each strain
  # since we are assuming that only one gain/loss event can happen per timestep
  difference_list = list()
  
  for(i in 1:nrow(bacteria)){
    
    #valid transitions must have the same parent type
    same_parents = which(bacteria[,"parent"] == bacteria[i,"parent"])
    
    #calculate the difference between each strain and every other strain
    differences = t(t(bacteria[,-c(ncol(bacteria)-1, ncol(bacteria))]) -
                      bacteria[i,-c(ncol(bacteria)-1, ncol(bacteria))])
    
    #only keep transitions with at most 1 difference in MGE profile with the strain we're currently looking at
    possible_transitions = which(rowsums(abs(differences)) < 2,)
    
    #add an "id" variable to use later for indexing
    differences = cbind(differences, "id" = c(1:nrow(differences)))
    
    #final possible transitions are those with same parent AND at most 1 difference in MGE profile
    possible_transitions = intersect(same_parents, possible_transitions)
    
    #only keep the section of the difference matrix that's for valid transitions
    differences = differences[possible_transitions,]
    
    #store in list
    difference_list[[i]] = differences
    
  }
  
  #parameters for growth
  #nb that's just a logistic deterministic function in the model at the moment
  Nmax = 1e7
  # fitness cost
  fitness_cost_all <- rowSums(t(t(bacteria[,1:(ncol(bacteria)-2)]) * fitness_costs))
  
  #summary matrix to store number of bacteria in each strain at each timepoint
  all_results = matrix(0, nrow = tsteps, ncol = nrow(bacteria))
  all_results[1,] = bacteria[,"freq"]
  
  #summary matrix to store MGE prevalence at each timepoint in each parent
  all_mge_prev = matrix(0, nrow = 2*tsteps, ncol = (ncol(bacteria)-2))
  
  # Final check on input parameters 
  if(any(fitness_cost_all > 1)){print("Fitness cost too large");  prev_predict = c(); totl_predict = c(); }else{
    
    for(t in 2:tsteps){
      
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
      MGE_prevalence_1 = Rfast::colsums(MGE_prevalence[1:nrow(bacteria)/2,])/sum(bacteria[1:nrow(bacteria)/2,"freq"]) # just in parent 1
      MGE_prevalence_2 = Rfast::colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) # just in parent 2
      
      #Store MGE prev
      all_mge_prev[(t-1),] = MGE_prevalence_1
      all_mge_prev[(1+tsteps) + (t-1),] = MGE_prevalence_2
      
      #for each strain present
      for(i in 1:nrow(bacteria)){
        
        #skip if the number of bacteria of that strain is 0
        if(bacteria[i,"freq"] == 0) next
        
        #extract strain profile
        bacteria_profile = bacteria[i,-c(ncol(bacteria)-1, ncol(bacteria))]
        
        #loss proba is either the loss probability or 0 (if the MGE is already absent in that strain)
        lose_probas = pmin(rate_loss, bacteria_profile)
        
        #gain proba is either a density dependent proba (rate*n_recipient*n_donors/all_bacteria)
        #   or 0 (if MGE is already present in that strain)
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
        probas = rowsums(probas)
        #okay, this is needed to set a probability of the strain NOT changing profile
        #it's not ideal, because sometimes the probabilities add up to more than 1, hence the max() function
        #something to look into...
        probas[rowsums(differences[,1:10])==0] = max(0, 1 - probas)
        
        #use multinomial sampling to decide what the bacteria from strain i now become,
        # then add that amount to the updated matrix of bacteria numbers
        #here's where the "id" column in "differences" is useful: to align the indexing
        # between "differences" (which only contains valid transitions for strain i) and
        # "new_bacteria" (which contains all 2048 possible strains)
        new_bacteria[differences[,"id"],"freq"] = new_bacteria[differences[,"id"],"freq"] +
          rmultinom(1, bacteria[i, "freq"], probas)
        
        #if some bacteria remain in their original strain i, they now grow
        #currently just a deterministic logistic calculation
        # Need to add in fitness cost of elements: assume additive atm 
        new_bacteria[i,"freq"] = new_bacteria[i,"freq"] +
          round(bacteria[i,"freq"] * (1 - fitness_cost_all[i]) * growth_rate * (1 - sum(bacteria[,"freq"])/Nmax))
      }
      
      #update main bacteria matrix
      bacteria = new_bacteria
      #store numbers for each strain at that timepoint
      all_results[t,] = bacteria[,"freq"]
      
    }
    all_results_out <- all_results # store to check later
    
    # Last MGE prevalence
    MGE_prevalence = bacteria[,c(1:10)]*bacteria[,"freq"]
    MGE_prevalence_all = colsums(MGE_prevalence)/tot_bacteria
    MGE_prevalence_1 = colsums(MGE_prevalence[1:nrow(bacteria)/2,])/sum(bacteria[1:nrow(bacteria)/2,"freq"]) # just in parent 1
    MGE_prevalence_2 = colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) # just in parent 2
    
    #Store MGE prev
    all_mge_prev[(tsteps),] = MGE_prevalence_1
    all_mge_prev[(2*tsteps),] = MGE_prevalence_2
    
    #clean up 
    all_results = as.data.frame(all_results_out)
    all_results$time = c(1:nrow(all_results))
    all_results = reshape2::melt(all_results, id.vars = "time")
    all_results$parent = c(rep(1, tsteps * nrow(bacteria)/2), rep(2, tsteps * nrow(bacteria)/2))
    
    # Totals - want prevalence of each parent strains
    totl_predict<- all_results %>% select(time,value, parent) %>%
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time, parent) %>% 
      summarise(total = sum(value),.groups = "drop") 
    
    # Prevalence of MGE at same time as data
    #clean up MGE
    all_mge_prev = as.data.frame(all_mge_prev)
    all_mge_prev$time = seq(1,tsteps,1)
    all_mge_prev$parent = c(rep(1, tsteps), rep(2, tsteps))
    all_mge_prev = reshape2::melt(all_mge_prev, id.vars = c("parent","time"))
    
    prev_predict <- all_mge_prev %>% 
      filter(time %in% c(4,48,96,288,384))
  }
  # Output
  return(list(all_results = all_results, 
              prev_predict = prev_predict, totl_predict = totl_predict))
}

piglet_mrsa_movement <- function(tsteps, parameters_in, bacteria, difference_list){
  # tsteps = number of time steps
  # parameters_in = parameters needed to run model
  # bacteria = initial distribution and all possible profiles
  # difference_list = possible transitions
  
  # Assign parameter values
  # If have rates for all 10 elements
  if(length(parameters_in) == 32){
    rate_loss <- as.numeric(parameters_in[1:10])
    rate_gain <- as.numeric(parameters_in[11:20] )
    fitness_costs <- as.numeric(parameters_in[21:30])
    growth_rate <- as.numeric(parameters_in[31])
    rel_fit_human <- as.numeric(parameters_in[32])
  }
  
  # If just looking at the elements that move
  if(length(parameters_in) < 32 && length(parameters_in) > 8){
    rate_loss <- as.numeric(c(0,parameters_in["mu2"],0,0,parameters_in["mu5"],parameters_in["mu6"],parameters_in["mu7"],parameters_in["mu8"],0,parameters_in["mu10"]))
    rate_gain <- as.numeric(c(0,parameters_in["gamma2"],0,0,parameters_in["gamma5"],parameters_in["gamma6"],parameters_in["gamma7"],parameters_in["gamma8"],0,parameters_in["gamma10"]))
    fitness_costs <- as.numeric(c(0,parameters_in["f2"],0,0,parameters_in["f5"],parameters_in["f6"],parameters_in["f7"],parameters_in["f8"],0,parameters_in["f10"]))
    growth_rate <- as.numeric(parameters_in["grow"])
    rel_fit_human <- as.numeric(parameters_in["rel_fit"])
  }
  
  # If fixed input - same rates for all elements
  if(length(parameters_in) == 5){
    if(is.null(names(parameters_in))){stop("No names for parameters_in")}
    rate_loss <- as.numeric(c(0,parameters_in["mu"],0,0,parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],0,parameters_in["mu"]))
    rate_gain <- as.numeric(c(0,parameters_in["gamma"],0,0,parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],0,parameters_in["gamma"]))
    fitness_costs <- as.numeric(c(0,parameters_in["f"],0,0,parameters_in["f"],parameters_in["f"],parameters_in["f"],parameters_in["f"],0,parameters_in["f"]))
    growth_rate <- as.numeric(parameters_in["grow"])
    rel_fit_human <- as.numeric(parameters_in["rel_fit"])
  }
  
  # If fixed input - same rates for all elements and no fitness cost 
  if(length(parameters_in) == 4){
    if(is.null(names(parameters_in))){stop("No names for parameters_in")}
    rate_loss <- as.numeric(c(0,parameters_in["mu"],0,0,parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],0,parameters_in["mu"]))
    rate_gain <- as.numeric(c(0,parameters_in["gamma"],0,0,parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],0,parameters_in["gamma"]))
    fitness_costs <- as.numeric(c(0,0,0,0,0,0,0,0,0,0))
    growth_rate <- as.numeric(parameters_in["grow"])
    rel_fit_human <- as.numeric(parameters_in["rel_fit"])
  }
  
  # If fixed input - same rates for phage vs plasmids
  if(length(parameters_in) == 8){
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    rate_loss <- as.numeric(c(0,parameters_in["mu_phage"],0,0,parameters_in["mu_phage"],parameters_in["mu_plasmid"],parameters_in["mu_plasmid"],parameters_in["mu_plasmid"],0,parameters_in["mu_plasmid"]))
    rate_gain <- as.numeric(c(0,parameters_in["gamma_phage"],0,0,parameters_in["gamma_phage"],parameters_in["gamma_plasmid"],parameters_in["gamma_plasmid"],parameters_in["gamma_plasmid"],0,parameters_in["gamma_plasmid"]))
    fitness_costs <- as.numeric(c(0,parameters_in["f_phage"],0,0,parameters_in["f_phage"],parameters_in["f_plasmid"],parameters_in["f_plasmid"],parameters_in["f_plasmid"],0,parameters_in["f_plasmid"]))
    growth_rate <- as.numeric(parameters_in["grow"])
    rel_fit_human <- as.numeric(parameters_in["rel_fit"])
  }
  
  # If fixed input - same loss/gain rates different fitness
  if(length(parameters_in) == 10){
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    rate_loss <- as.numeric(c(0,parameters_in["mu"],0,0,parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],parameters_in["mu"],0,parameters_in["mu"]))
    rate_gain <- as.numeric(c(0,parameters_in["gamma"],0,0,parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],parameters_in["gamma"],0,parameters_in["gamma"]))
    fitness_costs <- as.numeric(c(0,parameters_in["f2"],0,0,parameters_in["f5"],parameters_in["f6"],parameters_in["f7"],parameters_in["f8"],0,parameters_in["f10"]))
    growth_rate <- as.numeric(parameters_in["grow"])
    rel_fit_human <- as.numeric(parameters_in["rel_fit"])
  }
  
  # Errors
  if(length(rate_gain) < 10){stop("Not enough gain parameters")}
  if(length(rate_loss) < 10){stop("Not enough loss parameters")}
  if(length(fitness_costs) < 10){stop("Not enough fitness parameters")}
  
  #parameters for growth
  #nb that's just a logistic deterministic function in the model at the moment
  Nmax = 0.5 * 1e7
  # Death rate: timestep currently 1 hr, replicate in 20min in optimal... this captures immune system killing? 10% die? 
  death_rate <- 0.1
  
  # fitness cost - firstly without strain issues, then with human strain cost
  fitness_cost_all_noh <- rowSums(t(t(bacteria[,1:(ncol(bacteria)-2)]) * fitness_costs))
  fitness_cost_all <- fitness_cost_all_noh * c(rep(1, nrow(bacteria)/2),rep(rel_fit_human, nrow(bacteria)/2))
  
  #summary matrix to store number of bacteria in each strain at each timepoint
  all_results = matrix(0, nrow = tsteps, ncol = nrow(bacteria))
  all_results[1,] = bacteria[,"freq"]
  
  #summary matrix to store MGE prevalence at each timepoint in each parent
  all_mge_prev = matrix(0, nrow = 2*tsteps, ncol = (ncol(bacteria)-2))
  
  # Final check on input parameters 
  if(any(fitness_cost_all > 20) || growth_rate > 2.94){print(c("Input error (fitness or growth)",fitness_cost_all));  prev_predict = c(); totl_predict = c(); }else{
    
    for(t in 2:tsteps){
      
      #copy over the bacteria matrix to store new bacteria numbers as we go along
      new_bacteria = bacteria
      new_bacteria[,"freq"] = 0
      
      #we can already do some calculations here
      #total bacteria in environment currently:
      tot_bacteria = sum(bacteria[,"freq"])
      
      #MGE prevalence currently:
      #(remember first 10 columns are the MGE profiles)
      MGE_prevalence = bacteria[,c(1:10)]*bacteria[,"freq"]
      MGE_prevalence_all = colsums(MGE_prevalence)/tot_bacteria
      #MGE_prevalence_1 = colsums(MGE_prevalence[1:nrow(bacteria)/2,])/sum(bacteria[1:nrow(bacteria)/2,"freq"]) # just in parent 1
      #MGE_prevalence_2 = colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) # just in parent 2
      
      if(sum(bacteria[1:nrow(bacteria)/2,"freq"]) > 0){ MGE_prevalence_1 = colsums(MGE_prevalence[1:nrow(bacteria)/2,])/sum(bacteria[1:nrow(bacteria)/2,"freq"])}else{MGE_prevalence_1 = rep(0,10)}  # just in parent 1
      if(sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) > 0){MGE_prevalence_2 = colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"])}else{MGE_prevalence_2=rep(0,10)} # just in parent 2
      
      
      #Store MGE prev
      all_mge_prev[(t-1),] = MGE_prevalence_1
      all_mge_prev[(1+tsteps) + (t-1),] = MGE_prevalence_2
      
      #for each strain present
      present_strains <- sample(which(bacteria[,"freq"]>0)) # do sample to give random order otherwise prioritise parent 1 growth at stationary phase
      
      for(i in present_strains){
        
        # In order for new profiles to appear at the start of the timestep remove
        # some of the existing strains - death - which will be replaced if lots of these remain
        # Background constant: relative to growth rate which is fitted 
        new_bacteria[i,"freq"] <- (1 - death_rate) * new_bacteria[i,"freq"] 
        
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
        probas = rowsums(probas)
        #okay, this is needed to set a probability of the strain NOT changing profile
        #it's not ideal, because sometimes the probabilities add up to more than 1, hence the max() function
        #something to look into...
        probas[rowsums(differences[,1:10])==0] = max(0, 1 - probas)
        #print(c("probas",probas,"gain",gain_probas))
        
        #use multinomial sampling to decide what the bacteria from strain i now become,
        # then add that amount to the updated matrix of bacteria numbers
        #here's where the "id" column in "differences" is useful: to align the indexing
        # between "differences" (which only contains valid transitions for strain i) and
        # "new_bacteria" (which contains all 2048 possible strains)
        probas[is.na(probas)] <- 0 # got na errors in rmultinom - fix with this for now
        #if(sum(probas) == "NA"){print(parameters_in); break}
        
        new_bacteria[differences[,"id"],"freq"] = new_bacteria[differences[,"id"],"freq"] +
          pmin((Nmax - sum(new_bacteria[,"freq"]))/length(new_bacteria[i,]), rmultinom(1, bacteria[i, "freq"], probas)) # due to discrete time step, don't want to end up with negative bugs as more than Nmax at some point
        
        #if some bacteria remain in their original strain i, they now grow
        #currently just a deterministic logistic calculation
        # Need to add in fitness cost of elements: assume additive atm 
        new_bacteria[i,"freq"] = new_bacteria[i,"freq"] +
          min((Nmax - sum(new_bacteria[,"freq"])), round(bacteria[i,"freq"] * (1 - fitness_cost_all[i]) * growth_rate * (1 - sum(bacteria[,"freq"])/Nmax)))
        
        #print(c(sum(bacteria[,"freq"]),round(bacteria[i,"freq"] * (1 - fitness_cost_all[i]) * growth_rate * (1 - sum(bacteria[,"freq"])/Nmax))))
      }
      
      #update main bacteria matrix
      bacteria = new_bacteria
      #store numbers for each strain at that timepoint
      all_results[t,] = bacteria[,"freq"]
      
    }
    all_results_out <- all_results # store to check later
    
    # Last MGE prevalence
    MGE_prevalence = bacteria[,c(1:10)]*bacteria[,"freq"]
    MGE_prevalence_all = colsums(MGE_prevalence)/tot_bacteria
    if(sum(bacteria[1:nrow(bacteria)/2,"freq"]) > 0){ MGE_prevalence_1 = colsums(MGE_prevalence[1:nrow(bacteria)/2,])/sum(bacteria[1:nrow(bacteria)/2,"freq"])}else{MGE_prevalence_1 = rep(0,10)}  # just in parent 1
    if(sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) > 0){MGE_prevalence_2 = colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"])}else{MGE_prevalence_2=rep(0,10)} # just in parent 2
    
    #Store MGE prev
    all_mge_prev[(tsteps),] = MGE_prevalence_1
    all_mge_prev[(2*tsteps),] = MGE_prevalence_2
    
    #clean up 
    all_results = as.data.frame(all_results_out)
    colnames(all_results) <- seq(1:2048)
    all_results$time = c(1:nrow(all_results))
    all_results = reshape2::melt(all_results, id.vars = "time")
    all_results$parent = c(rep(1, tsteps * nrow(bacteria)/2), rep(2, tsteps * nrow(bacteria)/2))
    
    # Totals - want prevalence of each parent strains
    totl_predict<- all_results %>% select(time,value, parent) %>%
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time, parent) %>% 
      summarise(total = sum(value),.groups = "drop") 
    
    # Prevalence of MGE at same time as data
    #clean up MGE
    all_mge_prev = as.data.frame(all_mge_prev)
    all_mge_prev$time = seq(1,tsteps,1)
    all_mge_prev$parent = c(rep(1, tsteps), rep(2, tsteps))
    all_mge_prev = reshape2::melt(all_mge_prev, id.vars = c("parent","time"))
    
    prev_predict <- all_mge_prev %>% 
      filter(time %in% c(4,48,96,288,384))
  }
  
  # Output
  return(list(all_results = all_results,
              prev_predict = prev_predict, totl_predict = totl_predict))
}



### Calculate log posterior from model at theta parameter values
run_sim_logPosterior <- function(theta_in){
  tsteps = 384
  
  ## Run for these parameters
  out <- piglet_mrsa_movement(tsteps, theta_in, ini$bacteria, ini$difference_list)
  
  profile_end1.1 <- c(1,1,1,1,1,1,1,1,0,0,1) # pig - same as at start
  profile_end2.1 <- c(0,1,0,0,1,0,0,0,1,0,2) # human - gains phi6 / phi2 / loses p4 (2/5/10)
  profile_end2.2 <- c(0,1,0,0,1,0,1,0,1,0,2) # human
  profiles_needed_end <- c(row.match(profile_end1.1, as.data.frame(ini$bacteria)%>%select(-freq), nomatch = NA),
                           row.match(profile_end2.1, as.data.frame(ini$bacteria)%>%select(-freq), nomatch = NA),
                           row.match(profile_end2.2, as.data.frame(ini$bacteria)%>%select(-freq), nomatch = NA))
  
  ### Likelihood
  if(!is.null(out$prev_predict)){
    
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
      mutate(exp_miss = dnorm(prob_all, 1, 0.1)) %>% # experimental measure - 10% around a prob of 1
      summarise(sum(log(exp_miss))) #summarise(sum(log(prob_all)))
    
    #### total bugs output from model (b)
    model_outputt <- out$totl_predict
    
    likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% 
      #mutate(val_in = as.numeric(between(total,min,max))) %>% # instead of 1 / 0 make distance 
      mutate(val_in = dnorm(total, mean = (max - min)/2 - min, sd = 2 * (max - min))) %>% # instead of 1 / 0 make distance 
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
    likelihood_profile_end <- sum(log(dnorm(prof_end,mean = c(0.8,0.2,0.2), sd = 0.5)))
    
    # Don't make -Inf possible
    # if(nrow(profile_end) == 3 && total_end[2] > 0){
    #   likelihood_profile_end <- profile_end %>% 
    #     mutate(prop = value / c(total_end[1],total_end[2],total_end[2]),
    #            cutoff = pmax(0,prop - 0.05),#c(0.8,0.2,0.2)), # more the better # not using -- too strict
    #            likelihood = log(prop)) %>% summarise(sum(likelihood))} else{likelihood_profile_end <- -Inf}
    
    #### Compare to data 
    compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals + 10 * likelihood_profile_end # add in a 10* weight for profile_end as otherwise only a small contribution relatively
  }else{compare_dat <- -Inf}
  
  # return log likelihood
  as.numeric(compare_dat)
}

### Calculate log posterior from model at theta parameter values: for LaplaceDemon

run_sim_logPosterior_para <- function(theta_in){
  library(tidyverse)
  
  tsteps = 384
  
  ## Run for these parameters
  # Instead of this: out <- piglet_mrsa_movement(tsteps, theta_in, ini$bacteria, ini$difference_list)
  # put function below... 
  
  # tsteps = number of time steps
  # theta_in = parameters needed to run model
  # bacteria = initial distribution and all possible profiles
  # difference_list = possible transitions
  
  ### HAVE TO PUT INITIAL IN HERE
  bacteria = cbind(rbind(as.matrix(expand.grid(rep(list(0:1), 10))),
                         as.matrix(expand.grid(rep(list(0:1), 10)))),0)
  
  # Initial conditions
  Pinit = t(c(1,1,1,1,1,1,1,1,0,0,0)) # pig adapted
  Qinit = t(c(0,0,0,0,0,0,0,0,1,1,0)) # human adapted
  pw <- which(apply(bacteria, 1, function(x) return(all(x == Pinit))))
  qw <- which(apply(bacteria, 1, function(x) return(all(x == Qinit))))
  
  # min and max pw / qw as same row exists in parent 1 and parent 2
  bacteria[min(pw),ncol(bacteria)] <- 7 * 10^2 # start with all at P1
  bacteria[max(qw),ncol(bacteria)] <- 7 * 10^2 # start with all at Q1
  
  # Add in parent label
  bacteria = cbind(bacteria, c(rep(1, nrow(bacteria)/2), rep(2, nrow(bacteria)/2)))
  
  #set column names for easier reference later
  colnames(bacteria)[c(ncol(bacteria)-1, ncol(bacteria))] = c("freq", "parent")
  
  #first, we work out valid transitions for each strain
  # since we are assuming that only one gain/loss event can happen per timestep
  difference_list = list()
  
  for(i in 1:nrow(bacteria)){
    
    #valid transitions must have the same parent type
    same_parents = which(bacteria[,"parent"] == bacteria[i,"parent"])
    
    #calculate the difference between each strain and every other strain
    differences = t(t(bacteria[,-c(ncol(bacteria)-1, ncol(bacteria))]) -
                      bacteria[i,-c(ncol(bacteria)-1, ncol(bacteria))])
    
    #only keep transitions with at most 1 difference in MGE profile with the strain we're currently looking at
    possible_transitions = which(Rfast::rowsums(abs(differences)) < 2,)
    
    #add an "id" variable to use later for indexing
    differences = cbind(differences, "id" = c(1:nrow(differences)))
    
    #final possible transitions are those with same parent AND at most 1 difference in MGE profile
    possible_transitions = intersect(same_parents, possible_transitions)
    
    #only keep the section of the difference matrix that's for valid transitions
    differences = differences[possible_transitions,]
    
    #store in list
    difference_list[[i]] = differences
    
  }
  
  # Assign parameter values
  # If have rates for all 10 elements
  if(length(theta_in) == 32){
    rate_loss <- as.numeric(theta_in[1:10])
    rate_gain <- as.numeric(theta_in[11:20] )
    fitness_costs <- as.numeric(theta_in[21:30])
    growth_rate <- as.numeric(theta_in[31])
    rel_fit_human <- as.numeric(theta_in[32])
  }
  
  # If just looking at the elements that move
  if(length(theta_in) < 32 && length(theta_in) > 8){
    rate_loss <- as.numeric(c(0,theta_in["mu2"],0,0,theta_in["mu5"],theta_in["mu6"],theta_in["mu7"],theta_in["mu8"],0,theta_in["mu10"]))
    rate_gain <- as.numeric(c(0,theta_in["gamma2"],0,0,theta_in["gamma5"],theta_in["gamma6"],theta_in["gamma7"],theta_in["gamma8"],0,theta_in["gamma10"]))
    fitness_costs <- as.numeric(c(0,theta_in["f2"],0,0,theta_in["f5"],theta_in["f6"],theta_in["f7"],theta_in["f8"],0,theta_in["f10"]))
    growth_rate <- as.numeric(theta_in["grow"])
    rel_fit_human <- as.numeric(theta_in["rel_fit"])
  }
  
  # If fixed input - same rates for all elements
  if(length(theta_in) == 5){
    if(is.null(names(theta_in))){stop("No names for theta_in")}
    rate_loss <- as.numeric(c(0,theta_in["mu"],0,0,theta_in["mu"],theta_in["mu"],theta_in["mu"],theta_in["mu"],0,theta_in["mu"]))
    rate_gain <- as.numeric(c(0,theta_in["gamma"],0,0,theta_in["gamma"],theta_in["gamma"],theta_in["gamma"],theta_in["gamma"],0,theta_in["gamma"]))
    fitness_costs <- as.numeric(c(0,theta_in["f"],0,0,theta_in["f"],theta_in["f"],theta_in["f"],theta_in["f"],0,theta_in["f"]))
    growth_rate <- as.numeric(theta_in["grow"])
    rel_fit_human <- as.numeric(theta_in["rel_fit"])
  }
  
  # If fixed input - same rates for all elements and no fitness cost 
  if(length(theta_in) == 4){
    if(is.null(names(theta_in))){stop("No names for theta_in")}
    rate_loss <- as.numeric(c(0,theta_in["mu"],0,0,theta_in["mu"],theta_in["mu"],theta_in["mu"],theta_in["mu"],0,theta_in["mu"]))
    rate_gain <- as.numeric(c(0,theta_in["gamma"],0,0,theta_in["gamma"],theta_in["gamma"],theta_in["gamma"],theta_in["gamma"],0,theta_in["gamma"]))
    fitness_costs <- as.numeric(c(0,0,0,0,0,0,0,0,0,0))
    growth_rate <- as.numeric(theta_in["grow"])
    rel_fit_human <- as.numeric(theta_in["rel_fit"])
  }
  
  # If fixed input - same rates for phage vs plasmids
  if(length(theta_in) == 8){
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    rate_loss <- as.numeric(c(0,theta_in["mu_phage"],0,0,theta_in["mu_phage"],theta_in["mu_plasmid"],theta_in["mu_plasmid"],theta_in["mu_plasmid"],0,theta_in["mu_plasmid"]))
    rate_gain <- as.numeric(c(0,theta_in["gamma_phage"],0,0,theta_in["gamma_phage"],theta_in["gamma_plasmid"],theta_in["gamma_plasmid"],theta_in["gamma_plasmid"],0,theta_in["gamma_plasmid"]))
    fitness_costs <- as.numeric(c(0,theta_in["f_phage"],0,0,theta_in["f_phage"],theta_in["f_plasmid"],theta_in["f_plasmid"],theta_in["f_plasmid"],0,theta_in["f_plasmid"]))
    growth_rate <- as.numeric(theta_in["grow"])
    rel_fit_human <- as.numeric(theta_in["rel_fit"])
  }
  
  # If fixed input - same loss/gain rates different fitness
  if(length(theta_in) == 10){
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    rate_loss <- as.numeric(c(0,theta_in["mu"],0,0,theta_in["mu"],theta_in["mu"],theta_in["mu"],theta_in["mu"],0,theta_in["mu"]))
    rate_gain <- as.numeric(c(0,theta_in["gamma"],0,0,theta_in["gamma"],theta_in["gamma"],theta_in["gamma"],theta_in["gamma"],0,theta_in["gamma"]))
    fitness_costs <- as.numeric(c(0,theta_in["f2"],0,0,theta_in["f5"],theta_in["f6"],theta_in["f7"],theta_in["f8"],0,theta_in["f10"]))
    growth_rate <- as.numeric(theta_in["grow"])
    rel_fit_human <- as.numeric(theta_in["rel_fit"])
  }
  
  # Errors
  if(length(rate_gain) < 10){stop("Not enough gain parameters")}
  if(length(rate_loss) < 10){stop("Not enough loss parameters")}
  if(length(fitness_costs) < 10){stop("Not enough fitness parameters")}
  
  #parameters for growth
  #nb that's just a logistic deterministic function in the model at the moment
  Nmax = 0.5 * 1e7
  # Death rate: timestep currently 1 hr, replicate in 20min in optimal... this captures immune system killing? 10% die? 
  death_rate <- 0.1
  
  # fitness cost - firstly without strain issues, then with human strain cost
  fitness_cost_all_noh <- rowSums(t(t(bacteria[,1:(ncol(bacteria)-2)]) * fitness_costs))
  fitness_cost_all <- fitness_cost_all_noh * c(rep(1, nrow(bacteria)/2),rep(rel_fit_human, nrow(bacteria)/2))
  
  #summary matrix to store number of bacteria in each strain at each timepoint
  all_results = matrix(0, nrow = tsteps, ncol = nrow(bacteria))
  all_results[1,] = bacteria[,"freq"]
  
  #summary matrix to store MGE prevalence at each timepoint in each parent
  all_mge_prev = matrix(0, nrow = 2*tsteps, ncol = (ncol(bacteria)-2))
  
  # Final check on input parameters 
  if(any(fitness_cost_all > 20) || growth_rate > 2.94){print(c("Input error (fitness or growth)",fitness_cost_all));  prev_predict = c(); totl_predict = c(); }else{
    
    for(t in 2:tsteps){
      
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
      MGE_prevalence_1 = Rfast::colsums(MGE_prevalence[1:nrow(bacteria)/2,])/sum(bacteria[1:nrow(bacteria)/2,"freq"]) # just in parent 1
      MGE_prevalence_2 = Rfast::colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) # just in parent 2
      
      #Store MGE prev
      all_mge_prev[(t-1),] = MGE_prevalence_1
      all_mge_prev[(1+tsteps) + (t-1),] = MGE_prevalence_2
      
      #for each strain present
      present_strains <- sample(which(bacteria[,"freq"]>0)) # do sample to give random order otherwise prioritise parent 1 growth at stationary phase
      
      for(i in present_strains){
        
        # In order for new profiles to appear at the start of the timestep remove
        # some of the existing strains - death - which will be replaced if lots of these remain
        # Background constant: relative to growth rate which is fitted 
        new_bacteria[i,"freq"] <- (1 - death_rate) * new_bacteria[i,"freq"] 
        
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
        probas[Rfast::rowsums(differences[,1:10])==0] = max(0, 1 - probas)
        #print(c("probas",probas,"gain",gain_probas))
        
        #use multinomial sampling to decide what the bacteria from strain i now become,
        # then add that amount to the updated matrix of bacteria numbers
        #here's where the "id" column in "differences" is useful: to align the indexing
        # between "differences" (which only contains valid transitions for strain i) and
        # "new_bacteria" (which contains all 2048 possible strains)
        #probas[is.na(probas)] <- 0 # got na errors in rmultinom - fix with this for now
        if(sum(probas) == "NA"){print(theta_in); break}
        
        new_bacteria[differences[,"id"],"freq"] = new_bacteria[differences[,"id"],"freq"] +
          pmin((Nmax - sum(new_bacteria[,"freq"]))/length(new_bacteria[i,]), rmultinom(1, bacteria[i, "freq"], probas)) # due to discrete time step, don't want to end up with negative bugs as more than Nmax at some point
        
        #if some bacteria remain in their original strain i, they now grow
        #currently just a deterministic logistic calculation
        # Need to add in fitness cost of elements: assume additive atm 
        new_bacteria[i,"freq"] = new_bacteria[i,"freq"] +
          min((Nmax - sum(new_bacteria[,"freq"])), round(bacteria[i,"freq"] * (1 - fitness_cost_all[i]) * growth_rate * (1 - sum(bacteria[,"freq"])/Nmax)))
        
        #print(c(sum(bacteria[,"freq"]),round(bacteria[i,"freq"] * (1 - fitness_cost_all[i]) * growth_rate * (1 - sum(bacteria[,"freq"])/Nmax))))
      }
      
      #update main bacteria matrix
      bacteria = new_bacteria
      #store numbers for each strain at that timepoint
      all_results[t,] = bacteria[,"freq"]
      
    }
    all_results_out <- all_results # store to check later
    
    # Last MGE prevalence
    MGE_prevalence = bacteria[,c(1:10)]*bacteria[,"freq"]
    MGE_prevalence_all = Rfast::colsums(MGE_prevalence)/tot_bacteria
    MGE_prevalence_1 = Rfast::colsums(MGE_prevalence[1:nrow(bacteria)/2,])/sum(bacteria[1:nrow(bacteria)/2,"freq"]) # just in parent 1
    MGE_prevalence_2 = Rfast::colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) # just in parent 2
    
    #Store MGE prev
    all_mge_prev[(tsteps),] = MGE_prevalence_1
    all_mge_prev[(2*tsteps),] = MGE_prevalence_2
    
    #clean up 
    all_results = as.data.frame(all_results_out)
    colnames(all_results) <- seq(1:2048)
    all_results$time = c(1:nrow(all_results))
    all_results = reshape2::melt(all_results, id.vars = "time")
    all_results$parent = c(rep(1, tsteps * nrow(bacteria)/2), rep(2, tsteps * nrow(bacteria)/2))
    
    # Totals - want prevalence of each parent strains
    totl_predict<- all_results %>% select(time,value, parent) %>%
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time, parent) %>% 
      summarise(total = sum(value),.groups = "drop") 
    
    # Prevalence of MGE at same time as data
    #clean up MGE
    all_mge_prev = as.data.frame(all_mge_prev)
    all_mge_prev$time = seq(1,tsteps,1)
    all_mge_prev$parent = c(rep(1, tsteps), rep(2, tsteps))
    all_mge_prev = reshape2::melt(all_mge_prev, id.vars = c("parent","time"))
    
    prev_predict <- all_mge_prev %>% 
      filter(time %in% c(4,48,96,288,384))
  }
  
  # Output
  
  profile_end1.1 <- c(1,1,1,1,1,1,1,1,0,0,1) # pig - same as at start
  profile_end2.1 <- c(0,1,0,0,1,0,0,0,1,0,2) # human - gains phi6 / phi2 / loses p4 (2/5/10)
  profile_end2.2 <- c(0,1,0,0,1,0,1,0,1,0,2) # human
  profiles_needed_end <- c(prodlim::row.match(profile_end1.1, as.data.frame(bacteria)%>%select(-freq), nomatch = NA),
                           prodlim::row.match(profile_end2.1, as.data.frame(bacteria)%>%select(-freq), nomatch = NA),
                           prodlim::row.match(profile_end2.2, as.data.frame(bacteria)%>%select(-freq), nomatch = NA))
  
  ### Likelihood
  if(!is.null(prev_predict)){
    
    #### element prevalence from model (a)
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    model_outputp <- prev_predict %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
    model_outputp$prev <- round(model_outputp$value,2)
    
    # Check what distribution of n_colonies at this prevalence in the model pig
    distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, variable, prob_all, n_colonies_prev)
    # e.g. to check 
    #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
    
    # lookup the probability from this distribution for the data
    likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","variable","n_colonies_prev")) %>% summarise(sum(log(prob_all)))
    
    #### total bugs output from model (b)
    model_outputt <- totl_predict
    
    likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
      as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))
    
    #### ensure certain profiles present (c)
    total_end <- unlist(all_results %>% filter(time == tsteps) %>% group_by(parent) %>% summarise(total = sum(value)) %>%select(total))
    ## Need to be present at > 80% and > 20% for parent 1 and parent 2 respectively 
    profile_end <- all_results %>% filter(time == tsteps, variable %in% profiles_needed_end) 
    if(nrow(profile_end) == 3){
      likelihood_profile_end <- profile_end %>% 
        mutate(prop = value / c(total_end[1],total_end[2],total_end[2]),
               cutoff = pmax(0,prop - 0.05),#c(0.8,0.2,0.2)), # more the better # not using -- too strict
               likelihood = log(prop)) %>% summarise(sum(likelihood))} else{likelihood_profile_end <- -Inf}
    
    #### Compare to data 
    compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals + 10 * likelihood_profile_end # add in a 10* weight for profile_end as otherwise only a small contribution relatively
  }else{compare_dat <- -Inf}
  
  # return log likelihood
  as.numeric(compare_dat)
}

### Labelling function
getLabels <- function(df) {
  match( do.call("paste", c(df[, , drop = FALSE],
                            sep = "\\r")),
         do.call("paste", c(unique(df)[, , drop
                                       = FALSE], sep = "\\r")) )
}

# Plot output

plot_circles <- function(profiles, output_data,plot_name){
  # profiles = ini$bacteria
  # output_data = time series from out (out$all_results)
  # where to store plot
  
  # Build data 
  bugs <- as.data.frame(profiles)
  bugs$variable = seq(1:nrow(bugs))
  bugs$parent = c(rep(1, nrow(bugs)/2), rep(2, nrow(bugs)/2))
  
  results <- output_data %>% filter(time %in% c(0,4,48,72,288,384), value > 0) %>% select(variable, value, parent, time)
  results$variable <- as.numeric(results$variable)
  tot <- results %>% group_by(parent, time) %>% summarise(tot = sum(value)) 
  results <- left_join(results, tot, by = c("parent", "time")) %>% mutate(perc = 100 * value / tot)
  piggy_data <- left_join(results, bugs %>% select(-freq), by = c("variable", "parent")) %>% arrange(value)
  piggy_data$pig <- 1
  piggy_data <- rename(piggy_data, c("profile" = "variable","freq" = "value",
                                     "v1" = "Var1","v2" = "Var2","v3" = "Var3","v4" = "Var4","v5" = "Var5","v6" = "Var6","v7" = "Var7","v8" = "Var8","v9" = "Var9","v10" = "Var10"))
  piggy_data <- as.data.frame(piggy_data)
  
  ######################## plot
  minir = 0.25
  bigr = 1
  x_centrebig = 5
  y_centrebig = 5
  x_centre_plasmid = c()
  y_centre_plasmid = 3
  
  #x(t) = r cos(t) + j
  #y(t) = r sin(t) + k
  pt = c(c(60,30,0,330,270)*pi/180, # on the edge
         3.5, 4.5, 5.5, # plasmid
         c(240*pi/180), 6.5) # Q
  
  pigg_plot <- piggy_data %>% pivot_longer(cols = v1:v10, names_to = "name") %>% 
    mutate(label = ifelse(value == 1, name,0)) %>% 
    mutate(mge = as.numeric(substr(name,2,3)),
           plasmid = ifelse(name %in% c("v1","v2","v3","v4","v5","v9"),0,1),
           x_centre = ifelse(plasmid == 0, bigr*cos(pt[mge]) + x_centrebig,pt[mge]),
           y_centre = ifelse(plasmid == 0, bigr*sin(pt[mge]) + x_centrebig,y_centre_plasmid))
  
  pigg_plotp <- pigg_plot %>% filter(parent == 1)
  pigg_plotq <- pigg_plot %>% filter(parent == 2)
  
  mm <- max(as.numeric(pigg_plotp$profile))
  pigg_plotp$profile = factor(pigg_plotp$profile, levels=seq(mm,1,-1)) # reorder profiles to match McCarthy Fig2 JUST FOR P
  
  pigg_plotp$perc_lab <- paste0(signif(pigg_plotp$perc,2),"%")
  pigg_plotq$perc_lab <- paste0(signif(pigg_plotq$perc,2),"%")
  
  gp <- ggplot(pigg_plotp %>% filter(perc > 5), aes(x0=x_centre, y0 = y_centre, group = profile)) + geom_circle(aes(r = minir, col= label, fill  = label)) + 
    geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
    geom_text(aes(x = 5, y = 6.5, label = perc_lab)) +
    facet_grid(time ~ profile + parent) + 
    scale_fill_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"), 
                      values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8"),drop = FALSE)+ 
    scale_color_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"), 
                       values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8"),drop = FALSE) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  gq <- ggplot(pigg_plotq %>% filter(perc > 5), aes(x0=x_centre, y0 = y_centre, group = profile)) + geom_circle(aes(r = minir, col= label, fill  = label)) + 
    geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
    geom_text(aes(x = 5, y = 6.5, label = perc_lab)) +
    facet_grid(time ~ profile + parent) + 
    scale_fill_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"), 
                      values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8"),drop = FALSE)+ 
    scale_color_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"), 
                       values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8"),drop = FALSE) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  ## Only need to do first time
  g <- gp +  gq + plot_layout(guides = 'collect') +
    plot_layout(widths = c(1, 1))
  ggsave(paste0("plots/",plot_name,".pdf"), width = 45, height = 25)
}

plot_time_series <- function(out, plot_name){
  # out = all data from run 
  # plot_name = where to plot
  theme_set(theme_bw(base_size = 11))
  preva = out$prev_predict
  tots = out$totl_predict
  ALL = out$all_results
  
  if(!is.null(preva)){
    
    #### element prevalence from model (a)
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    model_outputp <- preva %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
    model_outputp$prev <- round(model_outputp$value,2)
    
    # Check what distribution of n_colonies at this prevalence in the model pig
    distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, variable, prob_all, n_colonies_prev)
    # e.g. to check 
    #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
    
    # lookup the probability from this distribution for the data
    likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","variable","n_colonies_prev")) %>% summarise(sum(log(prob_all)))
    
    #### total bugs output from model (b)
    model_outputt <- tots
    
    likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
      as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))
    
    #### ensure certain profiles present (c)
    total_end <- unlist(ALL %>% filter(time == tsteps) %>% group_by(parent) %>% summarise(total = sum(value)) %>%select(total))
    ## Need to be present at > 80% and > 20% for parent 1 and parent 2 respectively 
    profile_end <- ALL %>% filter(time == tsteps, variable %in% profiles_needed_end) 
    if(nrow(profile_end) == 3){
      likelihood_profile_end <- profile_end %>% 
        mutate(prop = value / c(total_end[1],total_end[2],total_end[2]),
               cutoff = pmax(0,prop - 0.01),#c(0.8,0.2,0.2)), # more the better
               likelihood = log(prop)) %>% summarise(sum(likelihood))} else{likelihood_profile_end <- -Inf}
    
    #### Compare to data 
    compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals + 10 * likelihood_profile_end # add in a 10* weight for profile_end as otherwise only a small contribution relatively
  }else{compare_dat <- -Inf}
  
  model_outputp <- rename(model_outputp, parent_strain = parent)
  model_outputp$name <- recode(model_outputp$variable, V2 = "phi6",V5 = "phi2",V6 = "p1",V7 = "p2",V8 = "p3",V10 ="p4")
  
  g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")),
               aes(x=time, y = sum_prop, group = interaction(name, pig))) +
    geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) +
    geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) +
    geom_point(aes(col = name),size = 1.5) +
    facet_wrap(name~parent_strain, ncol = 2) + 
    ggtitle(paste0("likelihood = ", round(compare_dat,3), "(",round(likelihood_lookup_elements,3),"+",
                                                                    round(likelihood_lookup_totals,3),"+",round(likelihood_profile_end,3),")"))
  
  g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
    geom_point(aes(col = factor(parent))) +
    geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") +
    scale_y_log10() +
    geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))
  
  g3 <- ggplot(ALL %>% filter(variable %in% profiles_needed_end), aes(x=time, y = value)) + geom_line(aes(col = variable)) + geom_vline(xintercept = tsteps) + 
    geom_hline(yintercept = c(0.05) * max(ALL$value)) + scale_color_manual(values = c("red","green","blue"), breaks = c(profiles_needed_end))
  
  ## Only need to do first time
  g <- g1 / (g2 + g3)
  ggsave(paste0("plots/",plot_name,".pdf"))
  
  
}
