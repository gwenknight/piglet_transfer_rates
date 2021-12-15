#### code from Quentin - using transition matrix 

### Function to move MGE

piglet_mrsa_movement <- function(times, parameters){
  # tsteps = number of time steps
  # theta = parameters
  
  # Assign parameter values
  # If have rates for all 10 elements
  if(length(parameters) == 31){
    rate_loss <- as.numeric(parameters[1:10])
    rate_gain <- as.numeric(parameters[11:20] )
    fitness_costs <- as.numeric(parameters[21:30])
    growth_rate <- as.numeric(parameters[31])
  }
  
  # If just looking at the elements that move
  if(length(parameters) < 31 && length(parameters) > 7){
    rate_loss <- as.numeric(c(0,parameters["mu2"],0,0,parameters["mu5"],parameters["mu6"],parameters["mu7"],parameters["mu8"],0,parameters["mu10"]))
    rate_gain <- as.numeric(c(0,parameters["gamma2"],0,0,parameters["gamma5"],parameters["gamma6"],parameters["gamma7"],parameters["gamma8"],0,parameters["gamma10"]))
    fitness_costs <- as.numeric(c(0,parameters["f2"],0,0,parameters["f5"],parameters["f6"],parameters["f7"],parameters["f8"],0,parameters["f10"]))
    growth_rate <- as.numeric(parameters["grow"])
  }
  
  # If fixed input - same rates for all elements
  if(length(parameters) == 4){
    if(is.null(names(parameters))){stop("No names for parameters")}
    rate_loss <- as.numeric(c(0,parameters["mu"],0,0,parameters["mu"],parameters["mu"],parameters["mu"],parameters["mu"],0,parameters["mu"]))
    rate_gain <- as.numeric(c(0,parameters["gamma"],0,0,parameters["gamma"],parameters["gamma"],parameters["gamma"],parameters["gamma"],0,parameters["gamma"]))
    fitness_costs <- as.numeric(c(0,parameters["f"],0,0,parameters["f"],parameters["f"],parameters["f"],parameters["f"],0,parameters["f"]))
    growth_rate <- as.numeric(parameters["grow"])
  }
  
  # If fixed input - same rates for phage vs plasmids
  if(length(parameters) == 7){
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    rate_loss <- as.numeric(c(0,parameters["mu_phage"],0,0,parameters["mu_phage"],parameters["mu_plasmid"],parameters["mu_plasmid"],parameters["mu_plasmid"],0,parameters["mu_plasmid"]))
    rate_gain <- as.numeric(c(0,parameters["gamma_phage"],0,0,parameters["gamma_phage"],parameters["gamma_plasmid"],parameters["gamma_plasmid"],parameters["gamma_plasmid"],0,parameters["gamma_plasmid"]))
    fitness_costs <- as.numeric(c(0,parameters["f_phage"],0,0,parameters["f_phage"],parameters["f_plasmid"],parameters["f_plasmid"],parameters["f_plasmid"],0,parameters["f_plasmid"]))
    growth_rate <- as.numeric(parameters["grow"])
  }
  
  # Errors
  if(length(rate_gain) < 10){stop("Not enough gain parameters")}
  if(length(rate_loss) < 10){stop("Not enough loss parameters")}
  if(length(fitness_costs) < 10){stop("Not enough fitness parameters")}
  if(sum(growth_in) > 1){stop("Fitness cost too large")}
  
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
  
  #parameters for growth
  #nb that's just a logistic deterministic function in the model at the moment
  Nmax = 1e6
  # fitness cost
  fitness_cost_all <- rowSums(t(t(bacteria[,1:(ncol(bacteria)-2)]) * fitness_costs))
  
  #summary matrix to store number of bacteria in each strain at each timepoint
  all_results = matrix(0, nrow = times, ncol = nrow(bacteria))
  all_results[1,] = bacteria[,"freq"]
  
  #summary matrix to store MGE prevalence at each timepoint in each parent
  all_mge_prev = matrix(0, nrow = 2*times, ncol = (ncol(bacteria)-2))
  
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
  
  
  for(t in 2:times){
    
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
    MGE_prevalence_1 = colsums(MGE_prevalence[1:nrow(bacteria)/2,])/sum(bacteria[1:nrow(bacteria)/2,"freq"]) # just in parent 1
    MGE_prevalence_2 = colsums(MGE_prevalence[(1+nrow(bacteria)/2):nrow(bacteria),])/sum(bacteria[(1+nrow(bacteria)/2):nrow(bacteria),"freq"]) # just in parent 2
    
    #Store MGE prev
    all_mge_prev[(t-1),] = MGE_prevalence_1
    all_mge_prev[(1+times) + (t-1),] = MGE_prevalence_2
    
    #for each strain:
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
  all_mge_prev[(times),] = MGE_prevalence_1
  all_mge_prev[(2*times),] = MGE_prevalence_2
  
  #clean up 
  all_results = as.data.frame(all_results_out)
  all_results$time = c(1:nrow(all_results))
  all_results = melt(all_results, id.vars = "time")
  all_results$parent = c(rep(1, times * nrow(bacteria)/2), rep(2, times * nrow(bacteria)/2))
  
  # Totals - want prevalence of each parent strains
  totl_predict<- all_results %>% select(time,value, parent) %>%
    filter(time %in% c(4,48,96,288,384)) %>% group_by(time, parent) %>% 
    summarise(total = sum(value),.groups = "drop") 
  
  # Prevalence of MGE at same time as data
  #clean up MGE
  all_mge_prev = as.data.frame(all_mge_prev)
  all_mge_prev$time = seq(1,times,1)
  all_mge_prev$parent = c(rep(1, times), rep(2, times))
  all_mge_prev = melt(all_mge_prev, id.vars = c("parent","time"))
  
  prev_predict <- all_mge_prev %>% 
    filter(time %in% c(4,48,96,288,384))
  
  # Output
  return(list(all_results = all_results, 
              prev_predict = prev_predict, totl_predict = totl_predict))
}


#piglet_mrsa_movement(384, theta)



### Calculate log posterior from model at theta parameter values
run_sim_logPosterior <- function(theta){
  tsteps = 384
  
  ## Run for these parameters
  out <- piglet_mrsa_movement(tsteps, theta)
  
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
  
  # return log likelihood
  as.numeric(compare_dat)
}
