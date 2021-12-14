#### code from Quentin - using transition matrix 
library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(Rfast) #for faster rowsums/colsums

start_time = Sys.time()

#matrix to store all possible MGE combinations (ie all possible strains)
#MGE1 presence / MGE2 presence / MGE3 presence / ... / freq / parent
bacteria = rbind(as.matrix(expand.grid(rep(list(0:1), 10))),
                 as.matrix(expand.grid(rep(list(0:1), 10))))
#just randomly deciding which 2 strains start at 1000 bacteria at the moment
bacteria = cbind(bacteria, sample(c(1000, 1000, rep(0, nrow(bacteria)-2)), nrow(bacteria), replace = F))
bacteria = cbind(bacteria, c(rep(1, nrow(bacteria)/2), rep(2, nrow(bacteria)/2)))

#set column names for easier reference later
colnames(bacteria)[c(ncol(bacteria)-1, ncol(bacteria))] = c("freq", "parent")

#parameters for growth
#nb that's just a logistic deterministic function in the model at the moment
growth_rate = 0.1
Nmax = 1e6

#rate of each MGE loss
#pretty random values currently
rate_loss= c(0.01,
             0.01,
             0.01,
             0.01,
             0.01,
             0.05,
             0.01,
             0.01,
             0.05,
             0.01)

#rate of MGE gain
rate_gain = c(0.05,
              0.01,
              0.01,
              0.01,
              0.01,
              0.01,
              0.01,
              0.01,
              0.01,
              0.05)

#times to repeat simulation
times = 384

#summary matrix to store number of bacteria in each strain at each timepoint
all_results = matrix(0, nrow = times, ncol = nrow(bacteria))
all_results[1,] = bacteria[,"freq"]


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
  MGE_prevalence = colsums(MGE_prevalence)/tot_bacteria
  
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
    gain_probas = 1-exp(-rate_gain*bacteria[i,"freq"]*MGE_prevalence)
    gain_probas = pmin(gain_probas, 1-bacteria_profile)
    
    #sum probas (since each MGE can either be gained or lost only)
    #now, these probas represent the probability for each MGE status to change during this timestep
    #ie an MGE already present can only be lost, and an MGE currently absent can only be gained
    probas = lose_probas + gain_probas
    
    #recover the difference matrix for our strain i
    differences = difference_list[[i]]
    
    #this aligns the probabilities with the other strains they correspond to
    #eg if there's a 0.5 proba to lose MGE 1, and that losing MGE 1 will change
    # the profile of strain i to the profile of strain 10, this will translate 
    # to a 0.5 proba of bacteria of strain i becoming strain 10 during this timestep
    #inversing the matrices twice is the fastest option to do this (I think?)
    probas = abs(t(t(differences[,-ncol(differences)]) * probas))
    
    #collapse to vector using rowsums (only 1 value per row will be greater than 0)
    probas = rowsums(probas)
    #okay, this is needed to set a probability of the strain NOT changing profile
    #it's not ideal, because sometimes the probabilities add up to more than 1, hence the max() function
    #something to look into...
    probas[probas == 0] = max(0, 1 - probas)
    
    #use multinomial sampling to decide what the bacteria from strain i now become,
    # then add that amount to the updated matrix of bacteria numbers
    #here's where the "id" column in "differences" is useful: to align the indexing
    # between "differences" (which only contains valid transitions for strain i) and
    # "new_bacteria" (which contains all 2048 possible strains)
    new_bacteria[differences[,"id"],"freq"] = new_bacteria[differences[,"id"],"freq"] +
      rmultinom(1, bacteria[i, "freq"], probas)
    
    #if some bacteria remain in their original strain i, they now grow
    #currently just a deterministic logistic calculation
    new_bacteria[i,"freq"] = new_bacteria[i,"freq"] +
      round(bacteria[i,"freq"] * growth_rate * (1 - sum(bacteria[,"freq"])/Nmax))
    
  }
  
  #update main bacteria matrix
  bacteria = new_bacteria
  #store numbers for each strain at that timepoint
  all_results[t,] = bacteria[,"freq"]
  
}


#clean up and plot
all_results = as.data.frame(all_results)
all_results$time = c(1:nrow(all_results))

all_results = melt(all_results, id.vars = "time")

all_results %>%
  filter(value > 0) %>%
  ggplot() +
  geom_line(aes(x = time, y = value, group = variable)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  theme_bw()

Sys.time() - start_time
