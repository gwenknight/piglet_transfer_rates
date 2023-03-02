########  Transition matrix function 

gain_rate <- 0.1
loss_rate <- 0.1
growth_rate <- 0.1

library(combinat)
# get the list of all permutations
my_list<-c(combinat::permn(c(1,0)),combinat::permn(c(1,1)),combinat::permn(c(0,0)))
#my_list<-c(combinat::permn(c(1,0,0)),combinat::permn(c(1,1,1)),combinat::permn(c(1,1,0)))
# my_list<-c(combinat::permn(c(1,0,0)),combinat::permn(c(1,1,1)),combinat::permn(c(1,1,0)))

# convert the list to a matrix
profiles<-do.call(rbind,my_list)
#take the unique rows
profiles<-unique(profiles)
profiles # ALL combinations of 1s and 0s

# Go through each profile? 
for(i in 1:nrow(profiles)){
  w_gains<-which(profiles[i,]==0)
  w_loss<-which(profiles[i,]==1)
  
  # Gain elements
  diff <- profiles - profiles[i,]
  w_gain_from <- which(diff == 1,arr.ind = TRUE)
  
  w_loss_to <- which(profiles[,w_loss] == 0,arr.ind = TRUE)
  
  
}

