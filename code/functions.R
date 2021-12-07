##### Gain function 
library(mgcv) # for uniquecombs function

# Functions to work out what MGE move


################################################################################################################################################################################################################################
########################################################************************ GAIN *******************************################################################################################################################
################################################################################################################################################################################################################################

# ORIGINAL
gain <- function(receiver, giver, mu_in){
  # receiver: bugs that are getting MGE
  # giver: bugs that are giving MGE
  # mu_in: rate of movement for each of the 10 MGE in this system
  
  ## Run through each row in receiver
  # Parallel
  #p11 <- foreach (j=1:dim(receiver)[1], .combine=rbind) %dopar% { # for each profile in receiver
  
  # Total number of bugs
  nbugs = sum(receiver$freq) + sum(giver$freq)  # total bugs
  
  # not parallel
  p11 <- c()
  for(j in 1:dim(receiver)[1]){
    
    pk <- c()
    for(k in 1:dim(giver)[1]){ # for each "other profile" in giver
      
      nh = round(receiver[j,"freq"]*giver[k,"freq"] / nbugs,0) # mass action: probability bump into this bug
      rph <- matrix(0,nh,10)
      
      ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
      if(nh * max(mu_in) < 5){
        # for each element - does transfer happen between this and the profiles? 
        for(ii in 1:10){ 
          if(giver[k,ii] - receiver[j,ii] > 0){ # so transfer can happen from profile k to profile j (doesn't matter if j = k, receiver = giver as this will be 0)
            rph[,ii] <- rbinom(n = nh, size = 1, mu_in[ii]) # for each bacteria, probability that transfer happens
          }
        }
      }else{
        # for each element - does transfer happen between this and the profiles? 
        for(ii in 1:10){ 
          if(giver[k,ii] - receiver[j,ii] > 0){ # so transfer can happen from profile k to profile j (doesn't matter if j = k, receiver = giver as this will be 0)
            nmovements = round(rnorm(n = 1, mean = nh*mu_in[ii],sd = sqrt(nh * mu_in[ii] * (1 - mu_in[ii]))),0)
            if(nmovements > 0){
              rph[runif(nmovements, 1, nh), ii] <- 1 # assign randomly to each of the nh bugs
            }
          }
        }
      }
      
      ## rph has all the profiles / moved elements
      if(sum(rph)>0){ # only do if there were some transfers
        rph <- as.data.frame(rph)
        colnames(rph) = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10") # 10 elements
        
        # Resulting new profiles and their count
        p1new <- plyr::count(rph, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")) %>% 
          mutate(Total = select(., v1:v10) %>% rowSums()) %>%  # count how many
          filter(Total>0) %>% select(-Total) # remove no transfer row
        
        # Have to add to existing profile
        if(dim(p1new)[1]>0){
          for(i in 1:dim(p1new)[1]){
            p1new[i,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")] <- 
              receiver[j,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")] + 
              p1new[i,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")]
          }
          p1new$j <- j # label which profile became this new one
          p1new$k <- k 
          p1new$ii <- ii 
          pk <- rbind(pk, p1new)
        }}
    }
    
    #pk # return pk from the parallel for each profile in receiver
    p11 <- rbind(p11, pk)
  }
  
  return(p11)
}


## Don't care who get from 
gain1 <- function(receiver, giver, mu_in){
  # receiver: bugs that are getting MGE
  # giver: bugs that are giving MGE
  # mu_in: rate of movement for each of the 10 MGE in this system
  
  ## Run through each row in receiver
  # Parallel
  #p11 <- foreach (j=1:dim(receiver)[1], .combine=rbind) %dopar% { # for each profile in receiver
  
  # Total number of bugs
  nbugs = sum(receiver$freq) + sum(giver$freq)  # total bugs
  
  # All combination 
  all_comb <- as.data.frame(expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1)))
  colnames(all_comb) <- mge
  
  # Giver - how many bugs have each MGE? 
  giver_mge <- colSums(giver[,mge]*giver$freq)
  pos_mge_giver <- as.numeric(giver_mge > 0)
  
  # Receiver - how many bugs have each MGE? 
  recvr_mge <- colSums(receiver[,mge]*receiver$freq)
  
  # not parallel
  p11 <- c()
  for(j in 1:dim(receiver)[1]){
    
    #left_join(all_comb, receiver) %>% filter(freq > 0)
    
    # What are the possible next profiles for receiver ones to gain mge to go to? 
    # need a 1 in all differences between all_comb and receiver
    #all_comb <- all_comb[1020:1024,] # for checks
    poss_next_gains <- all_comb[unique(which(sweep(all_comb, 2, as.numeric(receiver[j,mge])) == 1, arr.ind = TRUE)[,"row"]),] # Which profiles have a gain
    poss_next <- sweep(poss_next_gains,2, as.numeric(receiver[j,mge]))
    how_get_to_next <- poss_next[rowSums(poss_next < 0) == 0,] # Only those that are a gain, row numbers those of all_comb
    keep_profile <- c()
    for(ip in 1:dim(how_get_to_next)[1]){ifelse(all((pos_mge_giver - how_get_to_next[ip,]) > -0.1),keep_profile <- c(keep_profile,ip),0)}
    how_get_to_next <- how_get_to_next[keep_profile,]
    for(im in 1:10){how_get_to_next[,im] <- how_get_to_next[,im] * giver_mge[im]}
    
    
    pk <- c()
    for(k in 1:dim(giver)[1]){ # for each "other profile" in giver
      
      
      nh = round(receiver[j,"freq"]*giver[k,"freq"] / nbugs,0) # mass action: probability bump into this bug
      rph <- matrix(0,nh,10)
      
      ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
      if(nh * max(mu_in) < 5){
        # for each element - does transfer happen between this and the profiles? 
        for(ii in 1:10){ 
          if(giver[k,ii] - receiver[j,ii] > 0){ # so transfer can happen from profile k to profile j (doesn't matter if j = k, receiver = giver as this will be 0)
            rph[,ii] <- rbinom(n = nh, size = 1, mu_in[ii]) # for each bacteria, probability that transfer happens
          }
        }
      }else{
        # for each element - does transfer happen between this and the profiles? 
        for(ii in 1:10){ 
          if(giver[k,ii] - receiver[j,ii] > 0){ # so transfer can happen from profile k to profile j (doesn't matter if j = k, receiver = giver as this will be 0)
            nmovements = round(rnorm(n = 1, mean = nh*mu_in[ii],sd = sqrt(nh * mu_in[ii] * (1 - mu_in[ii]))),0)
            if(nmovements > 0){
              rph[runif(nmovements, 1, nh), ii] <- 1 # assign randomly to each of the nh bugs
            }
          }
        }
      }
      
      ## rph has all the profiles / moved elements
      if(sum(rph)>0){ # only do if there were some transfers
        rph <- as.data.frame(rph)
        colnames(rph) = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10") # 10 elements
        
        # Resulting new profiles and their count
        p1new <- plyr::count(rph, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")) %>% 
          mutate(Total = select(., v1:v10) %>% rowSums()) %>%  # count how many
          filter(Total>0) %>% select(-Total) # remove no transfer row
        
        # Have to add to existing profile
        if(dim(p1new)[1]>0){
          for(i in 1:dim(p1new)[1]){
            p1new[i,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")] <- receiver[j,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")] + p1new[i,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")]
          }
          p1new$j <- j # label which profile became this new one
          p1new$k <- k 
          p1new$ii <- ii 
          pk <- rbind(pk, p1new)
        }}
    }
    
    #pk # return pk from the parallel for each profile in receiver
    p11 <- rbind(p11, pk)
  }
  
  return(p11)
}

#### MGE based
gain_mge <- function(prnt1, prnt2, mu_in){
  # prnt/2: parent 1 or 2 profiles
  # mu_in: rate of movement for each of the 10 MGE in this system
  
  d1 <- dim(prnt1)[1]
  
  # total 
  prnt <- rbind(prnt1, prnt2)
  
  nprnt1 <- prnt1; nprnt2 <- prnt2;
  
  # Total number of bugs with each MGE
  nbugs_with_orig = colSums(prnt$freq * prnt[,1:10]) #sum(prnt$freq * prnt[,m])
  total_bugs = sum(prnt$freq)
  
  # random order
  order_mge <- sample(as.numeric(which(nbugs_with_orig > 0))) # ensures only do mge with > 0 bugs
  
  # As update total 
  for(m in order_mge){
    
    # Total number of bugs with each MGE
    nbugs_with = nbugs_with_orig[m]
    nbugs_without =  total_bugs - nbugs_with
    
    # Number transfer: mass action assumption
    nh <- nbugs_with * (nbugs_without / total_bugs)
    ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
    if(nh * as.numeric(mu_in[m]) < 5){
      n_transfer <- sum(rbinom(n = nh, size = 1, as.numeric(mu_in[m]))) # how many transfer = sum of 1s = successful binomial transfer events
    }else{
      n_transfer <- round(rnorm(n = 1, mean = nh*as.numeric(mu_in[m]),sd = sqrt(nh * as.numeric(mu_in[m]) * (1 - as.numeric(mu_in[m])))),0)
    }
    
    # If any transfers
    if(n_transfer > 0){
      # Blank to store these transfers
      new_prnt1 <- c(); new_prnt2 <- c();
      # Which nbugs_without profile get the MGE? 
      which_get = round(runif(n_transfer, 1,nbugs_without),0)
      profile_without_ind = intersect(which(prnt[,m] == 0), which(prnt[,"freq"]>0))
      allocated_bins <- table(cut(which_get, breaks = c(0, cumsum(prnt[profile_without_ind,"freq"]))))
      allocated <- as.numeric(allocated_bins) # how many to each profile?
      
      # Build new profiles
      alloc_prof <- which(allocated>0)
      for(i in 1:length(alloc_prof)){
        # new profile
        profile <- prnt[profile_without_ind[i],]
        profile[m] <- 1 # has a one at this mge
        profile["freq"] <- allocated[alloc_prof[i]]
        if(profile_without_ind[i] <= d1){ # if less than d1 then in first parent
          # New profiles come from old ones
          nprnt1[profile_without_ind[i],"freq"] <- nprnt1[profile_without_ind[i],"freq"] - allocated[i]
          new_prnt1 <- rbind(new_prnt1, profile)
        }else{
          nprnt2[profile_without_ind[i]-d1,"freq"] <- nprnt2[profile_without_ind[i]-d1,"freq"] - allocated[i]
          new_prnt2 <- rbind(new_prnt2, profile)
        }
      }
      
      nprnt1 <- rbind(nprnt1, new_prnt1)
      nprnt2 <- rbind(nprnt2, new_prnt2)
      # remove any zero profiles
      nprnt1 <- nprnt1[which(nprnt1$freq>0),]
      nprnt2 <- nprnt2[which(nprnt2$freq>0),]
    }
    
    
  }
  
  return(list(Pnew = nprnt1, Qnew = nprnt2))
}

#### MGE based
## with multivariate sampling
gain_mge_all <- function(prnt1, prnt2, mu_in){
  # prnt/2: parent 1 or 2 profiles
  # mu_in: rate of movement for each of the 10 MGE in this system
  
  d1 <- dim(prnt1)[1]
  
  # total 
  prnt <- rbind(prnt1, prnt2)
  
  nprnt1 <- prnt1; nprnt2 <- prnt2;
  
  # Total number of bugs with each MGE
  nbugs_with_orig = colSums(prnt$freq * prnt[,1:10]) #sum(prnt$freq * prnt[,m])
  total_bugs = sum(prnt$freq)
  nbugs_without_orig <- total_bugs - nbugs_with_orig
  
  # Total number of bugs with each MGE
  nbugs_with = nbugs_with_orig[m]
  nbugs_without =  total_bugs - nbugs_with
  
  # Number transfer: mass action assumption
  nh <- nbugs_with_orig * nbugs_without_orig / total_bugs #nbugs_with * (nbugs_without / total_bugs)
  
  if(sum(nh)>0){ # if any transfers to do! 
    
    pos <- which(nbugs_with_orig > 0)
    pos_norm <- which(nbugs_with_orig > 5) # these are normally sampled
    l_norm <- length(pos_norm)
    pos_binom <- setdiff(pos, pos_norm) # there are binomally sampled
    l_binom <- length(pos_binom)
    
    n_transfer <- matrix(0,1,10)
    
    ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
    if(l_binom > 0){ # if any np / nq < 5
      size_b = round(as.numeric(nh[pos_binom]),0)
      n_gain_binom <- rbinom( n = l_binom, size = size_b, prob = c(as.numeric(mu_in[pos_binom])))
      n_transfer[pos_binom] <- n_gain_binom # store
    }
    
    if(l_norm > 0){ # if any np / nq > 5
      size_n = round(as.numeric(mu_in[pos_norm]) * as.numeric(nh[pos_norm]),0)
      p_n = sqrt(size_n * (1-as.numeric(mu_in[pos_norm])))
      if(l_norm > 1){
        n_gain_norm <- round(rmvnorm( n = 1, mean = size_n, sigma = diag(p_n) ),0)
      } else {
        n_gain_norm <- round(rmvnorm( n = 1, mean = size_n, sd = p_n ),0)
      }
      n_transfer[pos_norm] <- n_gain_norm # store
    }
    
    
    # If any transfers
    if(sum(n_transfer) > 0){
      
      # random order
      order_mge <- sample(as.numeric(which(n_transfer > 0))) # ensures only do mge with > 0 bugs
      
      for(m in order_mge){
        # Blank to store these transfers
        new_prnt1 <- c(); new_prnt2 <- c();
        # Which nbugs_without profile get the MGE? 
        which_get = round(runif(n_transfer[m], 1,nbugs_without_orig[m]),0)
        profile_without_ind = intersect(which(prnt[,m] == 0), which(prnt[,"freq"]>0))
        allocated_bins <- table(cut(which_get, breaks = c(0, cumsum(prnt[profile_without_ind,"freq"]))))
        allocated <- as.numeric(allocated_bins) # how many to each profile?
        
        # Build new profiles
        alloc_prof <- which(allocated>0)
        for(i in 1:length(alloc_prof)){
          # new profile
          profile <- prnt[profile_without_ind[i],]
          profile[m] <- 1 # has a one at this mge
          profile["freq"] <- allocated[alloc_prof[i]]
          if(profile_without_ind[i] <= d1){ # if less than d1 then in first parent
            # New profiles come from old ones
            nprnt1[profile_without_ind[i],"freq"] <- nprnt1[profile_without_ind[i],"freq"] - allocated[i]
            new_prnt1 <- rbind(new_prnt1, profile)
          }else{
            nprnt2[profile_without_ind[i]-d1,"freq"] <- nprnt2[profile_without_ind[i]-d1,"freq"] - allocated[i]
            new_prnt2 <- rbind(new_prnt2, profile)
          }
        }
        
        nprnt1 <- rbind(nprnt1, new_prnt1)
        nprnt2 <- rbind(nprnt2, new_prnt2)
        # remove any zero profiles
        nprnt1 <- nprnt1[which(nprnt1$freq>0),]
        nprnt2 <- nprnt2[which(nprnt2$freq>0),]
      }
    }
    
    
  }
  
  return(list(Pnew = nprnt1, Qnew = nprnt2))
}

#### MGE based
## with multivariate sampling
## and allocation structure
gain_mge_all_new <- function(prnt1, prnt2, mu_in){
  # prnt/2: parent 1 or 2 profiles
  # mu_in: rate of movement for each of the 10 MGE in this system
  
  d1 <- dim(prnt1)[1]
  d2 <- dim(prnt2)[1]
  
  # total 
  prnt <- rbind(prnt1, prnt2)
  
  # New matrices
  nprnt1 <-  matrix(0, nrow = 15000, ncol = dim(prnt1)[2])
  colnames(nprnt1) <- colnames(prnt1)
  nprnt2 <-  matrix(0, nrow = 15000, ncol = dim(prnt1)[2])
  colnames(nprnt2) <- colnames(prnt1)
  
  nprnt1 <- rbind(prnt1, nprnt1) # if allocate get list conversion
  nprnt2 <- rbind(prnt2, nprnt2)
  
  row_for_new_1 <- d1 + 1
  row_for_new_2 <- d2 + 1 
  
  # Blank matrix 
  blank <- matrix(0,1000,dim(prnt1)[2]); colnames(blank) <- colnames(prnt1); blank <- as.data.frame(blank)
  
  # Total number of bugs with each MGE
  nbugs_with_orig = colSums(prnt$freq * prnt[,1:10]) #sum(prnt$freq * prnt[,m])
  total_bugs = sum(prnt$freq)
  nbugs_without_orig <- total_bugs - nbugs_with_orig
  
  # Number transfer: mass action assumption
  nh <- nbugs_with_orig * nbugs_without_orig / total_bugs #nbugs_with * (nbugs_without / total_bugs)
  
  if(sum(nh)>0){ # if any transfers to do! 
    
    pos <- which(nbugs_with_orig > 0)
    pos_norm <- which(nbugs_with_orig > 5) # these are normally sampled
    l_norm <- length(pos_norm)
    pos_binom <- setdiff(pos, pos_norm) # there are binomally sampled
    l_binom <- length(pos_binom)
    
    n_transfer <- matrix(0,1,10)
    
    ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
    if(l_binom > 0){ # if any np / nq < 5
      size_b = round(as.numeric(nh[pos_binom]),0)
      n_gain_binom <- rbinom( n = l_binom, size = size_b, prob = c(as.numeric(mu_in[pos_binom])))
      n_transfer[pos_binom] <- n_gain_binom # store
    }
    
    if(l_norm > 0){ # if any np / nq > 5
      size_n = round(as.numeric(mu_in[pos_norm]) * as.numeric(nh[pos_norm]),0)
      p_n = sqrt(size_n * (1-as.numeric(mu_in[pos_norm])))
      if(l_norm > 1){
        n_gain_norm <- round(rmvnorm( n = 1, mean = size_n, sigma = diag(p_n) ),0)
      } else {
        n_gain_norm <- round(rnorm( n = 1, mean = size_n, sd = p_n ),0)
      }
      n_transfer[pos_norm] <- n_gain_norm # store
    }
    
    
    # If any transfers
    if(sum(n_transfer) > 0){
      
      # random order
      wo <- which(n_transfer>0)
      if(length(wo)>1){
        order_mge <- sample(as.numeric(wo)) # ensures only do mge with > 0 bugs
      } else {order_mge <- wo}
      
      for(m in order_mge){
        # Blank to store these transfers
        new_prnt1 <- blank; new_prnt2 <- blank;
        
        index1 <- 0; index2 <- 0
        
        # Which nbugs_without profile get the MGE? 
        which_get = round(runif(n_transfer[m], 1,nbugs_without_orig[m]),0)
        profile_without_ind = intersect(which(prnt[,m] == 0), which(prnt[,"freq"]>0))
        allocated_bins <- table(cut(which_get, breaks = c(0, cumsum(prnt[profile_without_ind,"freq"]))))
        allocated <- as.numeric(allocated_bins) # how many to each profile?
        
        # Build new profiles
        alloc_prof <- which(allocated>0)
        
        profiles <- prnt[profile_without_ind[alloc_prof],]
        profiles[m] <- 1
        profiles["freq"] <- allocated[alloc_prof]
        
        p1_prof <- which(profile_without_ind[alloc_prof] <= d1) # these are parent 1 profiles
        p2_prof <- which(profile_without_ind[alloc_prof] > d1) # these are parent 1 profiles
        
        if(length(p1_prof)>0){
          nprnt1[profile_without_ind[alloc_prof[p1_prof]],"freq"] <- nprnt1[profile_without_ind[alloc_prof[p1_prof]], "freq"] - allocated[alloc_prof[p1_prof]]
          prnt[profile_without_ind[alloc_prof[p1_prof]],"freq"] <- prnt[profile_without_ind[alloc_prof[p1_prof]], "freq"] - allocated[alloc_prof[p1_prof]]
          new_prnt1[(1+index1):(index1 + length(p1_prof)),] <- profiles[1:(1+length(p1_prof)-1),] 
        }
        if(length(p2_prof)>0){
          nprnt2[(profile_without_ind[alloc_prof[p2_prof]]-d1),"freq"] <- nprnt2[(profile_without_ind[alloc_prof[p2_prof]]-d1), "freq"] - allocated[alloc_prof[p2_prof]]
          prnt[(profile_without_ind[alloc_prof[p2_prof]]),"freq"] <- prnt[(profile_without_ind[alloc_prof[p2_prof]]), "freq"] - allocated[alloc_prof[p2_prof]] # limit on how many can add to
          new_prnt2[(1+index2):(index2 + length(p2_prof)),] <- profiles[(1 + length(p1_prof)):(1 + length(p1_prof) + length(p2_prof)-1),]
        }
        # Store
        nprnt1[row_for_new_1:(row_for_new_1+999),] <- new_prnt1
        nprnt2[row_for_new_2:(row_for_new_2+999),] <- new_prnt2
        
        row_for_new_1 <- row_for_new_1 + 1000
        row_for_new_2 <- row_for_new_2 + 1000
        
      }
    }
  }
  
  # remove any zero profiles
  nprnt1 <- nprnt1[which(nprnt1$freq>0),]
  nprnt2 <- nprnt2[which(nprnt2$freq>0),]
  
  return(list(Pnew = nprnt1, Qnew = nprnt2))
}



################################################################################################################################################################################################################################
########################################################************************ LOSS *******************************################################################################################################################
################################################################################################################################################################################################################################

loss <- function(population, gamma_in){
  # population: bugs that are losing MGE
  # gamma_in: rate of loss for each of the 10 MGE in this system
  lp2 <- c()
  
  for(j in 1:dim(population)[1]){ # for each profile in the population
    
    lph <- matrix(0,population[j,"freq"],10)
    
    for(ii in 1:10){ # for each MGE
      
      if(population[j,ii] > 0){ # so loss can happen
        ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
        if(population[j,"freq"] * max(gamma_in) < 5){
          lph[,ii] <- rbinom(n = population[j,"freq"], size = 1, p = gamma_in[ii])
        }else{
          nmovements = rnorm(n = 1, mean = population[j,"freq"]*gamma_in[ii],sd = sqrt(population[j,"freq"] * gamma_in[ii] * (1 - gamma_in[ii])))
          if(nmovements > 0){
            lph[runif(nmovements, 1, population[j,"freq"]), ii] <- 1 # assign randomly to each of the population[j,"freq"] bugs
          }
        }
      }
      
    }
    if(sum(lph)>0){ # only do if some to lose
      lph <- as.data.frame(lph)
      colnames(lph) = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")
      
      # Count profiles
      lp <- plyr::count(lph, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")) %>% 
        mutate(Total = select(., v1:v10) %>% rowSums()) %>% 
        filter(Total>0) %>% select(-Total) # remove those with zero counts
      
      # If any, then take them off the population 
      if(dim(lp)[1]>0){
        for(ii in 1:dim(lp)[1]){
          lp[ii,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")] <- population[j,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")] - lp[ii,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")]
        }
        lp$j <- j
        lp2 <- rbind(lp2, lp) 
      }
    }
  }
  
  return(lp2)
}


#### MGE based
loss_mge <- function(prnt1, prnt2, gamma_in){
  # prnt/2: parent 1 or 2 profiles
  # mu_in: rate of movement for each of the 10 MGE in this system
  
  d1 <- dim(prnt1)[1]
  
  # total 
  prnt <- rbind(prnt1, prnt2)
  
  nprnt1 <- prnt1; nprnt2 <- prnt2;
  
  # Total number of bugs with each MGE
  nbugs_with_orig = colSums(prnt$freq * prnt[,1:10]) #sum(prnt$freq * prnt[,m])
  total_bugs = sum(prnt$freq)
  n_lost_orig <- nbugs_with_orig * gamma_in # on average
  
  ## Go thru MGE
  for(m in 1:10){  
    n_lost <- n_lost_orig[m]
    nbugs_with = nbugs_with_orig[m]
    
    if(nbugs_with > 0){
      
      # Number lost: 
      nlost <- nbugs_with * as.numeric(gamma_in[m]) # on average
      ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
      if(nlost < 5){
        n_lost <- sum(rbinom(n = nlost, size = 1, as.numeric(gamma_in[m]))) # how many transfer = sum of 1s = successful binomial transfer events
      }else{
        n_lost <- round(rnorm(n = 1, mean = nlost,sd = sqrt(nlost * (1 - as.numeric(gamma_in[m])))),0)
      }
      
      # If any transfers
      if(n_lost > 0){
        # Blank to store these transfers
        new_prnt1 <- c(); new_prnt2 <- c();
        # Which nbugs lose the MGE? 
        which_lose = round(runif(n_lost, 1,nbugs_with),0)
        profile_with_ind = intersect(which(prnt[,m] == 1), which(prnt[,"freq"]>0))
        allocated_bins <- table(cut(which_lose, breaks = c(0, cumsum(prnt[profile_with_ind,"freq"]))))
        allocated <- as.numeric(allocated_bins) # how many from each profile?
        
        # Build new profiles
        for(i in 1:length(profile_with_ind)){
          #print(c(n_lost, allocated[i]))
          if(allocated[i]>0){ # if allocate any to this profile
            if(profile_with_ind[i] <= d1){ # if less than d1 then in first parent
              profile <- prnt[profile_with_ind[i],]
              profile[m] <- 0
              profile["freq"] <- allocated[i]
              # New profiles come from old ones
              nprnt1[profile_with_ind[i],"freq"] <- nprnt1[profile_with_ind[i],"freq"] - allocated[i]
              new_prnt1 <- rbind(new_prnt1, profile)
            }else{
              profile <- prnt[profile_with_ind[i],]
              profile[m] <- 0
              profile["freq"] <- allocated[i]
              nprnt2[profile_with_ind[i]-d1,"freq"] <- nprnt2[profile_with_ind[i]-d1,"freq"] - allocated[i]
              new_prnt2 <- rbind(new_prnt2, profile)
            }
          }
        }
        
        nprnt1 <- rbind(nprnt1, new_prnt1)
        nprnt2 <- rbind(nprnt2, new_prnt2)
        # remove any zero profiles
        nprnt1 <- nprnt1[which(nprnt1$freq>0),]
        nprnt2 <- nprnt2[which(nprnt2$freq>0),]
        # New population
        d1 <- dim(nprnt1)[1]
      }
    }
    
  }
  
  return(list(Pnew = nprnt1, Qnew = nprnt2))
}


#### MGE based 
# try to do all at start - sped up by looking only at positive values / indices
loss_mge_all <- function(prnt1, prnt2, gamma_in){
  # prnt/2: parent 1 or 2 profiles
  # mu_in: rate of movement for each of the 10 MGE in this system
  
  d1 <- dim(prnt1)[1]
  d2 <- dim(prnt2)[1]
  
  # total 
  prnt <- rbind(prnt1, prnt2)
  
  # New matrices
  nprnt1 <-  matrix(0, nrow = 15000, ncol = dim(prnt1)[2])
  colnames(nprnt1) <- colnames(prnt1)
  nprnt2 <-  matrix(0, nrow = 15000, ncol = dim(prnt1)[2])
  colnames(nprnt2) <- colnames(prnt1)
  
  nprnt1 <- rbind(prnt1, nprnt1) # if allocate get list conversion
  nprnt2 <- rbind(prnt2, nprnt2)
  
  row_for_new_1 <- d1 + 1
  row_for_new_2 <- d2 + 1 
  
  # Blank matrix 
  blank <- matrix(0,1000,dim(prnt1)[2]); colnames(blank) <- colnames(prnt1); blank <- as.data.frame(blank)
  
  # Total number of bugs with each MGE
  nbugs_with_orig = colSums(prnt$freq * prnt[,1:10]) #sum(prnt$freq * prnt[,m])
  total_bugs = sum(prnt$freq)
  n_lost_orig <- nbugs_with_orig * gamma_in 
  
  if(sum(n_lost_orig>0)){ # if have to do anything (average chance - need to check)
    
    # which > 0? 
    pos <- which(n_lost_orig > 0)
    pos_norm <- which(n_lost_orig > 5) # these are normally sampled
    l_norm <- length(pos_norm)
    pos_binom <- setdiff(pos, pos_norm) # there are binomally sampled
    l_binom <- length(pos_binom)
    
    # Number lost: 
    n_lost <- matrix(0,1,10)
    
    ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
    if(l_binom > 0){ # if any np / nq < 5
      size_b = round(as.numeric(n_lost_orig[pos_binom]),0) # already gamma in n_lost_orig
      n_lost_binom <- rbinom( n = l_binom, size = size_b, prob = c(as.numeric(gamma_in[pos_binom])))
      n_lost[pos_binom] <- n_lost_binom # store
    }
    
    if(l_norm > 0){ # if any np / nq > 5
      size_n = round(as.numeric(n_lost_orig[pos_norm]),0) # already gamma in n_lost_orig
      p_n = sqrt(size_n * (1-as.numeric(gamma_in[pos_norm])))
      if(l_norm > 1){
        n_lost_norm <- round(rmvnorm( n = 1, mean = size_n, sigma = diag(p_n) ),0)
      } else {
        n_lost_norm <- round(rnorm( n = 1, mean = size_n, sd = p_n ),0)
      }
      n_lost[pos_norm] <- n_lost_norm # store
    }
    
    
    # If any transfers
    if(sum(n_lost) > 0){
      
      # which > 0? 
      pos <- which(n_lost > 0)
      
      for(m in pos){
        # Blank to store these transfers
        new_prnt1 <- blank; new_prnt2 <- blank;
        
        index1 <- 1; index2 <- 1
        # Which nbugs lose the MGE? 
        which_lose = round(runif(n_lost[m], 1, nbugs_with_orig[m]),0)
        profile_with_ind = intersect(which(prnt[,m] == 1), which(prnt[,"freq"]>0))
        allocated_bins <- table(cut(which_lose, breaks = c(0, cumsum(prnt[profile_with_ind,"freq"]))))
        allocated <- as.numeric(allocated_bins) # how many from each profile?
        
        alloc_pos <- which(allocated > 0) # only look at those that are positive
        p1_pos <- which(profile_with_ind[alloc_pos] <= d1) # parent 1
        p2_pos <- setdiff(seq(1,length(alloc_pos),1), p1_pos) # parent 2
        
        # Build new profiles
        #print(c(n_lost, allocated[i]))
        for(i in alloc_pos[p1_pos]){
          profile <- prnt[profile_with_ind[i],]
          profile[c(m,"freq")] <- c(0,allocated[i])
          # New profiles come from old ones
          nprnt1[profile_with_ind[i],"freq"] <- nprnt1[profile_with_ind[i],"freq"] - allocated[i]
          new_prnt1[index1,] <- profile
          index1 <- index1 + 1
        }
        for(i in alloc_pos[p2_pos]){
          profile <- prnt[profile_with_ind[i],]
          profile[m] <- 0
          profile["freq"] <- allocated[i]
          nprnt2[profile_with_ind[i]-d1,"freq"] <- nprnt2[profile_with_ind[i]-d1,"freq"] - allocated[i]
          new_prnt2[index2,] <- profile
          index2 <- index2 + 1
        }
        
        # Store
        nprnt1[row_for_new_1:(row_for_new_1+999),] <- new_prnt1
        nprnt2[row_for_new_2:(row_for_new_2+999),] <- new_prnt2
        
        row_for_new_1 <- row_for_new_1 + 1000
        row_for_new_2 <- row_for_new_2 + 1000
      }
    }
  }
  
  # remove any zero profiles
  nprnt1 <- nprnt1[which(nprnt1$freq>0),]
  nprnt2 <- nprnt2[which(nprnt2$freq>0),]
  
  return(list(Pnew = nprnt1, Qnew = nprnt2))
}


#### MGE based
### Trying to remove rbind - didn't work!
loss_mge_new <- function(prnt1, prnt2, gamma_in){
  # prnt/2: parent 1 or 2 profiles
  # gamma_in: rate of movement for each of the 10 MGE in this system
  
  # Size of parent 1
  d1 <- dim(prnt1)[1]
  
  # total 
  prnt <- rbind(prnt1, prnt2)
  
  # New matrices
  new <-  matrix(0, nrow = 15000, ncol = dim(prnt1)[2])
  colnames(new) <- colnames(prnt1)
  nprnt1 <- rbind(prnt1,new) # 10,000 for each new for each mge + 5,000 for initial
  nprnt2 <- rbind(prnt2,new)
  
  row_for_new_1 <- dim(prnt1)[1] + 1
  row_for_new_2 <- dim(prnt2)[1] + 1 
  
  # Blank matrix 
  blank <- matrix(0,1000,dim(prnt1)[2]); colnames(blank) <- colnames(prnt1); blank <- as.data.frame(blank)
  
  # Total number of bugs with each MGE
  nbugs_with_orig = colSums(prnt$freq * prnt[,1:10]) #sum(prnt$freq * prnt[,m])
  total_bugs = sum(prnt$freq)
  
  # Number lost: 
  n_lost_orig <- nbugs_with * gamma_in # on average
  
  ## Go thru MGE
  for(m in 1:10){  
    n_lost <- n_lost_orig[m]
    
    ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
    if(n_lost < 5){
      n_lost <- sum(rbinom(n = nlost, size = 1, gamma_in[m])) # how many transfer = sum of 1s = successful binomial transfer events
    }else{
      n_lost <- round(rnorm(n = 1, mean = nlost*gamma_in[m],sd = sqrt(nlost * gamma_in[m] * (1 - gamma_in[m]))),0)
    }
    
    nbugs_with <- nbugs_with_orig[m]
    
    if(nbugs_with > 0){  
      # If any transfers
      if(n_lost > 0){
        # Blank to store these transfers
        new_prnt1 <- blank; new_prnt2 <- blank;
        
        index1 <- 1; index2 <- 1
        # Which nbugs lose the MGE? 
        which_lose = round(runif(n_lost, 1,nbugs_with),0)
        profile_with_ind = intersect(which(prnt[,m] == 1), which(prnt[,"freq"]>0))
        allocated_bins <- table(cut(which_lose, breaks = c(0, cumsum(prnt[profile_with_ind,"freq"]))))
        allocated <- as.numeric(allocated_bins) # how many from each profile?
        
        # Build new profiles
        for(i in 1:length(profile_with_ind)){
          #print(c(n_lost, allocated[i]))
          if(allocated[i]>0){ # if allocate any to this profile
            if(profile_with_ind[i] <= d1){ # if less than d1 then in first parent
              profile <- prnt[profile_with_ind[i],]
              profile[m] <- 0
              profile["freq"] <- allocated[i]
              # New profiles come from old ones
              nprnt1[profile_with_ind[i],"freq"] <- nprnt1[profile_with_ind[i],"freq"] - allocated[i]
              new_prnt1[index1,] <- profile
              index1 <- index1 + 1
            }else{
              profile <- prnt[profile_with_ind[i],]
              profile[m] <- 0
              profile["freq"] <- allocated[i]
              nprnt2[profile_with_ind[i]-d1,"freq"] <- nprnt2[profile_with_ind[i]-d1,"freq"] - allocated[i]
              new_prnt2[index2,] <- profile
              index2 <- index2 + 1
            }
          }
        }
        
        # Store
        nprnt1[row_for_new_1:(row_for_new_1+999),] <- new_prnt1
        nprnt2[row_for_new_2:(row_for_new_2+999),] <- new_prnt2
        
        row_for_new_1 <- row_for_new_1 + 1000
        row_for_new_2 <- row_for_new_2 + 1000
      }
    }
  }
  
  # remove any zero profiles
  nprnt1 <- nprnt1[which(nprnt1$freq>0),]
  nprnt2 <- nprnt2[which(nprnt2$freq>0),]
  
  return(list(Pnew = nprnt1, Qnew = nprnt2))
}

#### MGE based 
# try to do all at start - sped up by looking only at positive values / indices
# and in the allocation?
loss_mge_all_new <- function(prnt1, prnt2, gamma_in){
  # prnt/2: parent 1 or 2 profiles
  # mu_in: rate of movement for each of the 10 MGE in this system
  
  d1 <- dim(prnt1)[1]
  d2 <- dim(prnt2)[1]
  
  # total 
  prnt <- rbind(prnt1, prnt2)
  
  # New matrices
  nprnt1 <-  matrix(0, nrow = 15000, ncol = dim(prnt1)[2])
  colnames(nprnt1) <- colnames(prnt1)
  nprnt2 <-  matrix(0, nrow = 15000, ncol = dim(prnt1)[2])
  colnames(nprnt2) <- colnames(prnt1)
  
  nprnt1 <- rbind(prnt1, nprnt1) # if allocate get list conversion
  nprnt2 <- rbind(prnt2, nprnt2)
  
  row_for_new_1 <- d1 + 1
  row_for_new_2 <- d2 + 1 
  
  # Blank matrix 
  blank <- matrix(0,1000,dim(prnt1)[2]); colnames(blank) <- colnames(prnt1); blank <- as.data.frame(blank)
  
  # Total number of bugs with each MGE
  nbugs_with_orig = colSums(prnt$freq * prnt[,1:10]) #sum(prnt$freq * prnt[,m])
  total_bugs = sum(prnt$freq)
  n_lost_orig <- nbugs_with_orig * gamma_in 
  
  if(sum(n_lost_orig>0)){ # if have to do anything (average chance - need to check)
    
    # which > 0? 
    pos <- which(n_lost_orig > 0)
    pos_norm <- which(n_lost_orig > 5) # these are normally sampled
    l_norm <- length(pos_norm)
    pos_binom <- setdiff(pos, pos_norm) # there are binomally sampled
    l_binom <- length(pos_binom)
    
    # Number lost: 
    n_lost <- matrix(0,1,10)
    
    ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
    if(l_binom > 0){ # if any np / nq < 5
      size_b = round(as.numeric(n_lost_orig[pos_binom]),0) # already gamma in n_lost_orig
      n_lost_binom <- rbinom( n = l_binom, size = size_b, prob = c(as.numeric(gamma_in[pos_binom])))
      n_lost[pos_binom] <- n_lost_binom # store
    }
    
    if(l_norm > 0){ # if any np / nq > 5
      size_n = round(as.numeric(n_lost_orig[pos_norm]),0) # already gamma in n_lost_orig
      p_n = sqrt(size_n * (1-as.numeric(gamma_in[pos_norm])))
      if(l_norm > 1){
        n_lost_norm <- round(rmvnorm( n = 1, mean = size_n, sigma = diag(p_n) ),0)
      } else {
        n_lost_norm <- round(rnorm( n = 1, mean = size_n, sd = p_n ),0)
      }
      n_lost[pos_norm] <- n_lost_norm # store
    }
    
    
    # If any transfers
    if(sum(n_lost) > 0){
      
      # which > 0? 
      pos <- which(n_lost > 0)
      
      for(m in pos){
        # Blank to store these transfers
        new_prnt1 <- blank; new_prnt2 <- blank;
        
        index1 <- 1; index2 <- 1
        # Which nbugs lose the MGE? 
        which_lose = round(runif(n_lost[m], 1, nbugs_with_orig[m]),0)
        profile_with_ind = intersect(which(prnt[,m] == 1), which(prnt[,"freq"]>0))
        allocated_bins <- table(cut(which_lose, breaks = c(0, cumsum(prnt[profile_with_ind,"freq"]))))
        allocated <- as.numeric(allocated_bins) # how many from each profile?
        
        alloc_pos <- which(allocated > 0) # only look at those that are positive
        p1_pos <- which(profile_with_ind[alloc_pos] <= d1) # parent 1
        p2_pos <- setdiff(seq(1,length(alloc_pos),1), p1_pos) # parent 2
        
        # Build new profiles
        #print(c(n_lost, allocated[i]))
        if(length(p1_pos) > 0){
          # parent 1
          profiles <- prnt[profile_with_ind[alloc_pos[p1_pos]],]
          profiles[,m] <- 0
          profiles[,"freq"] <- allocated[alloc_pos[p1_pos]]
          
          nprnt1[profile_with_ind[alloc_pos[p1_pos]],"freq"] <- nprnt1[profile_with_ind[alloc_pos[p1_pos]],"freq"] - allocated[alloc_pos[p1_pos]]
          prnt[profile_with_ind[alloc_pos[p1_pos]],"freq"] <- prnt[profile_with_ind[alloc_pos[p1_pos]],"freq"] - allocated[alloc_pos[p1_pos]] # take from original too else can get additional
          
          new_prnt1[index1:(index1+length(p1_pos)-1),] <- profiles
          index1 <- index1+length(p1_pos)
        } 
        # parent 2
        if(length(p2_pos)>0){
          profiles <- prnt[profile_with_ind[alloc_pos[p2_pos]],]
          profiles[,m] <- 0
          profiles[,"freq"] <- allocated[alloc_pos[p2_pos]]
          
          nprnt2[(profile_with_ind[alloc_pos[p2_pos]] - d1),"freq"] <- nprnt2[(profile_with_ind[alloc_pos[p2_pos]]-d1),"freq"] - allocated[alloc_pos[p2_pos]]
          prnt[(profile_with_ind[alloc_pos[p2_pos]]),"freq"] <- prnt[(profile_with_ind[alloc_pos[p2_pos]]),"freq"] - allocated[alloc_pos[p2_pos]] # not "- d1" as all tog in prnt
          
          new_prnt2[index2:(index2+length(p2_pos)-1),] <- profiles
          index2 <- index2+length(p2_pos)
        }
        
        # Store
        nprnt1[row_for_new_1:(row_for_new_1+999),] <- new_prnt1
        nprnt2[row_for_new_2:(row_for_new_2+999),] <- new_prnt2
        
        row_for_new_1 <- row_for_new_1 + 1000
        row_for_new_2 <- row_for_new_2 + 1000
      }
    }
  }
  
  # remove any zero profiles
  nprnt1 <- nprnt1[which(nprnt1$freq>0),]
  nprnt2 <- nprnt2[which(nprnt2$freq>0),]
  
  return(list(Pnew = nprnt1, Qnew = nprnt2))
}

################################################################################################################################################################################################################################
########################################################************************ RUN *******************************################################################################################################################
################################################################################################################################################################################################################################

run_sim <- function(tsteps, parameters){
  # mu_in = gain rates
  # gamma_in = loss rates
  # growth_in = fitness cost
  # grate_in = underlying growth rate
  
  # If have rates for all 10 elements
  if(length(parameters) == 31){
  mu_in <- as.numeric(parameters[1:10])
  gamma_in <- as.numeric(parameters[11:20] )
  growth_in <- as.numeric(parameters[21:30])
  grate_in <- as.numeric(parameters[31])
  }
  
  # If just looking at the elements that move
  if(length(parameters) < 31 && length(parameters) > 7){
    mu_in <- as.numeric(c(0,parameters["mu2"],0,0,parameters["mu5"],parameters["mu6"],parameters["mu7"],parameters["mu8"],0,parameters["mu10"]))
    gamma_in <- as.numeric(c(0,parameters["gamma2"],0,0,parameters["gamma5"],parameters["gamma6"],parameters["gamma7"],parameters["gamma8"],0,parameters["gamma10"]))
    growth_in <- as.numeric(c(0,parameters["f2"],0,0,parameters["f5"],parameters["f6"],parameters["f7"],parameters["f8"],0,parameters["f10"]))
    grate_in <- as.numeric(parameters["grow"])
  }
  
  # If fixed input - same rates for all elements
  if(length(parameters) == 4){
    mu_in <- as.numeric(c(0,parameters["mu"],0,0,parameters["mu"],parameters["mu"],parameters["mu"],parameters["mu"],0,parameters["mu"]))
    gamma_in <- as.numeric(c(0,parameters["gamma"],0,0,parameters["gamma"],parameters["gamma"],parameters["gamma"],parameters["gamma"],0,parameters["gamma"]))
    growth_in <- as.numeric(c(0,parameters["f"],0,0,parameters["f"],parameters["f"],parameters["f"],parameters["f"],0,parameters["f"]))
    grate_in <- as.numeric(parameters["grow"])
  }
  
  # If fixed input - same rates for phage vs plasmids
  if(length(parameters) == 7){
    #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
    mu_in <- as.numeric(c(0,parameters["mu_phage"],0,0,parameters["mu_phage"],parameters["mu_plasmid"],parameters["mu_plasmid"],parameters["mu_plasmid"],0,parameters["mu_plasmid"]))
    gamma_in <- as.numeric(c(0,parameters["gamma_phage"],0,0,parameters["gamma_phage"],parameters["gamma_plasmid"],parameters["gamma_plasmid"],parameters["gamma_plasmid"],0,parameters["gamma_plasmid"]))
    growth_in <- as.numeric(c(0,parameters["f_phage"],0,0,parameters["f_phage"],parameters["f_plasmid"],parameters["f_plasmid"],parameters["f_plasmid"],0,parameters["f_plasmid"]))
    grate_in <- as.numeric(parameters["grow"])
  }
  

  
  # carrying capacity
  K = 0.8 * 10^7
  
  ## initial conditions
  #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
  Pinit = t(c(0,1,0,0,1,1,1,1,0,0))
  Qinit = t(c(0,0,0,0,0,0,0,0,0,1))
  Pinit <- as.data.frame(Pinit)
  colnames(Pinit) = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")
  Qinit <- as.data.frame(Qinit)
  colnames(Qinit) = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")
  
  ## State matrix 
  Pinit$freq <- 7 * 10^2 # start with all at P1
  Qinit$freq <- 7 * 10^2 # start with all at Q1
  
  if(length(mu_in) < 10){stop("Not enough gain parameters")}
  if(length(gamma_in) < 10){stop("Not enough loss parameters")}
  if(length(growth_in) < 10){stop("Not enough fitness parameters")}
  if(sum(growth_in) > 1){stop("Fitness cost too large", return(list(P_all = Pinit, Q_all = Qinit, error = "", 
                                                                    prev_predict = c(0), totl_predict = c(0))))}
  
  # At start
  P = Pinit
  Q = Qinit
  P_all <- c()
  Q_all <- c()
  
  for(iit in 1:tsteps){
    #print(iit)
    
    ####### Events
    nbugs = sum(P$freq) + sum(Q$freq)  # total bugs
    
    ######################## (1) gain 
    #print("gain")
    new_after_gain <- gain_mge_all_new(P,Q,mu_in)
    P <- new_after_gain$Pnew
    Q <- new_after_gain$Qnew
    
    ######################## (2) loss
    #print("loss")
    new_after_loss <- loss_mge_all_new(P,Q,gamma_in)
    Ploss <- new_after_loss$Pnew
    Qloss <- new_after_loss$Qnew
    
    # Pull together - may be multiple rows with same profiles
    P <- plyr::count(Ploss, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"),wt_var = "freq")
    Q <- plyr::count(Qloss, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"),wt_var = "freq")
    
    ######################### (3) Growth & Death
    #print("growth")
    
    gp3 <- matrix(0,1,dim(P)[1]); gq3 <- matrix(0,1,dim(Q)[1]);
    # Add in fitness column for each profile (fitness cost = min(1,sum(P[j,1:10] * growth_in)))
    # Each can have a fitness cost: assume additive up to totally unfit = 1
    # If no element (v# == 0) then no fitness cost conferred 
    P <- P %>% rowwise() %>% mutate(fitness_cost = sum(v1*growth_in[1] + v2*growth_in[2] + v3*growth_in[3] + 
                                                         v4*growth_in[4] + v5*growth_in[5] + v6*growth_in[6] + 
                                                         v7*growth_in[7] + v8*growth_in[8] + v9*growth_in[9] + v10*growth_in[10]))
    Q <- Q %>% rowwise() %>% mutate(fitness_cost = sum(v1*growth_in[1] + v2*growth_in[2] + v3*growth_in[3] + 
                                                         v4*growth_in[4] + v5*growth_in[5] + v6*growth_in[6] + 
                                                         v7*growth_in[7] + v8*growth_in[8] + v9*growth_in[9] + v10*growth_in[10]))
    
    # Add in probability survive
    cps = (1-nbugs/K)
    P <- P %>% rowwise() %>% mutate(prob_survive = 0.5 * (1 + (1-fitness_cost) * grate_in * cps)) %>% ungroup()
    Q <- Q %>% rowwise() %>% mutate(prob_survive = 0.5 * (1 + (1-fitness_cost) * grate_in * cps)) %>% ungroup()
    # Correct list / tibble making: stops loss working
    P <- as.data.frame(P)
    Q <- as.data.frame(Q)
    
    ## For P 
    # which > 0? 
    survivors = P$freq * P$prob_survive
    pos <- which(survivors > 0)
    pos_norm <- which(survivors > 5) # these are normally sampled
    l_norm <- length(pos_norm)
    pos_binom <- setdiff(pos, pos_norm) # there are binomally sampled
    l_binom <- length(pos_binom)
    
    ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
    if(l_binom > 0){ # if any np / nq < 5
      size_b = 2 * round(as.numeric(P$freq[pos_binom]),0)
      n_grow_binom <- rbinom( n = l_binom, size = size_b, prob = c(as.numeric(P$prob_survive[pos_binom])))
      gp3[pos_binom] <- n_grow_binom # store
    }
    
    if(l_norm > 0){ # if any np / nq > 5
      size_n = 2 * round(as.numeric(P$prob_survive[pos_norm]) * as.numeric(P$freq[pos_norm]),0)
      p_n = sqrt(size_n * (1-as.numeric(P$prob_survive[pos_norm])))
      if(l_norm > 1){
        n_grow_norm <- round(rmvnorm( n = 1, mean = size_n, sigma = diag(p_n) ),0)
      } else {
        n_grow_norm <- round(rnorm( n = 1, mean = size_n, sd = p_n ),0)
      }
      gp3[pos_norm] <- n_grow_norm # store
    }
    
    # For Q
    # which > 0? 
    survivors = Q$freq * Q$prob_survive
    pos <- which(survivors > 0)
    pos_norm <- which(survivors > 5) # these are normally sampled
    l_norm <- length(pos_norm)
    pos_binom <- setdiff(pos, pos_norm) # there are binomally sampled
    l_binom <- length(pos_binom)
    
    ## If X ~ B(n, p) and if np or nq > 5, then X is approximately N(np, npq)
    if(l_binom > 0){ # if any np / nq < 5
      size_b = 2 * round(as.numeric(Q$freq[pos_binom]),0)
      n_grow_binom <- rbinom( n = l_binom, size = size_b, prob = c(as.numeric(Q$prob_survive[pos_binom])))
      gq3[pos_binom] <- n_grow_binom # store
    }
    
    if(l_norm > 0){ # if any np / nq > 5
      if(l_norm > 1){
        size_n = 2 * round(as.numeric(Q$prob_survive[pos_norm]) * as.numeric(Q$freq[pos_norm]),0)
        p_n = sqrt(size_n * (1-as.numeric(Q$prob_survive[pos_norm])))
        n_grow_norm <- round(rmvnorm( n = 1, mean = size_n, sigma = diag(p_n) ),0)
      } else {
        size_n = 2 * round(as.numeric(Q$prob_survive[pos_norm]) * as.numeric(Q$freq[pos_norm]),0)
        p_n = sqrt(size_n * (1-as.numeric(Q$prob_survive[pos_norm])))
        n_grow_norm <- round(rnorm( n = 1, mean = size_n, sd = p_n ),0)
      }
      gq3[pos_norm] <- n_grow_norm # store
    }
    
    ######################### Totals 
    P$freq <- as.numeric(gp3)
    Q$freq <- as.numeric(gq3)
    
    # # # Remove any profiles which no longer have any bugs with this profile
    P <- P[which(P$freq>0),]
    Q <- Q[which(Q$freq>0),]
    
    #########################  Store
    P_new <- P
    P_new$time <- iit
    Q_new <- Q
    Q_new$time <- iit
    
    P_all <- rbind(P_all, P_new)
    Q_all <- rbind(Q_all, Q_new)
    
    ### CHECKS
    error <- 0
    # (1) check grate not too quick: if 107 before 3 days or after 5 days then too slow
    total_bugs = sum(P$freq) + sum(Q$freq)
    if(total_bugs > 0.9 * K && iit < 72){error <- 11; break;} #print("Growth rate too fast"); 
    if(total_bugs < 0.9 * K && iit > 210){error <- 12;  break} #print("Growth rate too slow"); 
    
    # # (2) check not too many profiles: if more 100 stop!
    #print(dim(P))
    #print(dim(Q))
    #if(dim(P)[1] > 100){error <- 21; print("P: too many profiles"); break}
    #if(dim(Q)[1] > 100){error <- 22; print("Q: too many profiles"); break}
    
    # (3) if more than 15 at 5% then stop 
    if(dim(P)[1]>15){ # if more than 15 at 5% in general then check next calcs
      prop <- 100 * P$freq / sum(P$freq)
      n_15 = length(which(prop > 5))
      if(n_15 > 15){error <- 31;  break}#print("P: too many over 5%");
    }
    if(dim(Q)[1]>15){ # if more than 15 at 5% in general then check next calcs
      qprop <- 100 * Q$freq / sum(Q$freq)
      n_15 = length(which(qprop > 5))
      if(n_15 > 15){error <- 32;  break} #print("Q: too many over 5%");
    }
    # (4) if too few after 24hrs than stop 
    if(iit > 24){
      if(dim(P)[1]< 3){error <- 41;  break} # if fewer than 3 profiles at 2 days then stop
      if(dim(Q)[1]< 3){error <- 42; break} # if fewer than 3 profiles at 2 days then stop
    }
    
    # (5) Need both populations to grow evenly "colonising equally well" after 2 days
    if(iit > 48){
      if(sum(P$freq)/total_bugs < 0.01){error <- 51;  break} #print("P: less than 40%");
      if(sum(P$freq)/total_bugs > 0.99){error <- 52;  break} #print("P: greater than 40%");
    }
    
  }
  
  # Label profiles
  P_all$label <- getLabels(P_all[,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")])
  Q_all$label <- getLabels(Q_all[,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")])
  
  ## Convert to prevalence over time
  if(iit == tsteps){ # if reach end timepoint without error
    colnames(P_all) <- c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4","freq","fitness","prob_survive","time","label")
    colnames(Q_all) <- c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4","freq","fitness","prob_survive","time","label")
    pnew <- P_all %>% pivot_longer(cols = c("SCCmec":"p4")) %>% 
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% 
      mutate(total = sum(freq)/10) %>% group_by(time, name) %>% # /10 for total as 10 elements
      mutate(nbugs_with = value * freq) %>% 
      summarise(nbugs_with_total = sum(nbugs_with), total = min(total), # min here but could be max or anything - just need to move total into summary
                prev = ifelse(nbugs_with_total == 0, 0, nbugs_with_total / total),.groups = "drop") %>%
      select(time,name,prev) %>% mutate(parent = 1)
    
    qnew <-  Q_all %>% pivot_longer(cols = c("SCCmec":"p4")) %>% 
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% 
      mutate(total = sum(freq)/10) %>% group_by(time, name) %>% # /10 for total as 10 elements
      mutate(nbugs_with = value * freq) %>% 
      summarise(nbugs_with_total = sum(nbugs_with), total = min(total), # min here but could be max or anything - just need to move total into summary
                prev = ifelse(nbugs_with_total == 0, 0, nbugs_with_total / total),.groups = "drop") %>%
      select(time,name,prev) %>% mutate(parent = 2)
    
    # Totals 
    ptotals <- P_all %>% select(time,freq) %>%
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% 
      summarise(total = sum(freq),.groups = "drop") %>% mutate(parent = 1)
    qtotals <- Q_all %>% select(time,freq) %>%
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% 
      summarise(total = sum(freq),.groups = "drop") %>% mutate(parent = 2)
    
    # to match data
    prev_predict <- rbind(pnew,qnew)
    totl_predict <- rbind(ptotals, qtotals)
  } else {prev_predict <- c(); totl_predict <- c()}
  # Output
  return(list(P_all = P_all, Q_all = Q_all, error = error, 
              prev_predict = prev_predict, totl_predict = totl_predict))
}


run_sim_old <- function(tsteps, parameters){
  # mu_in = gain rates
  # gamma_in = loss rates
  # growth_in = fitness cost
  # grate_in = underlying growth rate
  
  mu_in <- as.numeric(parameters[1:10])
  gamma_in <- as.numeric(parameters[11:20] )
  growth_in <- as.numeric(parameters[21:30])
  grate_in <- as.numeric(parameters[31])
  
  # If just looking at the elements that move
  if(length(parameters) < 31){
    mu_in <- as.numeric(c(0,parameters["mu2"],0,0,parameters["mu5"],parameters["mu6"],parameters["mu7"],parameters["mu8"],0,parameters["mu10"]))
    gamma_in <- as.numeric(c(0,parameters["gamma2"],0,0,parameters["gamma5"],parameters["gamma6"],parameters["gamma7"],parameters["gamma8"],0,parameters["gamma10"]))
    growth_in <- as.numeric(c(0,parameters["f2"],0,0,parameters["f5"],parameters["f6"],parameters["f7"],parameters["f8"],0,parameters["f10"]))
    grate_in <- as.numeric(parameters["grow"])
  }
  
 
  
  # carrying capacity
  K = 0.8 * 10^7
  
  ## initial conditions
  #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
  Pinit = t(c(0,1,0,0,1,1,1,1,0,0))
  Qinit = t(c(0,0,0,0,0,0,0,0,0,1))
  Pinit <- as.data.frame(Pinit)
  colnames(Pinit) = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")
  Qinit <- as.data.frame(Qinit)
  colnames(Qinit) = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")
  
  ## State matrix 
  Pinit$freq <- 7 * 10^2 # start with all at P1
  Qinit$freq <- 7 * 10^2 # start with all at Q1
  
  if(length(mu_in) < 10){stop("Not enough gain parameters")}
  if(length(gamma_in) < 10){stop("Not enough loss parameters")}
  if(length(growth_in) < 10){stop("Not enough fitness parameters")}
  if(sum(growth_in) > 1){stop("Fitness cost too large", return(list(P_all = Pinit, Q_all = Qinit, error = "error", 
                                                                    prev_predict = c(0), totl_predict = c(0))))}
  
  # At start
  P = Pinit
  Q = Qinit
  P_all <- c()
  Q_all <- c()
  
  for(iit in 1:tsteps){
    #print(iit)
    
    ####### Events
    nbugs = sum(P$freq) + sum(Q$freq)  # total bugs
    
    ######################## (1) gain 
    # print("gain")
    new_after_gain <- gain_mge(P,Q,mu_in)
    P <- new_after_gain$Pnew
    Q <- new_after_gain$Qnew
    
    ######################## (2) loss
    # print("loss")
    new_after_loss <- loss_mge(P,Q,gamma_in)
    Ploss <- new_after_loss$Pnew
    Qloss <- new_after_loss$Qnew
    
    # Pull together - may be multiple rows with same profiles
    P <- plyr::count(Ploss, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"),wt_var = "freq")
    Q <- plyr::count(Qloss, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"),wt_var = "freq")
    
    ######################### (3) Growth & Death
    # print("growth")
    
    gp3 <- c(); gq3 <- c();
    # Add in fitness column for each profile (fitness cost = min(1,sum(P[j,1:10] * growth_in)))
    # Each can have a fitness cost: assume additive up to totally unfit = 1
    # If no element (v# == 0) then no fitness cost conferred 
    P <- P %>% rowwise() %>% mutate(fitness_cost = sum(v1*growth_in[1] + v2*growth_in[2] + v3*growth_in[3] + 
                                                         v4*growth_in[4] + v5*growth_in[5] + v6*growth_in[6] + 
                                                         v7*growth_in[7] + v8*growth_in[8] + v9*growth_in[9] + v10*growth_in[10]))
    Q <- Q %>% rowwise() %>% mutate(fitness_cost = sum(v1*growth_in[1] + v2*growth_in[2] + v3*growth_in[3] + 
                                                         v4*growth_in[4] + v5*growth_in[5] + v6*growth_in[6] + 
                                                         v7*growth_in[7] + v8*growth_in[8] + v9*growth_in[9] + v10*growth_in[10]))
    
    # Add in probability survive
    cps = (1-nbugs/K)
    P <- P %>% rowwise() %>% mutate(prob_survive = 0.5 * (1 + (1-fitness_cost) * grate_in * cps)) %>% ungroup()
    Q <- Q %>% rowwise() %>% mutate(prob_survive = 0.5 * (1 + (1-fitness_cost) * grate_in * cps)) %>% ungroup()
    # Correct list / tibble making: stops loss working
    P <- as.data.frame(P)
    Q <- as.data.frame(Q)
    
    
    for(j in 1:dim(P)[1]){ # for each profile in P
      
      ## How many will survive? they double (hence the 2 x P[j,"freq"]) then a proportion die so that on average
      # P * grate * fitness survive: have to have prob_survive < 1 to work.
      if(as.numeric(P[j,"freq"]*P[j,"prob_survive"]) < 6){
        gpnew = sum(rbinom(n = 2 * as.numeric(P[j,"freq"]), size = 1, as.numeric(P[j,"prob_survive"]))) # for each bacteria, probability that survive
      } else {
        gpnew = round(rnorm(n = 1, mean = 2 * as.numeric(P[j,"freq"] * P[j,"prob_survive"]),
                            sd = sqrt(2 * as.numeric(P[j,"freq"]) * as.numeric(P[j,"prob_survive"]) * (1 - as.numeric(P[j,"prob_survive"])))),0)
      }
      gp3 <- c(gp3, gpnew)
    }
    
    for(j in 1:dim(Q)[1]){ # for each profile in Q
      
      ## How many will survive? they double (hence the 2 x Q[j,"freq"]) then a proportion die so that on average
      # Q * grate * fitness survive: have to have prob_survive < 1 to work.
      if(as.numeric(Q[j,"freq"]*Q[j,"prob_survive"] )< 6){
        gqnew = sum(rbinom(n = 2 * as.numeric(Q[j,"freq"]), size = 1, as.numeric(Q[j,"prob_survive"]))) # for each bacteria, probability that transfer happens
      } else {
        gqnew = round(rnorm(n = 1, mean = 2 * as.numeric(Q[j,"freq"] * Q[j,"prob_survive"]),
                            sd = sqrt(2 * as.numeric(Q[j,"freq"] * Q[j,"prob_survive"]) * (1 - as.numeric(Q[j,"prob_survive"])))),0)
      }
      gq3 <- c(gq3, gqnew)
    }
    
    ######################### Totals 
    P$freq <- as.numeric(gp3)
    Q$freq <- as.numeric(gq3)
    
    # # # Remove any profiles which no longer have any bugs with this profile
    P <- P[which(P$freq>0),]
    Q <- Q[which(Q$freq>0),]
    
    #########################  Store
    P_new <- P
    P_new$time <- iit
    Q_new <- Q
    Q_new$time <- iit
    
    P_all <- rbind(P_all, P_new)
    Q_all <- rbind(Q_all, Q_new)
    
    ### CHECKS
    error <- 0
    # (1) check grate not too quick: if 107 before 3 days or after 5 days then too slow
    total_bugs = sum(P$freq) + sum(Q$freq)
    if(total_bugs > 0.9 * K && iit < 72){error <- 11; break;} #print("Growth rate too fast"); 
    if(total_bugs < 0.9 * K && iit > 210){error <- 12;  break} #print("Growth rate too slow"); 
    
    # # (2) check not too many profiles: if more 100 stop!
    #print(dim(P))
    #print(dim(Q))
    #if(dim(P)[1] > 100){error <- 21; print("P: too many profiles"); break}
    #if(dim(Q)[1] > 100){error <- 22; print("Q: too many profiles"); break}
    
    # (3) if more than 15 at 5% then stop 
    if(dim(P)[1]>15){ # if more than 15 at 5% in general then check next calcs
      prop <- 100 * P$freq / sum(P$freq)
      n_15 = length(which(prop > 5))
      if(n_15 > 15){error <- 31;  break}#print("P: too many over 5%");
    }
    if(dim(Q)[1]>15){ # if more than 15 at 5% in general then check next calcs
      qprop <- 100 * Q$freq / sum(Q$freq)
      n_15 = length(which(qprop > 5))
      if(n_15 > 15){error <- 32;  break} #print("Q: too many over 5%");
    }
    # (4) if too few after 24hrs than stop 
    if(iit > 24){
      if(dim(P)[1]< 3){error <- 41;  break} # if fewer than 3 profiles at 2 days then stop
      if(dim(Q)[1]< 3){error <- 42; break} # if fewer than 3 profiles at 2 days then stop
    }
    
    # (5) Need both populations to grow evenly "colonising equally well" after 2 days
    if(iit > 48){
      if(sum(P$freq)/total_bugs < 0.01){error <- 51;  break} #print("P: less than 40%");
      if(sum(P$freq)/total_bugs > 0.99){error <- 52;  break} #print("P: greater than 40%");
    }
    
  }
  
  # Label profiles
  P_all$label <- getLabels(P_all[,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")])
  Q_all$label <- getLabels(Q_all[,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")])
  
  ## Convert to prevalence over time
  if(iit == tsteps){ # if reach end timepoint without error
    colnames(P_all) <- c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4","freq","fitness","prob_survive","time","label")
    colnames(Q_all) <- c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4","freq","fitness","prob_survive","time","label")
    pnew <- P_all %>% pivot_longer(cols = c("SCCmec":"p4")) %>% 
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% 
      mutate(total = sum(freq)/10) %>% group_by(time, name) %>% # /10 for total as 10 elements
      mutate(nbugs_with = value * freq) %>% 
      summarise(nbugs_with_total = sum(nbugs_with), total = min(total), # min here but could be max or anything - just need to move total into summary
                prev = ifelse(nbugs_with_total == 0, 0, nbugs_with_total / total),.groups = "drop") %>%
      select(time,name,prev) %>% mutate(parent = 1)
    
    qnew <-  Q_all %>% pivot_longer(cols = c("SCCmec":"p4")) %>% 
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% 
      mutate(total = sum(freq)/10) %>% group_by(time, name) %>% # /10 for total as 10 elements
      mutate(nbugs_with = value * freq) %>% 
      summarise(nbugs_with_total = sum(nbugs_with), total = min(total), # min here but could be max or anything - just need to move total into summary
                prev = ifelse(nbugs_with_total == 0, 0, nbugs_with_total / total),.groups = "drop") %>%
      select(time,name,prev) %>% mutate(parent = 2)
    
    # Totals 
    ptotals <- P_all %>% select(time,freq) %>%
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% 
      summarise(total = sum(freq),.groups = "drop") %>% mutate(parent = 1)
    qtotals <- Q_all %>% select(time,freq) %>%
      filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% 
      summarise(total = sum(freq),.groups = "drop") %>% mutate(parent = 2)
    
    # to match data
    prev_predict <- rbind(pnew,qnew)
    totl_predict <- rbind(ptotals, qtotals)
  } else {prev_predict <- c(); totl_predict <- c()}
  # Output
  return(list(P_all = P_all, Q_all = Q_all, error = error, 
              prev_predict = prev_predict, totl_predict = totl_predict))
}

### Calculate log posterior from model at theta parameter values
run_sim_logPosterior <- function(theta){
  tsteps = 384
  
  ## Run for these parameters
  out <- run_sim(tsteps, theta)
  max(out$P_all$time)
  
  if(max(out$P_all$time)==tsteps){ # if get to end 
    #### element prevalence from model 
    model_outputp <- out$prev_predict %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
    model_outputp$prev <- round(model_outputp$prev,2)
    
    # Check what distribution of n_colonies at this prevalence in the model pig
    distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, name, prob_all, n_colonies_prev)
    # e.g. to check 
    #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
    
    # lookup the probability from this distribution for the data
    likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","name","n_colonies_prev")) %>% summarise(sum(log(prob_all)))
    
    #### total bugs output from model 
    model_outputt <- out$totl_predict
    
    likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
      as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))
    
    #### Compare to data 
    compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals
    
  }else{compare_dat <- as.numeric(-Inf)}
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


