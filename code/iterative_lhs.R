### Iterative LHS
# Takes in sample and limits - runs for this

lhs_build_run <- function(init_values, limit, save_place, nsamples = 5000){
  # init_values = initial values to build LHS around
  # limit = how far around the init_value
  # save_place = where to 
  
  # Make if need to 
  dir.create(file.path(paste0(here::here(),"/fits"), save_place), showWarnings = FALSE)
  setwd(file.path(paste0(here::here(),"/fits"), save_place))
  
  # What are the ranges? 
  param_ranges <- as.data.frame(cbind(init_values - abs(init_values)/limit,
                                      init_values + abs(init_values)/limit))
  colnames(param_ranges) <- c("min","max")
  
  # transform a Latin hypercube
  nparameters <- nrow(param_ranges)
  X <- randomLHS(nsamples, nparameters) # first = number of samples, 1000 maybe? second = number of parameters, here 2
  Y <- matrix(0, nrow=(nsamples+1), ncol=nparameters)
  # Assume parameters uniformly arranged over the range (could assume normal etc if have evidence...)
  for(ii in 1:nparameters){
    Y[1:nsamples,ii] <- qunif(X[,ii], min = param_ranges[ii,1], max = param_ranges[ii,2])
  }
  Y[(nsamples+1),]<-init_values
  write.csv(Y, "paraset.csv") # save parameter set
  
  #### Parallel
  numCores = 4
  registerDoParallel(numCores)
  
  foreach (ii=1:nsamples) %dopar% {
    
    if(length(Y[1,])==5){
      parameters = c(mu = Y[ii,1],
                     gamma = Y[ii,2],
                     f = Y[ii,3], 
                     grow = Y[ii,4], 
                     rel_fit = Y[ii,5])
    }
    if(length(Y[1,])==8){
      parameters <-  c(mu_phage = Y[ii,1],
                       mu_plasmid = Y[ii,2],
                       gamma_phage = Y[ii,3],
                       gamma_plasmid = Y[ii,4],
                       f_phage = Y[ii,5], 
                       f_plasmid = Y[ii,6], 
                       grow = Y[ii,7], 
                       rel_fit = Y[ii,8])
    }
    if(length(Y[1,])==10){
      parameters <- c(mu = Y[ii,1],
                      gamma = Y[ii,2],
                      f2 = Y[ii,3], f5 = Y[ii,4], f6 = Y[ii,5], 
                      f7 = Y[ii,6], f8 = Y[ii,7], f10 = Y[ii,8], 
                      grow = Y[ii,9], 
                      rel_fit = Y[ii,10])
    }
    
    if(length(Y[1,])==20){
      parameters <- c(mu2 = Y[ii,1], mu5 = Y[ii,2], mu6 = Y[ii,3],
                      mu7 = Y[ii,4], mu8 = Y[ii,5], mu10 = Y[ii,6],
                      gamma2 = Y[ii,7], gamma5 = Y[ii,8],
                      gamma6 = Y[ii,9], gamma7 = Y[ii,10],
                      gamma8 = Y[ii,11], gamma10 = Y[ii,12],
                      f2 = Y[ii,13], f5 = Y[ii,14], f6 = Y[ii,15],
                      f7 = Y[ii,16], f8 = Y[ii,17], f10 = Y[ii,18],
                      grow = Y[ii,19], rel_fit = Y[ii,20])
    }
    
    out <- run_sim_logPosterior(parameters)
    if(abs(out)<10000){write.csv(out,paste0(ii,".csv"))} # don't keep infinite values
  }
  
  stopImplicitCluster()
  
  # Grab those that worked
  which_para_csv10 <- list.files() 
  which_para10 <- as.numeric(sub("\\..*", "",which_para_csv10[-length(which_para_csv10)]))
  
  work10 <-list.files(pattern = "*.csv") %>% 
    map_df(~read_csv(.)) %>% filter(!is.na(x))
  
  if(length(which_para10)>0){
    if(length(which_para10) == 1){
      max_ll_val <- max(work10$x, na.rm = TRUE)
      max_ll_para <- Y[which_para10[which.max(work10$x)],]
      worked <- c(work10$x, Y[which_para10,])}else{
        max_ll_val <- max(work10$x, na.rm = TRUE)
        max_ll_para <- Y[which_para10[which.max(work10$x)],]
        worked <- cbind(work10$x, Y[which_para10,])}
  }else{max_ll_val = c(); max_ll_para = -Inf; worked = c()}
  
  return(list(max_ll_val = max_ll_val, max_ll_para = max_ll_para, worked = worked))
}