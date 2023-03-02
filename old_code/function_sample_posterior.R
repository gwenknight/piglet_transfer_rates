#### Function to sample posterior from trace


posterior_sample_trace <- function(tra, n_samples = 100, break_size = 0.01, plotit = FALSE, plotname = "here"){
  # tra = trace from MCMC
  # n_samples = how many samples from posterior?
  # break_size = size of bins
  # plotit = want to output posterior plot? 
  # plotname = name for plot in plots/
  
  # Posteriors
  n_parameters = dim(tra)[2]
  n_runs = dim(tra)[1]
  prob_vals = matrix(0,n_parameters, length(seq(-1,1.5,break_size))-1)
  for(i in 1:n_parameters){
    prob_vals[i,] <- hist(tra[,i], breaks = seq(-1,1.5,break_size), plot = FALSE)$density#summary(cut(tra[,i], breaks=seq(-1,1,break_size), include.lowest=TRUE, right=FALSE))/n_runs
  }
  colnames(prob_vals) <- hist(tra[,i], breaks = seq(-1,1.5,break_size), plot = FALSE)$mids
  
  if(plotit==TRUE){
    prob_vals <- as.data.frame(prob_vals)
    prob_vals$para_name <- colnames(tra)
    pp <- prob_vals %>% pivot_longer(cols = colnames(prob_vals)[1]:colnames(prob_vals)[ncol(prob_vals)-1])
    pp$name <- as.numeric(pp$name)
    
    g <- ggplot(pp, aes(x=name, y = value)) + geom_line() + facet_wrap(~ para_name, scales = "free")
    ggsave(g, paste0("plots/",plotname,".pdf"))
  }
  
  # Sample
  samples <- matrix(0, n_parameters, n_samples)
  for(i in 1:n_parameters){
    samples[i,] <- sample(hist(tra[,i], breaks = seq(-1,1.5,break_size), plot = FALSE)$mids,100,replace = TRUE, prob = prob_vals[i,])
  }
  
  if(plotit==TRUE){
    sampls <- as.data.frame(samples)
    sampls$para_name <- colnames(tra)
    pp <- sampls %>% pivot_longer(cols = colnames(sampls)[1]:colnames(sampls)[ncol(sampls)-1])
    
    g <- ggplot(pp, aes(x = value)) + geom_histogram(binwidth = break_size) + facet_wrap(~ para_name, scales = "free")
    ggsave(g, paste0("plots/",plotname,"_sample.pdf"))
  }
  
  return(list(probs = prob_vals, vals = colnames(prob_vals), samples = samples))
}

