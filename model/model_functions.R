
#model utility functions

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
  
  return(list(difference_list = difference_list, bacteria = bacteria)) 
}


plot_circles <- function(profiles, output_data, save = F, plot_name = ""){
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
  
  gp <- ggplot(pigg_plotp %>% filter(perc > 5), aes(x0=x_centre, y0 = y_centre, group = profile)) +
    geom_circle(aes(r = minir, col= label, fill  = label)) + 
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
  
  gq <- ggplot(pigg_plotq %>% filter(perc > 5), aes(x0=x_centre, y0 = y_centre, group = profile)) +
    geom_circle(aes(r = minir, col= label, fill  = label)) + 
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
  g <- gp + gq + plot_layout(guides = 'collect') +
    plot_layout(widths = c(1, 1))
  
  if(save){
    ggsave(paste0("plots/",plot_name,".pdf"), width = 45, height = 25)
  }
  
  plot(g)
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
    geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.6) +
    geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name)),size = 1) +
    geom_point(aes(col = name),size = 1.5) +
    facet_wrap(name~parent_strain, ncol = 2) + 
    scale_color_discrete("MGE") + 
    scale_linetype_discrete("Piglet") + 
    scale_y_continuous("Proportion of population", lim = c(0,1.1)) + 
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
