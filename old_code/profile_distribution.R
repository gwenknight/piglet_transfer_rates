#### distribution of profiles at end

#piglet, parent, profile, percentage
prof_dist_data_end <- as.data.frame(rbind(c(1,1,1,100),
      c(1,2,1,45),
      c(1,2,2,25),
      c(1,2,3,10),
      c(1,2,4,5),
      c(1,2,5,15),
      c(2,1,1,95),
      c(2,1,2,5),
      c(2,2,1,40),
      c(2,2,6,5),
      c(2,2,2,25),
      c(2,2,7,20),
      c(2,2,8,10),
      c(3,1,2,5),
      c(3,1,1,95),
      c(3,2,2,20),
      c(3,2,4,5),
      c(3,2,1,45),
      c(3,2,7,5),
      c(3,2,5,5),
      c(4,1,1,100),
      c(4,2,2,70),
      c(4,2,1,30)))
colnames(prof_dist_data_end) <- c("piglet", "parent", "profile", "percentage")

ggplot(prof_dist_data_end, aes(x=percentage, group_by = interaction(parent, profile))) + 
  geom_histogram() + 
  facet_grid(parent~profile)

### In parent 1 there is really only one profile.. 
prof_dist_data_end %>% filter(parent == 1, profile == 1)
mean(c(100,95,95,100))
sd(c(100,95,95,100))

## Uniform over 0.9 - 1 to give 20% leeway? too harsh... set to normal 
plot(seq(0,1,0.01),log(dunif(seq(0,1,0.01), min = 0.9, max = 1)))
plot(seq(0,1,0.01),(dnorm(seq(0,1,0.01), mean = 0.975, sd = 0.0288)))
plot(seq(0,1,0.01),log(dnorm(seq(0,1,0.01), mean = 0.975, sd = 0.0288)))


### These two profiles in parent 2, give over 65% of the profiles - fit to this instead? 
par2 <- prof_dist_data_end %>% filter(parent == 2, profile < 3) %>% 
  group_by(parent, piglet) %>% summarise(s = sum(percentage))

mean(par2$s)
sd(par2$s)

plot(seq(0,100,1),log(dnorm(seq(0,100,1), mean = 75, sd = 17)))
plot(seq(0,100,1),(dtruncnorm(seq(0,100,1), a = 0, b = 100, mean = 75, sd = 17)))
plot(seq(0,100,1),log(dtruncnorm(seq(0,100,1), a = 0, b = 100, mean = 75, sd = 17)))
## Do a truncated normal around the mean + sd from data. Up to 100 = 1 = max proportion. Symmetric

### Check likelihood contribution 
log(dunif(0.91, 0.9,1)) + log(dtruncnorm(0.5,a = 0.5, b = 1, mean = 0.75, sd = 0.17))



