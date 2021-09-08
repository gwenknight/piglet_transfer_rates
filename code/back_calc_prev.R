##### Back calc prevalence

### Use combinatorics to work out the distribution of expected bugs given what observed

# needed libraries
library(ggridges) # for ridge plot
library(tidyverse)
theme_set(theme_ridges())
theme_set(theme_bw(base_size = 26))
library(zoo) # for interpolation

############################### If piglet has a prevalence of x, what is the probability of seeing 1/2/3...20 colonies on the plate? ####################################################################################################################

# A = profile name
# x = prevalence in pig
# c = prevalence on plate
# p = probability that prevalence in pig is x given c

# e.g. 
# c = 1/20 for A 
# what is the p(A at x% in pig AND 1/20 = A AND 19/20 not A)? 

pA = c()
for(j in 0:20){
  for(i in 0:100){
    
    # Given that A = x%, what's the probability that 1/20 = A and 19/20 not A? 
    pA = rbind(pA, c(i, j, choose(20, j) * ((i/100)^j) * ((100-i)/100)^(20-j)))
    
  }
}
pA <- as.data.frame(pA)
colnames(pA) <- c("prev","n_colonies","prob")
ggplot(pA %>% filter(prev == 50), aes(x=n_colonies, y = prob)) + geom_line()
ggplot(pA %>% filter(prev == 0), aes(x=n_colonies, y = prob)) + geom_line() # near edges higher chance get 0 than 1 but 1/20 = 0.05 so makes sense
ggplot(pA %>% filter(prev == 1), aes(x=n_colonies, y = prob)) + geom_line() # 80% chance 0 colonies if prevalence in pig is 1, 16% chance 1 colony, 1.58% chance 2 colonies...
# prev = Percentage of colonies with this profile in the pigs
# n_colonies = number of colonies with this profile out of 20
# prob = probability that n_colonies = n_colonies if prev = prev

pA %>% filter(prev == 1) %>% summarise(sum(prob)) # sum over a single prevalence = 1. 


############################### If have 1 colony, what is the prevalence in pig of this profile?  ####################################################################################################################

## ????? HARD especially at 0 colonies on plate => could be some in pig. Can't use the above method... 
#### WRONG BELOW
# # A = profile name
# # x = prevalence in pig
# # c = prevalence on plate
# # p = probability that prevalence in pig is x given c
# 
# # e.g. 
# # c = 1/20 for A 
# # what is the p(1/20 = A on plate AND A at x% in pig)? 
# 
# pA = c()
# for(j in 0:20){
#   for(i in 0:100){
#     
#     # Given that 1/20 = A,  what's the probability that x% = A and 1-x% =  not A? 
#     pA = rbind(pA, c(i, j, choose(100, i) * ((j/20)^i) * ((20-j)/20)^(100-i)))
#     
#   }
# }
# pA <- as.data.frame(pA)
# colnames(pA) <- c("prev","n_colonies","prob")
# ggplot(pA %>% filter(prev == 50), aes(x=n_colonies, y = prob)) + geom_line()
# ggplot(pA %>% filter(n_colonies == 0), aes(x=n_colonies, y = prob)) + geom_line() # near edges higher chance get 0 than 1 but 1/20 = 0.05 so makes sense
# ggplot(pA %>% filter(prev == 1), aes(x=n_colonies, y = prob)) + geom_line() # 80% chance 0 colonies if prevalence in pig is 1, 16% chance 1 colony, 1.58% chance 2 colonies...
# # prev = Percentage of colonies with this profile in the pigs
# # n_colonies = number of colonies with this profile out of 20
# # prob = probability that n_colonies = n_colonies if prev = prev


################################# Explore above data: if piglet has a prevalence of x, what is the probability of seeing 1/2/3...20 colonies on the plate? ################################################
#### Look at prevalence and see what prob this n colonies with this prev
pA <- pA %>% group_by(prev) %>% mutate(prob_this_n_colonies = prob/sum(prob)) %>% ungroup()


ggplot(pA %>% filter(prev %in% c(0,seq(0,100, by = 10))), aes(x = n_colonies, y = prob_this_n_colonies, group = prev)) + geom_line(aes(col=factor(prev))) + 
#  facet_wrap(~prev, scales = "free") +
  geom_vline(aes(xintercept = (prev/100)*20, col = factor(prev))) + 
  scale_x_continuous("Number of colonies with this profile") + 
  scale_color_discrete("Prevalence in pig") + 
  scale_y_continuous("Probability that n_colonies = n_colonies if prev = prev") + 
  ggtitle("If prev in pig is colour, prev on plate is x")
ggsave("plots/Probability_prev_col.pdf", width = 8)

##### Says:
#### if prevalence is x, what is the most likely number of positive colonies on a plate? = LIKELIHOOD
#### Probability get this data given this model output

# WRONG WAY
pA %>% filter(n_colonies == 5) %>% summarise(sum(prob)) # > 1
# ggplot(pA %>% select(prev, n_colonies, prob), 
#        aes(x = prev, y = n_colonies, group = n_colonies,height = prob)) + #geom_line(aes(group = n_colonies))
#   geom_ridgeline(aes(fill = n_colonies))  + 
#   scale_y_continuous("Number of positive colonies out of 20") + 
#   scale_x_continuous("Prevalence in piglet")
# ggsave("plots/ridge_plot.pdf")

#### PLOT correct
ggplot(pA %>% select(prev, n_colonies, prob) %>% filter(prev %in% seq(0,100,by = 5)), 
       aes(x = n_colonies, y = prev, group =prev,height = prob)) + #geom_line(aes(group = n_colonies))
  geom_ridgeline(aes(fill = prev))  + 
  scale_x_continuous("Number of positive colonies out of 20") + 
  scale_y_continuous("Prevalence in piglet")
ggsave("plots/ridge_plot.pdf")


ggplot(pA %>% select(prev, n_colonies, prob) %>% filter(prev %in% seq(0,100,by = 5), prob > 0.01), 
       aes(x = n_colonies, y = prev, group =prev,height = prob)) + #geom_line(aes(group = n_colonies))
  geom_ridgeline(aes(fill = prev))  + 
  scale_x_continuous("Number of positive colonies out of 20") + 
  scale_y_continuous("Prevalence in piglet") 
ggsave("plots/ridge_plot_0.01.pdf")

ggplot(pA %>% select(prev, n_colonies, prob) %>% filter(prev %in% seq(0,100,by = 5), prob > 0.01), 
       aes(x = n_colonies, y = prev, group =prev,height = prob)) + #geom_line(aes(group = n_colonies))
  geom_ridgeline(aes(fill = prev))  + 
  scale_x_continuous("Number of positive colonies out of 20") + 
  scale_y_continuous("Prevalence in piglet") + 
  geom_text(x=5, y=90, label= "e.g. prevalence in pig 35%\nprobability see data (number of colonies on plate)\ngiven by distribution") +  # NO! "e.g. 1/20 = this profile (y axis), \nv unlikely prevalence in piglet of 25% (x-axis)\nbut likely that 3-8% do (dashed lines)") 
  geom_hline(yintercept = c(33,38), linetype = "dashed")
ggsave("plots/ridge_plot_0.01.pdf")

### Save
pA$prev <- pA$prev / 100
pA$n_colonies_prev <- pA$n_colonies/20
pA <- pA %>% select(prev, prob, n_colonies_prev)

ggplot(pA, aes(x = n_colonies_prev, y = prev, group = prev,height = prob)) + #geom_line(aes(group = n_colonies))
  geom_ridgeline(aes(fill = prev)) ## WEIRD!!

### Interpolate: for each prevalence, need to fill in probability for each n_colonies_prev
pA_interp <- pA %>%
  group_by(prev) %>%
  arrange(prev, n_colonies_prev) %>%
  complete(n_colonies_prev = seq(0, 1, len = 10001)) %>%
  mutate(prob_all_t = na.approx(prob)) %>% # Need to divide by new total - no longer just a bar at say n_colonies = 0.1, but now have a bit of that at 0.09, 0.95 etc
  mutate(total = sum(prob_all_t), 
         prob_all = prob_all_t/total) %>% 
  select(-c(prob, prob_all_t)) # remove non interpolated

## check prob = 1
pA_interp %>% filter(prev == 0.53) %>% summarise(sum(prob_all))

write.csv(pA_interp, "seen_predicted.csv")






####### BELOW WRONG WAY TO THINK ABOUT IT

## if find 1 positive colony then probably X wrong way to think of it: if prev is something, this is most likely number of colonies on the plate with this profile
pA %>% filter(n_colonies == 1) %>% mutate(max = max(prob), equalmax = ifelse(prob == max, 1, 0)) %>% filter(equalmax == 1)
# a prevalence of 5%, makes sense as 1/20 = 5%

# ## if find 0 positive colony then probably
# pA %>% filter(n_colonies == 0) %>% mutate(max = max(prob), equalmax = ifelse(prob == max, 1, 0)) %>% filter(equalmax == 1)
# # a prevalence of 0
# pA %>% filter(n_colonies == 0) # but some variation -> 
# 
# ## if find 6  positive colonies then probably
# pA %>% filter(n_colonies == 6) %>% mutate(max = max(prob), equalmax = ifelse(prob == max, 1, 0)) %>% filter(equalmax == 1)
# # a prevalence of 30






#### OLD

# ##### (1) ######
# ##### Now flip: instead say, if see this many colonies - what is the probability of that prevalence?WRONG
# seen_predicted = pA %>% group_by(n_colonies) %>% mutate(max = max(prob), equalmax = ifelse(prob == max, 1, 0), n_cols_prev = 100*n_colonies / 20) 
# 
# # Check maximum probability matches prevalence of underlying - this is the most likely thing to count
# ggplot(seen_predicted %>% 
#          filter(equalmax == 1) %>% dplyr::select(n_colonies, n_cols_prev, prev)
#        , aes(x=n_cols_prev, y = prev)) + geom_point() + 
#   geom_line(data = as.data.frame(cbind(x = c(0:100), y = c(0:100))), aes(x=x, y = y)) + 
#   scale_x_continuous("Prevalence predicted by proportion of 20 that are this profile") + 
#   scale_y_continuous("Prevalence in pig population") + 
#   ggtitle("Filtered on max probability") + 
#   geom_text(x=30, y=90, label="e.g. 4/20 = 25% of colonies have this profile (x-axis):\nimplies 25% prevalence in pig (y axis)")
# 
# 
# 
# # Not so good as stages on y axis >> 1 so can't see the distribution so well
# ggplot(seen_predicted %>% select(prev, n_colonies, prob)%>% filter(prob > 0.01), 
#        aes(x = prev, y = 100*n_colonies/20, group = n_colonies,height = prob)) + #geom_line(aes(group = n_colonies))
#   geom_ridgeline(aes(fill = n_colonies))  + 
#   scale_y_continuous("Prevalence in sample") + 
#   scale_x_continuous("Prevalence in piglet")
# ggsave("plots/ridge_plot_0.01_asperc.pdf")
# 
# ## For single number of colonies out of 20
# ncol = 16 # how many with this profile? 
# seen_predicted %>% filter(n_colonies == ncol) %>% ggplot(aes(x=prev, y = prob_this_n_colonies)) + 
#   geom_line() + geom_vline(xintercept = 100 *ncol / 20, linetype = "dashed")  + 
#   scale_y_continuous("Probablility this many colonies with the profile, given prevalence in piglet") + 
#   scale_x_continuous("Prevalence in piglet")
# ggsave("plots/16col_pos_prob.pdf")
# 
# write.csv(seen_predicted, "seen_predicted.csv")
# 
# ###### (2) 
# # Need to change around: if prevalence in pig is x what is the chance that see y on plate? 
# 
# #p (A = x% AND 1/20 = A AND 19/20 not A)
# pA = c()
# 
# for(i in 0:20){
#   for(j in 0:100){
#     
#     # Given that A = x%, what's the probability that 1/20 = A and 19/20 not A? 
#     # Given that 1/20 = A on plate, what's the probability that x% = A and 1-x not A?
#     # Have to calculate combinatorics??!? 
#     pA = rbind(pA, c(i, j, choose(100, j) * ((i/20)^j) * ((20-i)/20)^(100-j)))
#     
#   }
# }
# pA <- as.data.frame(pA)
# colnames(pA) <- c("n_colonies","prev","prob")
# ggplot(pA %>% filter(prev == 20), aes(x=n_colonies, y = prob)) + geom_line() # more spiky than other pA as only limited number can be positive out of 20 colonies
# # prev = Percentage of colonies with this profile in the pigs
# # n_colonies = number of colonies with this profile out of 20
# # prob = probability that prev = prev if n_colonies = n_colonies
# 
# ### Note now that
# pA %>% filter(prev == 2) %>% summarise(sum(prob)) # < 1
# pA %>% filter(n_colonies == 2) %>% summarise(sum(prob)) # this is 1/ If n colonies is 2, distributino of prevalence
# ggplot(pA %>% filter(n_colonies == 2), aes(x=prev, y = prob)) + geom_line() # this distribution 
# 
# theme_set(theme_bw(base_size = 26))
# ggplot(pA %>% filter(prev %in% c(0,seq(0,100, by = 5))), aes(x = n_colonies, y = prob_this_prevalence, group = prev)) + geom_line(aes(col=factor(prev))) + 
#   #  facet_wrap(~prev, scales = "free") +
#   geom_vline(aes(xintercept = (prev/100)*20, col = factor(prev))) + 
#   scale_x_continuous("Number of colonies with this profile") + 
#   scale_color_discrete("Prevalence in pig") + 
#   scale_y_continuous("Probability that  prev = prev if n_colonies = n_colonies") + 
#   ggtitle("if number of colonies with this profile is x, what is the most likely prevalence in the pig?")
# ggsave("plots/Probability_prev_col.pdf", width = 8)

