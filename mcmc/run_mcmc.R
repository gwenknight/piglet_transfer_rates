
# script to run mcmc

dist_like <- read.csv("data/seen_predicted.csv")[,-1]
data <- read.csv("data/data_to_fit.csv")[,-1]
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
totals <- read.csv("data/totals_bug.csv")[,-1] %>%
  select(time,value, parent, lim) %>%
  pivot_wider(names_from = lim) %>%
  mutate(weight = 1/(max - min))
