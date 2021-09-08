Rprof()

p11 <- gain2(P,P, mu_in)

Rprof(NULL)  

summaryRprof()


system.time(p1new <- plyr::count(rph, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")) %>%
              mutate(Total = select(., v1:v10) %>% rowSums()) %>%  # count how many
              filter(Total>0) %>% select(-Total) )

system.time(rph2 <- rph[-which(rowSums(rph)==0),])
system.time(p1new <- plyr::count(rph2, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")) %>%
              mutate(Total = select(., v1:v10) %>% rowSums()) %>%  # count how many
              filter(Total>0) %>% select(-Total) )

system.time(rph2 <- rph %>% remove_empty("rows")) # remove no transfer rows
system.time(p1new <- plyr::count(rph2, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")) %>%
  mutate(Total = select(., v1:v10) %>% rowSums()))

1.408
0.578 + 0.823
0.549 + 1.428


system.time(p1new <- plyr::count(rph, vars = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")))
system.time(p1new$Total <- select(p1new, v1:v10) %>% rowSums())
system.time(p1new <- p1new %>% filter(Total>0) %>% select(-Total))


system.time(for(ijh in 1:20){rph %>% remove_empty("rows")})
system.time(for(ijh in 1:20){rph[-which(rowSums(rph)==0),]})

system.time(for(ijh in 1:10){p11 <- gain(P,Q, mu_in)})
system.time(for(ijh in 1:10){p11 <- gain2(P,Q, mu_in)})
system.time(for(ijh in 1:10){p11 <- gain_mge(P,Q, mu_in)})
