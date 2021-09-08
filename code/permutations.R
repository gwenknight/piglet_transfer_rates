### First pass at basic system to generate range of profiles to match Figure 2 of McCarthy et al
#### How to run basic code and 

library(tidyverse)
library(janitor)
library(ggforce)
library(doParallel)
library(RColorBrewer)
library(patchwork)
library(lhs)
theme_set(theme_bw(base_size = 8))

# numCores <- detectCores() - 2
# numCores
# registerDoParallel(numCores)
# stopCluster(cl)

## Functions
source("code/functions.R")

### All 
all_comb <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
mge <- c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

## Rate of transfer: 10 elements
mu = rep(0.0000000004,10)
## Rate of loss: 10 elements
gamma = rep(0.00000001,10)
## Rate of growth: cost per element if present
growth = rep(0.001,10)
grate = 0.6 # growth rate
K = 1 * 10^7 # carrying capacity

## Time step = 1 hr
tsteps = 16 * 24


##works: 
mu = rep(0.0000004,10) 
gamma = rep(0.00001,10) 
growth = rep(0.001,10)
grate = 0.1

parameters <- c(mu2 = 0.000001,mu5 = 0.000001,mu6 = 0.000001,
                mu7 = 0.000001,mu8 = 0.000001,mu10 = 0.000001,
                gamma2 = 0.000001,gamma5 = 0.000001,gamma6 = 0.000001,
                gamma7 = 0.000001,gamma8 = 0.000001,gamma10 = 0.000001,
                f2 = 0.01,f5 = 0.01,f6 = 0.01,
                f7 = 0.01,f8 = 0.01,f10 = 0.01,
                grow = 0.01)

#out <- run_sim(384, c(mu, gamma, growth, grate))
#out <- run_sim(50, c(1000*mu, 1000*gamma, growth, 0.1))
out <- run_sim(384, parameters)
P_all <- out$P_all
Q_all <- out$Q_all

P <- P_all %>% filter(time == 384)
pprop <- 100 * P$freq / sum(P$freq)
length(which(pprop > 5))

Q <- Q_all %>% filter(time == 384)
qprop <- 100 * Q$freq / sum(Q$freq)
length(which(qprop > 5))


#Save
write.csv(P_all,"first_P_all.csv")
write.csv(Q_all,"first_Q_all.csv")



### How summarise results? 

### Totals over time (384 hrs = 16 days)
P_all %>% group_by(time) %>% summarise(totalP = sum(freq), .groups = 'drop') %>% ggplot(aes(x=time, y = totalP)) + geom_line() + 
  geom_line(data = Q_all %>% group_by(time) %>% summarise(totalQ = sum(freq), .groups = 'drop'),aes(x=time, y = totalQ),col = "red") + 
  geom_vline(xintercept = 96) + geom_hline(yintercept = K/2) + 
  geom_vline(xintercept = c(72, 210), linetype = "dashed")  + geom_hline(yintercept = K/2) +
  scale_y_continuous("Number of bugs")
ggsave("plots/first_over_time.pdf")

### How many profiles? 
P_sample = P_all %>% filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% arrange(time, -freq) %>% mutate(total = sum(freq), prop = freq / total, value = 1, profile = cumsum(value)) %>% select(-value)
Q_sample = Q_all %>% filter(time %in% c(4,48,96,288,384)) %>% group_by(time) %>% arrange(time, -freq) %>% mutate(total = sum(freq), prop = freq / total, value = 1, profile = cumsum(value)) %>% select(-value)

P_sample$parent = 0
Q_sample$parent = 0
for(i in 1:dim(P_sample)[1]){
  if(all.equal(as.numeric(P_sample[i,1:10]), as.numeric(Pinit[,1:10])) == "TRUE"){P_sample[i,"parent"] = 1}
}
for(i in 1:dim(Q_sample)[1]){
  if(all.equal(as.numeric(Q_sample[i,1:10]), as.numeric(Qinit[,1:10])) == "TRUE"){Q_sample[i,"parent"] = 1}
}

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

colnames(P_sample) <- c(c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"),colnames(P_sample)[11:19])
colnames(Q_sample) <- c(c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"),colnames(Q_sample)[11:19])

P_plot <- as.data.frame(P_sample %>% pivot_longer(cols = v1:v10)) %>% mutate(label = ifelse(value == 1, name,0)) %>% 
  mutate(mge = as.numeric(substr(name,2,3)),
         plasmid = ifelse(name %in% c("v1","v2","v3","v4","v5","v9"),0,1),
         x_centre = ifelse(plasmid == 0, bigr*cos(pt[mge]) + x_centrebig,pt[mge]),
         y_centre = ifelse(plasmid == 0, bigr*sin(pt[mge]) + x_centrebig,y_centre_plasmid))

Q_plot <- as.data.frame(Q_sample %>% pivot_longer(cols = v1:v10)) %>% mutate(label = ifelse(value == 1, name,0)) %>% 
  mutate(mge = as.numeric(substr(name,2,3)),
         plasmid = ifelse(name %in% c("v1","v2","v3","v4","v5","v9"),0,1),
         x_centre = ifelse(plasmid == 0, bigr*cos(pt[mge]) + x_centrebig,pt[mge]),
         y_centre = ifelse(plasmid == 0, bigr*sin(pt[mge]) + x_centrebig,y_centre_plasmid))

mm <- max(P_plot$profile,Q_plot$profile)
P_plot$profile = factor(P_plot$profile, levels=seq(mm,1,-1)) # reorder profiles to match McCarthy Fig2

# how add a rectangle for parent???
#d=data.frame(x1=rep(4.2,5), x2=rep(6.7,5), y1=rep(6.2,5), y2=rep(6.2,5), label = "parent", profile = 1, time= c(4,48,96,288,384)) 
#geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = label), alpha=0.5) + 

P_plot$perc <- paste0(100 * signif(P_plot$prop,2),"%")
Q_plot$perc <- paste0(100 * signif(Q_plot$prop,2),"%")
gp <- ggplot(P_plot, aes(x0=x_centre, y0 = y_centre, group = profile)) + geom_circle(aes(r = minir, col= label, fill  = label)) + 
  geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
  geom_text(aes(x = 5, y = 6.5, label = perc)) +
  facet_grid(time ~ profile) + 
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

gq <- ggplot(Q_plot, aes()) + geom_circle(aes(x0=x_centre, y0 = y_centre,r = minir, col= label, fill  = label, group = profile)) + 
  geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
  geom_text(aes(x = 5, y = 6.5, label = perc)) + 
  facet_grid(time ~ profile) + 
  scale_fill_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","parent"), 
                    values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8",'slategray2'),drop = FALSE) + 
  scale_color_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","parent"), 
                     values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8",'slategray2'),drop = FALSE) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

gp + gq + plot_layout(guides = 'collect') + 
  plot_layout(widths = c(2, 1))
ggsave("plots/first_pass_pig.pdf", width = 45, height = 8)

### filter > 5%
gp2 <- ggplot(P_plot %>% filter(prop > 0.05), aes(x0=x_centre, y0 = y_centre, group = profile)) + geom_circle(aes(r = minir, col= label, fill  = label)) + 
  geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
  geom_text(aes(x = 5, y = 6.5, label = perc)) +
  facet_grid(time ~ profile) + 
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

gq2 <- ggplot(Q_plot %>% filter(prop > 0.05), aes()) + geom_circle(aes(x0=x_centre, y0 = y_centre,r = minir, col= label, fill  = label, group = profile)) + 
  geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
  geom_text(aes(x = 5, y = 6.5, label = perc)) + 
  facet_grid(time ~ profile) + 
  scale_fill_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","parent"), 
                    values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8",'slategray2'),drop = FALSE) + 
  scale_color_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","parent"), 
                     values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8",'slategray2'),drop = FALSE) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

gp2 + gq2 + plot_layout(guides = 'collect') + 
  plot_layout(widths = c(1, 1))
ggsave("plots/first_pass_pig_filter.pdf")

gp2 <- ggplot(P_plot %>% filter(prop > 0.01), aes(x0=x_centre, y0 = y_centre, group = profile)) + geom_circle(aes(r = minir, col= label, fill  = label)) + 
  geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
  geom_text(aes(x = 5, y = 6.5, label = perc)) +
  facet_grid(time ~ profile) + 
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

gq2 <- ggplot(Q_plot %>% filter(prop > 0.01), aes()) + geom_circle(aes(x0=x_centre, y0 = y_centre,r = minir, col= label, fill  = label, group = profile)) + 
  geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
  geom_text(aes(x = 5, y = 6.5, label = perc)) + 
  facet_grid(time ~ profile) + 
  scale_fill_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","parent"), 
                    values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8",'slategray2'),drop = FALSE) + 
  scale_color_manual("MGE",breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","parent"), 
                     values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8",'slategray2'),drop = FALSE) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

gp2 + gq2 + plot_layout(guides = 'collect') + 
  plot_layout(widths = c(1, 1))
ggsave("plots/first_pass_pig_filter1.pdf")




# P_plot$start <- 1:10
# P_plot$end <- 2:11
# 
# ggplot(P_plot, aes(xmin = start, xmax = end, ymin = 4, ymax = 5, fill = label, )) + 
#   geom_rect(aes(color = factor(parent))) + scale_y_continuous(lim = c(0,5)) + 
#   coord_polar() + facet_grid(time + plasmid ~ profile) + 
#   scale_fill_manual(breaks= c("0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"), 
#                     values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#984EA3", "#FF7F00","#E41A1C","#377EB8")) + 
#   scale_color_manual("Parent profile",breaks = c(0,1), values = c("grey","black"))
# 
# 
# 
# P_plot_notplasmid <- as.data.frame(P_sample %>% pivot_longer(cols = v1:v10)) %>% mutate(label = ifelse(value == 1, name,0)) %>% mutate(plasmid = ifelse(name %in% c("v1","v2","v3","v4","v5","v6"),0,1))
# P_plot_notplasmid$start <- 1:6
# P_plot_notplasmid$end <- 2:7
# 
# P_plot_plasmid <- as.data.frame(P_sample %>% pivot_longer(cols = v1:v10)) %>% mutate(label = ifelse(value == 1, name,0)) %>% filter(name %in% c("v7","v8","v9","v10"))
# P_plot_plasmid$start <- 1:4
# P_plot_plasmid$end <- 2:5
# 
# 
# g1 <- ggplot(P_plot_notplasmid, aes(xmin = start, xmax = end, ymin = 4, ymax = 5, fill = label)) + geom_rect() + scale_y_continuous(lim = c(0,5)) + 
#   coord_polar() + facet_grid(time ~ profile) + 
#   scale_fill_manual(values = c("white","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"))# = "Set1"[1:7],direction=-1)
# 
# g2 <- ggplot(P_plot_plasmid, aes(xmin = start, xmax = end, ymin = 4, ymax = 5, fill = label)) + geom_rect() + scale_y_continuous(lim = c(0,5)) + 
#   facet_grid(time ~ profile) + 
#   scale_fill_manual(breaks= c("0","v7","v8","v9","v10"), values = c("white","#984EA3", "#FF7F00","#E41A1C","#377EB8"))# = "Set1"[1:7],direction=-1)
# 
# g1/g2