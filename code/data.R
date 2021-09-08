### Data

#### (1) Read in data from the main bubble figure and plot
#### (2) Read in data on total number of bugs and plot

# Load needed libraries 
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggforce)
library(ggridges)
library(mgcv) # for uniquecombs function
theme_set(theme_bw(base_size = 11))

# Loads in basic functions for this model (need getLabels)
source("code/functions.R")

######################################################################################### Main data ####################################################################################################################################################################################
## initial conditions
Pinit = t(c(1,1,1,1,1,1,1,1,0,0))
Qinit = t(c(0,0,0,0,0,0,0,0,1,1))
Pinit <- as.data.frame(Pinit)
colnames(Pinit) = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")
Qinit <- as.data.frame(Qinit)
colnames(Qinit) = c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

mge <- c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

### Total numbers
t4 = 10^2
t48 = 5 * 10^5
t72 = 10^7
t288 = 10^7 
t384 = 10^7

#      pig parent v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 n time

p1data <- as.data.frame(rbind(c(1,1,1,1,1,1,1,1,1,1,0,0,10000,0),
                              c(1,1,1,1,1,1,1,1,1,1,0,0,0.7*t4,4),
                              c(1,1,1,1,1,1,1,0,1,1,0,0,0.1*t4,4),
                              c(1,1,1,1,1,1,1,1,1,0,0,0,0.15*t4,4),
                              c(1,1,1,1,1,1,1,0,0,0,0,0,0.05*t4,4),
                              c(1,1,1,1,1,1,1,1,1,1,0,0,0.95*t48,48),
                              c(1,1,1,1,1,1,1,1,0,1,0,0,0.05*t48,48),
                              c(1,1,1,1,1,1,1,1,1,1,0,0,0.85*t72,72),
                              c(1,1,1,1,1,1,1,0,1,1,0,0,0.15*t72,72),
                              c(1,1,1,1,1,1,1,1,1,1,0,0,1*t288,288),
                              c(1,1,1,1,1,1,1,1,1,1,0,0,1*t384,384),
                              
                              c(1,2,0,0,0,0,0,0,0,0,1,1,10000,0),
                              c(1,2,0,0,0,0,0,0,0,0,1,1,0.25*t4,4),
                              c(1,2,0,0,0,0,0,0,0,0,1,0,0.7*t4,4),
                              c(1,2,0,0,0,0,0,0,1,0,1,1,0.05*t4,4),
                              c(1,2,0,0,0,0,0,0,0,0,1,1,0.05*t48,48),
                              c(1,2,0,1,0,0,1,0,1,0,1,0,0.25*t48,48),
                              c(1,2,0,0,0,0,1,0,1,0,1,0,0.05*t48,48),
                              c(1,2,0,1,0,0,0,0,0,0,1,0,0.05*t48,48),
                              c(1,2,0,1,0,0,0,1,1,0,1,0,0.05*t48,48),
                              c(1,2,0,1,0,0,1,0,1,0,1,1,0.20*t48,48),
                              c(1,2,0,0,0,0,1,0,0,0,1,1,0.05*t48,48),
                              c(1,2,0,0,0,0,0,0,1,0,1,1,0.20*t48,48),
                              c(1,2,0,1,0,0,1,0,0,0,1,0,0.10*t48,48),
                              c(1,2,0,0,0,0,1,0,0,0,1,0,0.05*t72,72),
                              c(1,2,0,1,0,0,1,0,1,0,1,0,0.40*t72,72),
                              c(1,2,0,0,0,0,1,0,1,0,1,0,0.10*t72,72),
                              c(1,2,0,0,0,0,0,0,0,0,1,0,0.20*t72,72),
                              c(1,2,0,1,0,0,0,0,1,0,1,0,0.05*t72,72),
                              c(1,2,0,0,0,0,0,0,1,0,1,0,0.20*t72,72),
                              c(1,2,0,0,0,0,1,0,0,0,1,0,0.30*t288,288),
                              c(1,2,0,1,0,0,1,0,0,0,1,0,0.45*t288,288),
                              c(1,2,0,1,0,0,1,1,0,0,1,0,0.05*t288,288),
                              c(1,2,0,0,0,0,0,0,0,0,1,0,0.10*t288,288),
                              c(1,2,0,1,0,0,0,0,0,0,1,0,0.10*t288,288),
                              c(1,2,0,1,0,0,1,1,0,0,1,0,0.10*t384,384),
                              c(1,2,0,1,0,0,1,0,0,0,1,0,0.10*t384,384),
                              c(1,2,0,1,0,0,0,1,0,0,1,0,0.10*t384,384),
                              c(1,2,0,0,0,0,1,1,0,0,1,0,0.10*t384,384),
                              c(1,2,0,1,0,0,0,0,0,0,1,0,0.10*t384,384)))


p2data <-  as.data.frame(rbind(c(2,1,1,1,1,1,1,1,1,1,0,0,10000,0),
                               c(2,1,1,1,1,1,1,1,1,1,0,0,0.75*t4,4),
                               c(2,1,1,1,1,1,1,1,1,0,0,0,0.05*t4,4),
                               c(2,1,1,1,1,1,1,0,0,0,0,0,0.30*t4,4),
                               c(2,1,1,1,1,1,1,1,1,1,0,0,t48,48),
                               c(2,1,1,1,1,1,1,1,1,1,0,0,t72,72),
                               c(2,1,1,1,1,1,1,1,1,1,0,0,0.2*t288,288),
                               c(2,1,1,1,1,1,1,0,1,1,0,0,0.8*t288,288),
                               c(2,1,1,1,1,1,1,1,1,1,0,0,0.95*t384,384),
                               c(2,1,1,1,1,1,1,1,1,0,0,0,0.05*t384,384),
                               
                               c(2,2,0,0,0,0,0,0,0,0,1,1,10000,0),
                               c(2,2,0,0,0,0,0,0,0,0,1,1,0.3*t4,4),
                               c(2,2,0,0,0,0,0,0,0,0,1,0,0.7*t4,4),
                               c(2,2,0,1,0,0,0,0,0,0,1,0,0.65*t48,48),
                               c(2,2,0,0,0,0,0,0,0,0,1,0,0.35*t48,48),
                               c(2,2,0,1,0,0,0,0,0,0,1,0,0.35*t72,72),
                               c(2,2,0,0,0,0,0,0,0,0,1,0,0.05*t72,72),
                               c(2,2,0,1,0,0,1,0,0,0,1,0,0.50*t72,72),
                               c(2,2,0,0,0,0,1,0,0,0,1,0,0.10*t72,72),
                               c(2,2,0,1,0,0,1,0,0,0,1,0,0.80*t288,288),
                               c(2,2,0,0,0,0,1,0,0,0,1,0,0.20*t288,288),
                               
                               c(2,2,0,1,0,0,1,0,0,0,1,0,0.25*t384,384),
                               c(2,2,0,0,0,0,1,0,0,0,1,0,0.20*t384,384),
                               c(2,2,0,1,0,0,1,0,1,0,1,0,0.40*t384,384),
                               c(2,2,0,1,0,0,1,1,1,0,1,0,0.05*t384,384),
                               c(2,2,0,0,0,0,0,1,0,0,1,0,0.10*t384,384)))

p3data <-  as.data.frame(rbind(c(3,1,1,1,1,1,1,1,1,1,0,0,10000,0),
                               c(3,1,1,1,1,1,1,1,1,1,0,0,0.95*t4,4),
                               c(3,1,1,1,1,1,1,0,1,1,0,0,0.05*t4,4),
                               c(3,1,1,1,1,1,1,1,1,1,0,0,t48,48),
                               c(3,1,1,1,1,1,1,1,1,1,0,0,t72,72),
                               c(3,1,1,1,1,1,1,1,1,1,0,0,t288,288),
                               c(3,1,1,1,1,1,1,1,1,1,0,0,0.95*t384,384),
                               c(3,1,1,1,1,1,1,1,1,0,0,0,0.05*t384,384),
                               
                               c(3,2,0,0,0,0,0,0,0,0,1,1,10000,0),
                               c(3,2,0,0,0,0,0,0,0,0,1,1,t4,4),
                               c(3,2,0,1,0,0,0,0,0,0,1,0,0.35*t48,48),
                               c(3,2,0,0,0,0,0,0,0,0,1,0,0.6*t48,48),
                               c(3,2,0,1,0,0,0,1,0,1,1,0,0.05*t48,48),
                               c(3,2,0,0,0,0,1,0,0,0,1,0,0.25*t72,72),
                               c(3,2,0,0,0,0,0,0,0,0,1,0,0.15*t72,72),
                               c(3,2,0,0,0,0,1,0,1,0,1,0,0.35*t72,72),
                               c(3,2,0,0,0,0,0,0,1,0,1,0,0.20*t72,72),
                               c(3,2,0,0,0,0,1,0,1,1,1,0,0.05*t72,72),
                               
                               c(3,2,0,1,0,0,1,0,0,0,1,0,0.30*t288,288),
                               c(3,2,0,1,0,0,1,0,0,1,1,0,0.20*t288,288),
                               c(3,2,0,1,0,0,1,0,1,0,1,0,0.10*t288,288),
                               c(3,2,0,1,0,0,1,0,1,1,1,0,0.35*t288,288),
                               c(3,2,0,1,0,0,0,0,1,1,1,0,0.05*t288,288),
                               
                               c(3,2,0,1,0,0,1,0,0,0,1,0,0.20*t384,384),
                               c(3,2,0,0,0,0,1,0,1,0,1,0,0.05*t384,384),
                               c(3,2,0,1,0,0,1,0,1,0,1,0,0.65*t384,384),
                               c(3,2,0,0,0,0,1,0,0,0,1,0,0.05*t384,384),
                               c(3,2,0,1,0,0,0,0,0,0,1,0,0.05*t384,384)))

p4data <-  as.data.frame(rbind(c(4,1,1,1,1,1,1,1,1,1,0,0,10000,0),
                               c(4,1,1,1,1,1,1,1,1,1,0,0,t4,4),
                               c(4,1,1,1,1,1,1,1,1,1,0,0,0.95*t48,48),
                               c(4,1,1,1,1,1,1,0,1,1,0,0,0.05*t48,48),
                               c(4,1,1,1,1,1,1,1,1,1,0,0,t72,72),
                               c(4,1,1,1,1,1,1,1,1,1,0,0,t288,288),
                               c(4,1,1,1,1,1,1,1,1,1,0,0,t384,384),
                               
                               c(4,2,0,0,0,0,0,0,0,0,1,1,10000,0),
                               c(4,2,0,0,0,0,0,0,0,0,1,1,t4,4),
                               c(4,2,0,0,0,0,0,0,0,0,1,1,t48,48),
                               
                               c(4,2,0,0,0,0,1,0,1,0,1,0,0.20*t72,72),
                               c(4,2,0,1,0,0,0,0,1,0,1,0,0.10*t72,72),
                               c(4,2,0,1,0,0,1,0,1,0,1,0,0.05*t72,72),
                               c(4,2,0,0,0,0,0,0,0,0,1,0,0.20*t72,72),
                               c(4,2,0,1,0,0,0,0,0,0,1,0,0.10*t72,72),
                               c(4,2,0,0,0,0,0,0,1,0,1,0,0.35*t72,72),
                               
                               c(4,2,0,1,0,0,1,0,0,0,1,0,0.70*t288,288),
                               c(4,2,0,1,0,0,1,0,1,1,1,0,0.05*t288,288),
                               c(4,2,0,1,0,0,1,0,1,0,1,0,0.05*t288,288),
                               c(4,2,0,0,0,0,0,0,0,0,1,0,0.05*t288,288),
                               c(4,2,0,1,0,0,0,0,0,0,1,0,0.05*t288,288),
                               c(4,2,0,0,0,0,1,0,0,0,1,0,0.10*t288,288),
                               
                               c(4,2,0,1,0,0,1,0,0,0,1,0,0.70*t384,384),
                               c(4,2,0,1,0,0,1,0,1,0,1,0,0.30*t384,384)))

colnames(p1data)<- c("pig","parent_strain","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","freq","time")
colnames(p2data)<- c("pig","parent_strain","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","freq","time")
colnames(p3data)<- c("pig","parent_strain","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","freq","time")
colnames(p4data)<- c("pig","parent_strain","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","freq","time")

pig_data <- rbind(p1data, p2data, p3data, p4data) %>% group_by(pig, parent_strain, time) %>% 
  arrange(time, -freq) %>% mutate(total = sum(freq), prop = freq / total, value = 1, profile = cumsum(value)) %>% select(-value)

pig_data$parent = 0
for(i in 1:dim(pig_data)[1]){
  if(pig_data[i,"parent_strain"] == 1){
    if(all.equal(as.numeric(pig_data[i,mge]), as.numeric(Pinit[,mge])) == "TRUE"){pig_data[i,"parent"] = 1}
  }else{
    if(all.equal(as.numeric(pig_data[i,mge]), as.numeric(Qinit[,mge])) == "TRUE"){pig_data[i,"parent"] = 1}
  }
}

# Checks for unique combinations of the MGE
Xu <- uniquecombs(as.matrix(pig_data[,c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")]))
ind <- attr(Xu,"index") # Assigns them a profile number
pig_data$profile_number = ind # put this into the pig data


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

pigg_plot <- as.data.frame(pig_data %>% pivot_longer(cols = v1:v10)) %>% 
  mutate(label = ifelse(value == 1, name,0)) %>% 
  mutate(mge = as.numeric(substr(name,2,3)),
         plasmid = ifelse(name %in% c("v1","v2","v3","v4","v5","v9"),0,1),
         x_centre = ifelse(plasmid == 0, bigr*cos(pt[mge]) + x_centrebig,pt[mge]),
         y_centre = ifelse(plasmid == 0, bigr*sin(pt[mge]) + x_centrebig,y_centre_plasmid))


# how add a rectangle for parent???
#d=data.frame(x1=rep(4.2,5), x2=rep(6.7,5), y1=rep(6.2,5), y2=rep(6.2,5), label = "parent", profile = 1, time= c(4,48,96,288,384)) 
#geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = label), alpha=0.5) + 

pigg_plotp <- pigg_plot %>% filter(parent_strain == 1)
pigg_plotq <- pigg_plot %>% filter(parent_strain == 2)

mm <- max(pigg_plotp$profile)
pigg_plotp$profile = factor(pigg_plotp$profile, levels=seq(mm,1,-1)) # reorder profiles to match McCarthy Fig2 JUST FOR P

pigg_plotp$perc <- paste0(100 * signif(pigg_plotp$prop,2),"%")
pigg_plotq$perc <- paste0(100 * signif(pigg_plotq$prop,2),"%")

gp <- ggplot(pigg_plotp, aes(x0=x_centre, y0 = y_centre, group = profile)) + geom_circle(aes(r = minir, col= label, fill  = label)) + 
  geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
  geom_text(aes(x = 5, y = 6.5, label = perc)) +
  facet_grid(pig + time ~ profile) + 
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

gq <- ggplot(pigg_plotq, aes(x0=x_centre, y0 = y_centre, group = profile)) + geom_circle(aes(r = minir, col= label, fill  = label)) + 
  geom_circle(aes(x0 = x_centrebig, y0 = y_centrebig, r = bigr)) + 
  geom_text(aes(x = 5, y = 6.5, label = perc)) +
  facet_grid(pig + time ~ profile + parent_strain) + 
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
gp +  gq + plot_layout(guides = 'collect') +
 plot_layout(widths = c(1, 2))
ggsave("plots/data_pig.pdf", width = 45, height = 25)

####### Number of profiles over time
profile_o_time = pig_data %>% select(pig, time, profile_number) %>%
  group_by(pig, parent_strain,time) %>% summarise(nprof = n_distinct(profile_number))

ggplot(profile_o_time, aes(x=time, y= nprof, group = interaction(pig, parent_strain))) + 
  geom_line(aes(col = factor(parent_strain), linetype = factor(pig))) + 
  geom_point(aes(col = factor(parent_strain))) + 
  scale_y_continuous("Number of profiles") + 
  scale_x_continuous("Time (hrs)") + 
  scale_color_discrete("Parent strain") + 
  scale_linetype_discrete("Pig") + 
  facet_wrap(~pig)
ggsave("plots/profile_over_time_all_pigs.pdf")


profile_o_time = pig_data %>% select(pig, time, profile_number) %>%
  group_by(pig, parent_strain,time) %>% summarise(nprof = n_distinct(profile_number)) %>%
  group_by(parent_strain, time) %>% summarise(meanprof = mean(nprof), sdprof = sd(nprof))

ggplot(profile_o_time, aes(x=time, y = meanprof, group = parent_strain)) + 
  geom_ribbon(aes(ymin = meanprof - sdprof, ymax = meanprof + sdprof, fill = factor(parent_strain)), alpha= 0.2) + 
  geom_point(aes(col = factor(parent_strain))) + 
  geom_line(aes(col = factor(parent_strain))) + 
  scale_y_continuous("Number of profiles") + 
  scale_x_continuous("Time (hrs)") + 
  scale_color_discrete("Parent strain") + 
  guides(fill = FALSE)
ggsave("plots/profile_over_time_mean_sd.pdf")

#### Profiles over time
pig_data$label <- getLabels(pig_data[,c("parent_strain","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")])

ggplot(pig_data, aes(x=time, y = freq, group = interaction(pig,parent_strain,label))) + 
  geom_line(aes(linetype = factor(pig), col = factor(label))) + 
  facet_wrap(~parent_strain)
ggsave("plots/labels_over_time.pdf")

###### Elements over time
pigg_elements <- pigg_plot %>% group_by(pig, parent_strain, time, name) %>% mutate(prop_mge = value * prop) %>% summarise(sum_prop = sum(prop_mge))

pigg_elements$name = factor(pigg_elements$name, levels=mge, labels = c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4"))

theme_set(theme_bw(base_size = 8))
ggplot(pigg_elements, aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
  geom_point(aes(col = name),size = 1.5) + 
  facet_wrap(name~parent_strain, ncol = 2) 
ggsave("plots/mge_over_time.pdf")

ggplot(pigg_elements, aes(x=sum_prop, y = time, group = interaction(name, pig))) + 
  geom_line(aes(col = factor(pig), linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
  geom_point(aes(col = factor(pig)),size = 1.5) + 
  facet_wrap(name~parent_strain, ncol = 2) 
ggsave("plots/mge_over_prop.pdf")

ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")), aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
  geom_point(aes(col = name),size = 1.5) + 
  facet_wrap(name~parent_strain, ncol = 2) 
ggsave("plots/mge_over_time_moving.pdf")

write.csv(pigg_elements, "pigg_elements.csv")



##################################################################################### Totals of bug ####################################################################################################
time <- c(4,48,72,288,384)
## Not modelling any death rate - just net growth as constant... so keep final time point flat
total_p1_min <- c(50, 10^5, 10^6, 10^6, 10^6)
total_p1_max <- c(3000, 10^6, 10^8, 10^7, 10^7)
total_p2_min <- c(1000, 10^5, 10^6, 10^6, 10^6)
total_p2_max <- c(3000, 10^6, 10^8, 10^7, 10^7)

totals <- as.data.frame(cbind(time, total_p1_min, total_p1_max, total_p2_min, total_p2_max)) %>% 
  pivot_longer(cols = total_p1_min:total_p2_max)
totals$parent <- substr(totals$name,8,8)
totals$lim <- substr(totals$name,10,13)
totals$limparent <- paste0(totals$parent, totals$lim)
ggplot(totals, aes(x=time, y = value, group = limparent)) +
  geom_point(aes(col = limparent)) + 
  geom_line(aes(col = limparent)) + 
  scale_y_log10()

totals_err <- totals %>% pivot_wider(id_cols = c(time,  value, parent), names_from = lim)
ggplot(totals_err, aes(x=time, group = parent)) + 
  geom_ribbon((aes(ymin = min, ymax = max, fill = factor(parent))),alpha = 0.5) + 
  scale_y_log10() +
  scale_fill_discrete("Parent")
ggsave("plots/totals_area_to_hit.pdf")

write.csv(totals, "totals_bug.csv")

