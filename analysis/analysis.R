#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)
# library(nlme)
library(lme4)
library(car)
# library(AICcmodavg)

# import data
# add z-score continuous variables
met_long <- read_csv("data/clean/met_long.csv") %>%
  mutate(midge_z = (midge_treat - mean(midge_treat, na.rm=T))/(2*sd(midge_treat, na.rm=T)),
         temp_z = (temperature - mean(temperature, na.rm=T))/(2*sd(temperature, na.rm=T)),
         par_z = (par - mean(par, na.rm=T))/(2*sd(par, na.rm=T)),
         date = as.factor(sampledate))

# set theme
theme_set(theme_bw() %+replace% 
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=10),
                  legend.text = element_text(size=10),
                  axis.text=element_text(size=10, color="black"),
                  axis.title.y=element_text(angle = 90 ,margin=margin(0,15,0,0)),
                  axis.title.x=element_text(margin=margin(15,0,0,0))))

# define site colors
site_colors <- c("dodgerblue","firebrick","black")
names(site_colors) <- c("E2","E3","E5")



#==========
#========== GPP
#==========

# gpp_dat
gpp_dat <- met_long %>% filter(type == "gpp")

# fit model
gpp_m <- lmer(met ~ temp_z + par_z + (trial + midge_z + site)^2 + I(midge_z^2) + (1|rack/core), 
              data = gpp_dat)
Anova(gpp_m, type = 3, test.statistic = "F")
  
# random effects
summary(gpp_m)$varcor

# coefficients
summary(gpp_m)$coefficients %>% round(3)

# plot
gpp_dat %>%
  ggplot(aes(midge_density,met, color = site))+
  facet_wrap(~sampledate)+
  geom_point()+
  scale_color_manual(values=site_colors)














