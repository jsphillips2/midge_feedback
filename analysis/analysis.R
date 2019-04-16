#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lme4)
library(car)
library(AICcmodavg)

# import data
# z-score continuous variables
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
gpp_m <- lmer(met ~ temp_z + par_z + (trial + midge_z + site)^2 + I(midge_z^2) + I(midge_z^3)+
                (1|rack/core), 
              data = gpp_dat)
Anova(gpp_m, type = 3, test.statistic = "F")

  
# random effects
summary(gpp_m)$varcor

# coefficients
coefs <- summary(gpp_m)$coefficients %>% round(3)
coefs

# standardize observations
gpp_dat = gpp_dat %>%
  mutate(temp_cor = coefs[which(rownames(coefs)=="temp_z"), 1]*temp_z,
         met_z = met - temp_cor + mean(temp_cor))

# generate predicted data
gpp_nd <- gpp_dat %>% 
  tidyr::expand(trial, site, midge_z = seq(min(midge_z), max(midge_z), length.out = 100), 
                temp_z = 0, par_z = 0, rack = 1, coried = 1) %>%
  mutate(midge_treat = 2*midge_z*sd(gpp_dat$midge_treat) + mean(gpp_dat$midge_treat),
         midge_density = midge_treat*unique(gpp_dat$midge_density/gpp_dat$midge_treat))
gpp_preds <- predictSE.merMod(gpp_m, newdata = gpp_nd, REForm = NA, print.matrix = T)
gpp_nd$pred <- gpp_preds[,1]
gpp_nd$se <- gpp_preds[,2]


# plot
gpp_nd %>%
  ggplot(aes(midge_density/1000,pred, color = site))+
  facet_wrap(~trial)+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se, fill = site),
              alpha = 0.15, linetype = 0)+
  geom_line(aes(y = pred), size = 0.8)+
  geom_point(data = gpp_dat, aes(y = met_z), size = 2, alpha = 0.7)+
  scale_color_manual(values=site_colors)+
  scale_fill_manual(values=site_colors)+
  scale_y_continuous(expression("Gross Primary Production (g "*O[2]~m^{-2}~h^{-1}*")"))+
  scale_x_continuous(expression("Midge Density ("*1000~m^{-2}*")"),
                     breaks = c(15,50,85))





