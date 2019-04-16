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
adults <- read_csv("data/clean/adults.csv") %>%
  mutate(midge_z = (midge_treat - mean(midge_treat, na.rm=T))/(2*sd(midge_treat, na.rm=T)))

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

# p-values
Anova(gpp_m, type = 3, test.statistic = "F")
Anova(update(gpp_m, .~. - trial:midge_z - trial:site - midge_z:site), type = 3, test.statistic = "F")

  
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




#==========
#========== Adults
#==========

# select final data of adults
adults_0 <- adults  %>%
  mutate(duration = as.numeric(max(sampledate) - min(sampledate)))%>%
  filter(sampledate == max(sampledate)) %>%
  filter(midge_treat >= 50)

adults_m <- glmer(cbind(cumem, midge_treat - cumem) ~ midge_z*site + (1|rack/core),
                  data = adults_0, family = "binomial")

# p-values
Anova(adults_m, type = 3)
Anova(update(adults_m, .~.-midge_z:site))

# random effects
summary(adults_m)$varcor

# coefficients
coefs <- summary(adults_m)$coefficients %>% round(3)
coefs

# generate predicted values
adults_nd <- adults_0 %>%
  tidyr::expand(site, midge_z = seq(min(midge_z)-0.4, max(midge_z), length.out = 100), 
                rack = 1, coried = 1) %>%
  mutate(midge_treat = 2*midge_z*sd(adults_0$midge_treat) + mean(adults_0$midge_treat),
         midge_density = midge_treat*unique(adults_0$midge_density/adults_0$midge_treat))
adults_preds<- predictSE.merMod(adults_m, newdata = adults_nd, REForm = NA, print.matrix = T)
adults_nd$pred <- adults_preds[,1]
adults_nd$se <- adults_preds[,2]

# plot
adults_nd %>%
  ggplot(aes(midge_density/1000, pred, color = site))+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se, fill = site),
              alpha = 0.15, linetype = 0)+
  geom_line(aes(y = pred), size = 0.8)+
  geom_jitter(data = adults_0, aes(y = prop), size = 2, alpha = 0.7, width = 2)+
  scale_color_manual(values=site_colors)+
  scale_fill_manual(values=site_colors)+
  scale_y_continuous("Proportion Emerged")+
  scale_x_continuous(expression("Midge Density ("*1000~m^{-2}*")"))





#==========
#========== Adults ~ GPP
#==========

# calculate average GPP and average predicted GPP through time
adults_new <- adults_0 %>%
  full_join(gpp_dat %>%
              filter(midge_treat > 0)%>%
              group_by(core) %>%
              summarize(met_z = mean(met_z),
                        midge_density = unique(midge_density)) %>%
              mutate(met_per = 24*13*met_z/midge_density,
                     met_per_z = (met_per - mean(met_per))/sd(met_per))) 
# fit_model
adults_gpp <- glmer(cbind(cumem, midge_treat - cumem) ~ met_per_z + (1|rack/core),
                  data = adults_new, family = "binomial")

# p-values
Anova(adults_gpp)

# random effects
summary(adults_gpp)$varcor

# coefficients
coefs <- summary(adults_gpp)$coefficients %>% round(3)
coefs

# generate predicted values
adults_gpp_nd <- adults_new %>%
  tidyr::expand(nesting(met_per, met_per_z))
adults_gpp_preds<- predictSE.merMod(adults_gpp, newdata = adults_gpp_nd, REForm = NA, print.matrix = T)
adults_gpp_nd$pred <- adults_gpp_preds[,1]
adults_gpp_nd$se <- adults_gpp_preds[,2]

# plot
adults_gpp_nd %>%
  ggplot(aes(met_per, pred))+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se),
              alpha = 0.15, linetype = 0)+
  geom_line(aes(y = pred), size = 0.8)+
  geom_point(data = adults_new, aes(y = prop, color = site), size = 2, alpha = 0.7)+
  scale_color_manual(values=site_colors)+
  scale_fill_manual(values=site_colors)+
  scale_y_continuous("Proportion Emerged")+
  scale_x_continuous(expression("g C"*~midge^{-1}*""))





#==========
#========== Feedback
#==========

# calculate predicted gpp for observed levels
gpp_mod <- gpp_dat %>% 
  filter(midge_treat >= 50) %>%
  tidyr::expand(trial, site, midge_z = seq(min(midge_z), max(midge_z)+0.3, length.out = 100), 
                temp_z = 0, par_z = 0, rack = 1, coried = 1) %>%
  mutate(midge_treat = 2*midge_z*sd(gpp_dat$midge_treat) + mean(gpp_dat$midge_treat),
         midge_density = midge_treat*unique(gpp_dat$midge_density/gpp_dat$midge_treat)) %>%
  na.omit()

gpp_mod$pred <- predict(gpp_m, newdata = gpp_mod, re.form=NA) 

gpp_mod_sum <- gpp_mod %>%
  filter(midge_treat > 0) %>%
  group_by(site, midge_treat, midge_density) %>%
  summarize(pred = mean(pred))%>%
  mutate(pred_per = 24*13*pred/midge_density,
         pred_per_z = (pred_per - mean(adults_new$met_per))/sd(adults_new$met_per))

gpp_mod_pred <-predictSE.merMod(adults_gpp, newdata = gpp_mod_sum %>%
                                       mutate(met_per_z = pred_per_z), 
                                     REForm = NA, print.matrix = T)
gpp_mod_sum$pred <- gpp_mod_pred[,1]
gpp_mod_sum$se <- gpp_mod_pred[,2]

gpp_null <- gpp_dat %>% 
  filter(midge_treat == 0) %>%
  tidyr::expand(trial, site, midge_z, temp_z = 0, par_z = 0, rack = 1, coried = 1)

gpp_null$pred <- predict(gpp_m, newdata = gpp_null, re.form=NA) 

gpp_null_sum <- gpp_null %>%
  group_by(site) %>%
  summarize(pred = mean(pred)) %>%
  tidyr::expand(nesting(site,pred), 
                midge_density = {gpp_mod %>% filter(midge_treat > 0)}$midge_density+1500) %>%
  mutate(pred_per = 24*13*pred/midge_density,
         pred_per_z = (pred_per - mean(adults_new$met_per))/sd(adults_new$met_per))


gpp_null_pred <-predictSE.merMod(adults_gpp, newdata = gpp_null_sum %>%
                                  mutate(met_per_z = pred_per_z), 
                                REForm = NA, print.matrix = T)
gpp_null_sum$pred <- gpp_null_pred[,1]
gpp_null_sum$se <- gpp_null_pred[,2]


# plot
gpp_mod_sum %>%
  ggplot(aes(midge_density/1000, pred, color = site))+
  facet_wrap(~site, nrow = 1)+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se,fill = site),
              alpha = 0.15, linetype = 0)+
  geom_line(aes(y = pred), size = 0.8)+
  geom_line(data = gpp_null_sum,  aes(y = pred), size = 0.8, linetype = 2)+
  geom_ribbon(data = gpp_null_sum,aes(ymin = pred - se, ymax = pred + se,fill = site),
              alpha = 0.05, linetype = 0)+
  scale_color_manual(values=site_colors)+
  scale_fill_manual(values=site_colors)+
  scale_y_continuous("Proportion Emerged")+
  scale_x_continuous(expression("Midge Density ("*1000~m^{-2}*")"),
                     breaks = c(15,50,85))



